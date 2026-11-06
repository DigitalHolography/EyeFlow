"""Isolated input strategies using EyeFlow's real spectral spatial reducer.

Run with PYTHONPATH=src. Sources are read-only; every trial uses fresh processes
and private local scratch. Results include reduced-output HDF5 publication, NOT
the complete EyeFlow DAG/artifact workflow. No OS/server caches are purged.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from contextlib import ExitStack
import json
import multiprocessing as mp
from multiprocessing.shared_memory import SharedMemory
from pathlib import Path
import platform
import random
import shutil
import tempfile
from time import perf_counter
import traceback

import h5py
import numpy as np
import psutil

from pipeline_engine.data.storage import ProductStorage
from scratch_location import is_local_path
from pipelines.low_freq_spectra.spatial_reduction import SpatialReducer, spectral_batches

MODES = ("npy", "npy-direct", "ram", "ram-direct", "thread-full", "process-full",
         "slab-sync", "slab-thread", "local-h5")


class MeasuredSource:
    def __init__(self, source):
        self.source = source
        self.spans = []
        self.logical_bytes = 0

    def __getattr__(self, name):
        return getattr(self.source, name)

    def __getitem__(self, selection):
        start = perf_counter()
        value = self.source[selection]
        self.spans.append((start, perf_counter()))
        self.logical_bytes += value.nbytes
        return value

    def read_direct(self, destination, source_sel=None, dest_sel=None):
        start = perf_counter()
        self.source.read_direct(destination, source_sel=source_sel, dest_sel=dest_sel)
        self.spans.append((start, perf_counter()))
        self.logical_bytes += destination[dest_sel if dest_sel is not None else ()].nbytes


def describe(path):
    """Metadata only, including filters/cache; never materialize source values."""
    datasets = {}
    with h5py.File(path, "r") as source:
        def visit(name, node):
            if not isinstance(node, h5py.Dataset):
                return
            properties = node.id.get_create_plist()
            datasets[name] = dict(
                shape=node.shape, dtype=str(node.dtype), decoded_bytes=node.nbytes,
                stored_bytes=node.id.get_storage_size(), chunks=node.chunks,
                filters=[repr(properties.get_filter(i))
                         for i in range(properties.get_nfilters())],
                chunk_cache=node.id.get_access_plist().get_chunk_cache(),
            )
        source.visititems(visit)
    return dict(path=str(path), file_bytes=Path(path).stat().st_size, datasets=datasets)


def shared_reader(path, dataset, name, shape, dtype, connection, cache):
    """One independently opened spawn reader; sends only small control data."""
    memory = None
    array = None
    try:
        connection.send(("started", perf_counter()))
        memory = SharedMemory(name=name)
        array = np.ndarray(shape, dtype=np.dtype(dtype), buffer=memory.buf)
        with h5py.File(path, "r", **cache) as source:
            measured = MeasuredSource(source[dataset])
            measured.read_direct(array)
            connection.send(("ready", measured.spans, measured.logical_bytes))
    except BaseException:
        connection.send(("error", traceback.format_exc()))
    finally:
        array = None
        if memory is not None:
            memory.close()
        connection.close()


def stage_direct(storage, source):
    """Experimental single HDF5 read into an existing local replay mapping."""
    from pipeline_engine.data.storage import ArrayReference
    if shutil.disk_usage(storage.directory).free < source.nbytes + 4096:
        raise OSError("Insufficient local scratch space")
    staged = storage.path(".npy")
    output = None
    try:
        output = np.lib.format.open_memmap(staged, mode="w+", dtype=source.dtype, shape=source.shape)
        source.read_direct(output)
        output.flush()
    except BaseException:
        if output is not None:
            output._mmap.close()
        staged.unlink(missing_ok=True)
        raise
    else:
        output._mmap.close()
    return ArrayReference(str(staged), source.shape, source.dtype)


def prefetch_batches(source, buffer_bytes):
    """One owned block being consumed plus at most one pending block.

    spectral_batches reuses its buffer, so advancing it requires copying before
    handing off. Use independently owned slice results instead: no recycling and
    no mutation of a published block. The producer is joined on all exit paths.
    """
    frame_bytes = int(np.prod(source.shape[1:])) * source.dtype.itemsize
    count = min(source.shape[0], max(1, buffer_bytes // frame_bytes))
    if source.chunks and source.chunks[0] <= count:
        count = count // source.chunks[0] * source.chunks[0]
    bounds = iter(range(0, source.shape[0], count))
    with ThreadPoolExecutor(max_workers=1) as executor:
        start = next(bounds)
        pending = executor.submit(source.__getitem__, np.s_[start:start + count])
        while True:
            block = pending.result()
            following = next(bounds, None)
            pending = (executor.submit(source.__getitem__, np.s_[following:following + count])
                       if following is not None else None)
            yield start, block
            if pending is None:
                break
            start = following


def trial(config, directory):
    mode = config["mode"]
    started = perf_counter()
    metrics = dict(mode=mode, source=config["source"], buffer_mib=config["buffer_mib"],
                   network_bytes=None, physical_disk_bytes=None, startup_s=0.0,
                   scratch_bytes=0, file_copy_s=0.0)
    cache = {k: config[k] for k in ("rdcc_nbytes", "rdcc_nslots") if config.get(k) is not None}
    with ExitStack() as stack:
        path = config["source"]
        if mode == "local-h5":
            copied = Path(directory) / "source.h5"
            before = perf_counter()
            shutil.copyfile(path, copied)
            metrics.update(file_copy_s=perf_counter() - before, scratch_bytes=copied.stat().st_size)
            path = str(copied)
        source = MeasuredSource(stack.enter_context(h5py.File(path, "r", **cache))[config["dataset"]])
        if source.ndim != 4 or source.dtype.kind not in "iuf" or min(source.shape) < 1:
            raise ValueError("Select a nonempty real numeric (time,freq,y,x) dataset.")
        shape, dtype = source.shape, source.dtype
        if mode in {"ram", "ram-direct", "thread-full", "process-full"}:
            if source.nbytes > config["max_ram_mib"] * 2**20:
                raise MemoryError("Full-array candidate exceeds --max-ram-mib; no silent fallback.")
        wait_spans = []
        if mode in {"npy", "npy-direct"}:
            storage = stack.enter_context(ProductStorage(directory=Path(directory) / "arrays"))
            before = perf_counter()
            if mode == "npy":
                reference = storage.store_array(source)
            else:
                reference = stage_direct(storage, source)
            data = storage.read(reference)
            wait_spans.append((before, perf_counter()))
            metrics["scratch_bytes"] = Path(reference.path).stat().st_size
        elif mode in {"ram", "ram-direct"}:
            before = perf_counter()
            if mode == "ram":
                data = source[()]
            else:
                data = np.empty(shape, dtype=dtype)
                source.read_direct(data)
            wait_spans.append((before, perf_counter()))
        elif mode == "thread-full":
            executor = stack.enter_context(ThreadPoolExecutor(max_workers=1))
            before = perf_counter()
            future = executor.submit(source.__getitem__, ())
            data = future.result()
            wait_spans.append((before, perf_counter()))
        elif mode == "process-full":
            memory = SharedMemory(create=True, size=source.nbytes)
            # Cleanup order: terminate/join reader, drop ndarray, close/unlink owner.
            try:
                parent, child = mp.get_context("spawn").Pipe(duplex=False)
                process = mp.get_context("spawn").Process(
                    target=shared_reader,
                    args=(path, config["dataset"], memory.name, shape, dtype.str, child, cache))
                before = perf_counter()
                process.start()
                child.close()
                try:
                    if not parent.poll(config["timeout"]):
                        raise TimeoutError("Reader startup timed out")
                    message = parent.recv()
                    if message[0] != "started":
                        raise RuntimeError(message)
                    metrics["startup_s"] = message[1] - before
                    if not parent.poll(config["timeout"]):
                        raise TimeoutError("Reader timed out")
                    message = parent.recv()
                    if message[0] != "ready":
                        raise RuntimeError(message)
                    source.spans, source.logical_bytes = message[1:]
                    wait_spans.append((before, perf_counter()))
                    data = np.ndarray(shape, dtype=dtype, buffer=memory.buf)
                    data.flags.writeable = False
                    result = consume(data, config, started, wait_spans)
                finally:
                    if process.is_alive():
                        process.join(5)
                    if process.is_alive():
                        process.terminate()
                    process.join()
                    parent.close()
                    data = None
            finally:
                memory.close()
                memory.unlink()
        else:
            data = source
        if mode != "process-full":
            result = consume(data, config, started, wait_spans)
        outputs, compute_spans, first_block_s = result
        before = perf_counter()
        output_path = Path(directory) / "reduced.h5"
        with h5py.File(output_path, "x") as output:
            for index, values in enumerate(outputs):
                output.create_dataset(str(index), data=values, compression="gzip", compression_opts=9)
        metrics.update(
            publication_s=perf_counter() - before,
            total_through_publication_s=perf_counter() - started,
            first_block_s=first_block_s,
            hdf5_read_s=sum(b-a for a, b in source.spans),
            source_calls=len(source.spans), logical_bytes=source.logical_bytes,
            input_wait_s=sum(b-a for a, b in wait_spans),
            compute_s=sum(b-a for a, b in compute_spans),
            read_compute_overlap_s=sum(max(0, min(b,d)-max(a,c))
                                       for a,b in source.spans for c,d in compute_spans),
        )
        np.savez(Path(directory) / "values.npz", *outputs)
    return metrics


def consume(data, config, started, wait_spans):
    buffer_bytes = config["buffer_mib"] * 2**20
    batches = (prefetch_batches(data, buffer_bytes) if config["mode"] == "slab-thread"
               else spectral_batches(data, read_buffer_bytes=buffer_bytes))
    outputs = [np.empty(data.shape[:2], dtype=np.float64) for _ in range(3)]
    # Fixed, nontrivial gap support. It exercises the production gathered median.
    gaps = np.flatnonzero(np.indices(data.shape[-2:]).sum(axis=0).ravel() % 3 == 0)
    spans = []
    first = None
    with SpatialReducer(workers=config["workers"]) as reducer:
        iterator = iter(batches)
        try:
            while True:
                before = perf_counter()
                try:
                    start, block = next(iterator)
                except StopIteration:
                    break
                wait_spans.append((before, perf_counter()))
                if first is None:
                    first = perf_counter() - started
                before = perf_counter()
                values = reducer(block, gaps)
                for output, value in zip(outputs, values):
                    output[start:start + len(block)] = value
                spans.append((before, perf_counter()))
        finally:
            iterator.close()
    for output in outputs:
        np.maximum(output, 0, out=output)
    return outputs, spans, first


def trial_entry(config, directory, result_path):
    try:
        result = trial(config, directory)
    except BaseException:
        result = {"error": traceback.format_exc(), "mode": config["mode"]}
    Path(result_path).write_text(json.dumps(result, indent=2), encoding="utf-8")


def run_trial(config, root):
    """Sample trial + reader RSS from the supervisor, outside the read path."""
    with tempfile.TemporaryDirectory(prefix="trial-", dir=root) as directory:
        result_path = Path(directory) / "result.json"
        process = mp.get_context("spawn").Process(target=trial_entry, args=(config, directory, str(result_path)))
        before = perf_counter()
        process.start()
        peak = 0
        try:
            while process.is_alive():
                try:
                    owner = psutil.Process(process.pid)
                    peak = max(peak, sum(p.memory_info().rss for p in [owner, *owner.children(recursive=True)]))
                except psutil.Error:
                    pass
                if perf_counter() - before > config["timeout"]:
                    # Includes native-read hangs. Stop descendants before allocator.
                    for child in psutil.Process(process.pid).children(recursive=True):
                        child.kill()
                    process.terminate()
                    raise TimeoutError(f"Trial exceeded {config['timeout']} seconds")
                process.join(0.02)
        finally:
            if process.is_alive():
                process.terminate()
            process.join()
        if not result_path.exists():
            raise RuntimeError(f"Trial exited {process.exitcode} without a result")
        result = json.loads(result_path.read_text(encoding="utf-8"))
        if "error" in result:
            raise RuntimeError(result["error"])
        result["sampled_process_tree_rss_bytes"] = peak
        result["supervised_wall_s"] = perf_counter() - before
        reference = Path(root) / "reference.npz"
        values_path = Path(directory) / "values.npz"
        if reference.exists():
            with np.load(reference) as expected, np.load(values_path) as values:
                for name in expected.files:
                    np.testing.assert_array_equal(values[name], expected[name], strict=True)
        else:
            shutil.copyfile(values_path, reference)
        result["exact_output_equivalence"] = True
        return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--dataset", default="spectrograms")
    parser.add_argument("--output", type=Path, required=True, help="New local result directory")
    parser.add_argument("--modes", nargs="+", choices=MODES, default=list(MODES))
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--buffer-mib", type=int, default=64)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--max-ram-mib", type=int, default=4096)
    parser.add_argument("--rdcc-nbytes", type=int)
    parser.add_argument("--rdcc-nslots", type=int)
    parser.add_argument("--timeout", type=float, default=600)
    parser.add_argument("--cache-label", default="uncontrolled; no cache purge")
    args = parser.parse_args()
    if min(args.repeats, args.buffer_mib, args.workers, args.max_ram_mib, args.timeout) <= 0:
        parser.error("Counts, limits and timeout must be positive")
    if not is_local_path(args.output):
        parser.error("--output must be on local storage")
    args.output.mkdir(parents=True, exist_ok=False)
    metadata = describe(args.source)
    import subprocess
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    report = dict(revision=revision, metadata=metadata, environment=dict(python=platform.python_version(),
                  h5py=h5py.__version__, hdf5=h5py.version.hdf5_version,
                  numpy=np.__version__, platform=platform.platform()),
                  cache_label=args.cache_label, trials=[])
    config = vars(args).copy()
    config["source"] = str(args.source.resolve())
    config.pop("output")
    report["configuration"] = config
    plan = [(repeat, mode) for repeat in range(args.repeats) for mode in args.modes]
    random.Random(2026).shuffle(plan)
    (args.output / "results.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    for repeat, mode in plan:
        result = run_trial({**config, "mode": mode}, args.output)
        result["repeat"] = repeat
        report["trials"].append(result)
        (args.output / "results.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
        print(json.dumps(result), flush=True)


if __name__ == "__main__":
    main()
