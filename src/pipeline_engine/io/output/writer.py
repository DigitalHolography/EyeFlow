"""Exclusive HDF5 writer process; it never accepts scientific transactions."""

import hashlib
import pickle
import shutil
from math import prod
from pathlib import Path
from time import perf_counter, process_time

import h5py
import numpy as np
import psutil

from input_output.writers.h5_publication import write_record, reconcile_attributes
from runtime_limits import configure_native_threads
from scratch_location import is_local_path

from .protocol import WriteAcknowledgement


def _required_free_bytes(record):
    """Reserve a whole next write, retaining native sparse HDF5 savings.

    Dense input files may compress poorly, so checking just the next buffer is
    insufficient. Native HDF5 copying preserves allocated chunks and filters;
    its physical backing is the appropriate bound instead of its logical shape.
    """
    from pipeline_engine.data.storage import H5ArrayReference

    value = record.value
    physical = Path(value.path).stat().st_size if hasattr(value, "path") else getattr(value, "nbytes", 0)
    amount = physical
    if not isinstance(value, H5ArrayReference):
        shape, dtype = getattr(value, "shape", None), getattr(value, "dtype", None)
        if shape is not None and dtype is not None:
            itemsize = max(np.dtype(dtype).itemsize, np.dtype(record.options.get("dtype", dtype)).itemsize)
            amount = max(amount, prod(shape) * itemsize)
            chunks = record.options.get("chunks")
            if isinstance(chunks, (tuple, list)) and len(chunks) == len(shape) and all(chunks):
                count = prod((size + chunk - 1) // chunk for size, chunk in zip(shape, chunks))
                # Edge chunks are allocated at full chunk size; allow index and
                # filter overhead even when explicit chunks are unusually small.
                amount = max(amount, count * (prod(chunks) * itemsize + 128))
    return amount + amount // 20 + 16 * 2**20


def _process_sample(process):
    """OS process counters include cached I/O; they are not device telemetry."""
    sample = {"rss_bytes": process.memory_info().rss}
    try:
        counters = process.io_counters()
    except (AttributeError, OSError, psutil.Error):
        return sample
    for name in ("read_bytes", "write_bytes", "read_count", "write_count"):
        sample[name] = getattr(counters, name)
    return sample


def _load_manifest(command):
    path = Path(command.manifest)
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(2**20):
            digest.update(block)
        if digest.hexdigest() != command.digest:
            raise ValueError("Output command manifest checksum differs")
        stream.seek(0)
        return pickle.load(stream)


def writer_main(commands, acknowledgements, local_path, cache_bytes):
    """Open once, execute accepted commands in order, close before final ACK."""
    configure_native_threads(1)
    local = Path(local_path)
    output = None
    sequence = 0
    completed = {}
    datasets = {}
    try:
        if not is_local_path(local.parent):
            raise OSError(f"EyeFlow output scratch must be local: {local.parent}")
        local.parent.mkdir(parents=True, exist_ok=True)
        process = psutil.Process()
        output = h5py.File(local, "w", rdcc_nbytes=cache_bytes)
        acknowledgements.send(WriteAcknowledgement(0))
        idle_at = perf_counter()
        while True:
            command = commands.recv()
            started = perf_counter()
            cpu_started = process_time()
            before = _process_sample(process)
            metrics = {'output_writer_idle_seconds':started-idle_at}
            messages = []
            if command.accepted_at:
                metrics['output_commit_to_write_start_seconds'] = started-command.accepted_at
            try:
                previous = completed.get(command.sequence)
                if previous is not None:
                    if previous != (command.digest, command.final):
                        raise ValueError("Conflicting repeated output sequence")
                    acknowledgements.send(WriteAcknowledgement(command.sequence))
                    continue
                if command.sequence != sequence + 1:
                    raise ValueError("Out-of-order output sequence")
                contents = _load_manifest(command)
                if command.final:
                    groups, expected = contents
                    if set(expected) != set(datasets):
                        raise ValueError("Accepted output dataset inventory differs from writer inventory")
                    if any(path not in output or not isinstance(output[path], h5py.Dataset) for path in expected):
                        raise ValueError("A serialized scientific dataset is missing")
                    for path, record in sorted(groups.items(), key=lambda item: (item[0].count('/'), item[0])):
                        group = output.require_group(path) if path else output
                        reconcile_attributes(group, record)
                    flushed = perf_counter()
                    output.flush()
                    metrics["output_flush_seconds"] = perf_counter() - flushed
                    closing = perf_counter()
                    output.close()
                    metrics['output_close_seconds'] = perf_counter()-closing
                    output = None
                    metrics["output_completed_file_bytes"] = local.stat().st_size
                    metrics["output_final_metadata_seconds"] = perf_counter() - started
                else:
                    for path, record in contents:
                        if not path or path != path.strip('/') or '\\' in path or record.kind != "dataset" or not record.committed:
                            raise ValueError("Writer received an invalid or uncommitted scientific dataset")
                        signature = hashlib.sha256(pickle.dumps(record, protocol=5)).hexdigest()
                        if path in datasets:
                            if datasets[path] != signature:
                                raise ValueError(f"Conflicting accepted output dataset '{path}'")
                            continue
                        value = record.value
                        physical = Path(value.path).stat().st_size if hasattr(value, "path") else getattr(value, "nbytes", 0)
                        if shutil.disk_usage(local.parent).free < _required_free_bytes(record):
                            raise OSError("Insufficient local scratch space for the next output dataset")
                        written = perf_counter()
                        node = write_record(output, path, record, metrics=metrics)
                        duration = perf_counter() - written
                        metrics["output_dataset_write_seconds"] = metrics.get("output_dataset_write_seconds", 0) + duration
                        metrics["output_datasets_written"] = metrics.get("output_datasets_written", 0) + 1
                        metrics["output_source_backing_bytes"] = metrics.get("output_source_backing_bytes", 0) + physical
                        datasets[path] = signature
                        stored = node.id.get_storage_size()
                        metrics['output_dataset_storage_bytes'] = metrics.get('output_dataset_storage_bytes',0)+stored
                        messages.append(f"output_dataset path={path} shape={node.shape} dtype={node.dtype} logical_bytes={node.size*node.dtype.itemsize} storage_bytes={stored} accepted_at={command.accepted_at:.6f} at={written:.6f} end={written+duration:.6f} seconds={duration:.6f} options={dict(record.options)}")
                    flushed = perf_counter()
                    output.flush()
                    metrics['output_progress_flush_seconds'] = perf_counter()-flushed
                sequence = command.sequence
                completed[sequence] = command.digest, command.final
                metrics["output_writer_operation_seconds"] = perf_counter() - started
                metrics["output_writer_cpu_seconds"] = process_time() - cpu_started
                after = _process_sample(process)
                metrics["output_writer_sampled_rss_bytes"] = max(before["rss_bytes"], after["rss_bytes"])
                for name in ("read_bytes", "write_bytes", "read_count", "write_count"):
                    if name in before and name in after:
                        metrics[f"output_writer_os_{name}"] = max(0, after[name] - before[name])
                written_paths = () if command.final else tuple(path for path, _ in contents)
                acknowledgements.send(WriteAcknowledgement(sequence, metrics=metrics, paths=written_paths, messages=tuple(messages)))
                idle_at = perf_counter()
                if command.final:
                    return
            except BaseException as exc:
                acknowledgements.send(WriteAcknowledgement(command.sequence, f"{type(exc).__name__}: {exc}"))
                return
    except BaseException as exc:
        acknowledgements.send(WriteAcknowledgement(sequence + 1, f"{type(exc).__name__}: {exc}"))
    finally:
        if output is not None:
            output.close()
