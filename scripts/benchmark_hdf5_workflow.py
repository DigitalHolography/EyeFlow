"""Run the real DAG through final publication with an isolated staging prototype.

No installed settings are changed. Run each trial in a fresh interpreter.
Use --settings to pin the identical exported JSON configuration for SSD/NAS.
Only --strategy npy-direct patches HDF5-to-replay staging for this process;
all ndarray/product storage and scientific consumers retain production behavior.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import json
import subprocess
import shutil
from pathlib import Path
from itertools import product
from time import perf_counter
from unittest.mock import patch

import h5py
import numpy as np
import psutil

from benchmark_hdf5_input_io import stage_direct
from app_settings import AppSettingsStore
from runtime_limits import configure_numeric_threads
from pipeline_engine import execute_run, resolve_run_spec, selectable_pipeline_registry
from pipeline_engine.data.storage import ProductStorage
from scratch_location import is_local_path
from pipeline_engine.output.output_model import OutputModel
from pipelines import load_pipeline_catalog
from utils.logger import Logger
from input_output import resolve_selected_run_layouts


def _dataset_blocks(dataset, target_bytes=16 * 1024**2):
    if not dataset.shape or not dataset.size:
        yield ()
        return
    block = list(dataset.shape)
    for axis, size in enumerate(block):
        other_bytes = int(np.prod(block)) // size * dataset.dtype.itemsize
        # Keep chunk boundaries so successive blocks do not decompress the same
        # profile chunk repeatedly when its leading axis spans a whole segment.
        chunk = dataset.chunks[axis] if dataset.chunks else 1
        count = max(chunk, (target_bytes // other_bytes) // chunk * chunk)
        block[axis] = min(size, count)
        if int(np.prod(block)) * dataset.dtype.itemsize <= target_bytes:
            break
    for starts in product(*(range(0, size, step) for size, step in zip(dataset.shape, block))):
        yield tuple(slice(start, min(start + step, size))
                    for start, step, size in zip(starts, block, dataset.shape))


def compare_hdf5(first, second, *, allowed_root_differences=()):
    """Check every value, scientific attribute and node, allowing source locations."""
    ignored = {"source_file", "holodoppler_source_file", "doppler_vision_source_file",
               "input_read_direct"}
    ignored.update(allowed_root_differences)
    checked = dict(datasets=0, numeric_datasets=0, attributes=0, numeric_bytes=0)
    with h5py.File(first, "r") as left, h5py.File(second, "r") as right:
        root_differences = []
        for key in set(left.attrs) | set(right.attrs):
            try:
                assert key in left.attrs and key in right.attrs
                np.testing.assert_equal(left.attrs[key], right.attrs[key])
            except AssertionError:
                root_differences.append(key)
        nodes = [""]
        other = [""]
        left.visit(nodes.append)
        right.visit(other.append)
        if nodes != other:
            raise AssertionError("Output node paths differ")
        for path in nodes:
            a, b = (left[path], right[path]) if path else (left, right)
            if isinstance(a, h5py.Dataset) != isinstance(b, h5py.Dataset):
                raise AssertionError(f"Node type differs: {path}")
            keys = set(a.attrs) - (ignored if not path else set())
            if keys != set(b.attrs) - (ignored if not path else set()):
                raise AssertionError(f"Attribute keys differ: {path}")
            for key in keys:
                np.testing.assert_equal(a.attrs[key], b.attrs[key], err_msg=f"{path}:{key}")
                checked["attributes"] += 1
            if isinstance(a, h5py.Dataset):
                if a.shape != b.shape or a.dtype != b.dtype:
                    raise AssertionError(f"Shape/dtype differs: {path}")
                for block in _dataset_blocks(a):
                    av, bv = a[block], b[block]
                    np.testing.assert_array_equal(av, bv, strict=True, err_msg=path)
                    if a.dtype.kind in "biufc" and np.asarray(av).tobytes() != np.asarray(bv).tobytes():
                        raise AssertionError(f"Numeric bits differ (including NaNs/zero signs): {path}")
                if a.dtype.kind in "biufc":
                    checked["numeric_datasets"] += 1
                    checked["numeric_bytes"] += a.nbytes
                checked["datasets"] += 1
    return dict(**checked, equivalent=True, ignored_root_attributes=sorted(ignored),
                root_attribute_differences=sorted(root_differences))


@contextmanager
def staging_strategy(strategy, counters):
    original = ProductStorage.store_array

    def store(storage, source):
        if not isinstance(source, h5py.Dataset):
            return original(storage, source)
        process = psutil.Process()
        io_before = process.io_counters()
        net_before = psutil.net_io_counters().bytes_recv
        before = perf_counter()
        result = (stage_direct(storage, source) if strategy == "npy-direct"
                  and source.ndim and source.dtype.kind in "biufc" and source.size
                  else original(storage, source))
        counters.append(dict(path=source.name, decoded_bytes=source.nbytes,
                             staging_s=perf_counter()-before,
                             process_read_api_bytes=process.io_counters().read_bytes-io_before.read_bytes,
                             host_network_receive_bytes=psutil.net_io_counters().bytes_recv-net_before))
        return result

    with patch.object(ProductStorage, "store_array", store):
        yield


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--holo", type=Path)
    parser.add_argument("--output", type=Path, required=True, help="New local trial directory")
    parser.add_argument("--settings", type=Path)
    parser.add_argument("--compare", type=Path, nargs=2, metavar=("BASELINE_H5", "CANDIDATE_H5"))
    parser.add_argument("--allow-runtime-differences", action="store_true",
                        help="Allow stage layout/product diagnostics and timing metadata to differ.")
    parser.add_argument("--strategy", choices=("npy", "npy-direct", "local-h5"), default="npy")
    args = parser.parse_args()
    if not is_local_path(args.output):
        parser.error("Use local benchmark output to isolate source I/O")
    args.output.mkdir(parents=True, exist_ok=False)
    if args.compare:
        result = compare_hdf5(*args.compare, allowed_root_differences=(
            ("stage_status", "stage_products", "stage_timings", "artifact_timings")
            if args.allow_runtime_differences else ()
        ))
        (args.output / "equivalence.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
        print(json.dumps(result, indent=2))
        return
    if args.holo is None or args.settings is None:
        parser.error("--holo and --settings are required for workflow runs")
    configure_numeric_threads()
    settings = AppSettingsStore(path=args.settings)
    # Snapshot the effective settings; never normalize/save the user's active file.
    raw = settings.load()
    (args.output / "settings.json").write_text(json.dumps(raw, indent=2), encoding="utf-8")
    available, missing = load_pipeline_catalog()
    catalog = (*available, *missing)
    visible = selectable_pipeline_registry(catalog)
    targets = tuple(name for name in visible if raw.get("pipeline_visibility", {}).get(name, False))
    options = {name: tuple(key for key, enabled in values.items() if enabled)
               for name, values in raw.get("pipeline_options", {}).items()}
    before = perf_counter()
    network_before = psutil.net_io_counters().bytes_recv
    input_path = args.holo
    copied_bytes = 0
    copy_seconds = 0.0
    if args.strategy == "local-h5":
        layout, = resolve_selected_run_layouts([args.holo])
        input_path = args.output / "inputs" / layout.stem
        copy_started = perf_counter()
        for source in (layout.hd_h5, layout.dv_h5):
            destination = input_path / source.parent.parent.name
            (destination / "h5").mkdir(parents=True)
            shutil.copyfile(source, destination / "h5" / source.name)
            sidecars = source.parent.parent / "json"
            if sidecars.is_dir():
                shutil.copytree(sidecars, destination / "json")
        copied_bytes = sum(path.stat().st_size for path in input_path.rglob("*") if path.is_file())
        copy_seconds = perf_counter()-copy_started
    spec = resolve_run_spec(
        input_paths=[input_path], target_names=targets, pipelines=catalog, pipeline_options=options,
        band_ratio_frequency_scale_hz=settings.load_band_ratio_frequency_scale_hz(),
        band_limited_harmonics=settings.load_band_limited_harmonics(),
        copy_hd_passthrough=settings.load_copy_hd_passthrough(),
        velocity_estimation_method=settings.load_velocity_estimation_method(),
        output_root=args.output / "results",
    )
    staging = []
    phases = []
    publication = []
    resources = []
    publish = OutputModel.publish

    def measured_publish(model):
        start = perf_counter()
        try:
            return publish(model)
        finally:
            publication.append(perf_counter()-start)

    def capture_resource(message):
        if "[RESOURCES]" in message:
            resources.append(message)

    with (staging_strategy(args.strategy, staging),
          patch.object(OutputModel, "publish", measured_publish),
          Logger.capture_messages(capture_resource)):
        result = execute_run(spec, on_pipeline_start=lambda name, *_: phases.append(
            dict(pipeline=name, start_s=perf_counter()-before)))
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    report = dict(revision=revision, strategy=args.strategy, holo=str(args.holo), targets=spec.plan.targets,
                  total_through_publication_s=perf_counter()-before,
                  source_staging=staging, pipeline_starts=phases,
                  final_hdf5_publication_s=sum(publication),
                  file_copy_s=copy_seconds, copied_input_bytes=copied_bytes,
                  host_network_receive_bytes=psutil.net_io_counters().bytes_recv-network_before,
                  resource_messages=resources,
                  outputs=[str(path) for path in result.outputs],
                  failures=[str(failure) for failure in result.failures])
    (args.output / "timing.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2), flush=True)
    if result.failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
