"""Synchronized host-to-host timings on synthetic cross-section movies.

Run from the repository with CuPy installed:
    python benchmarks/benchmark_cross_section_gpu.py --output benchmark.json
Optionally pass --baseline-file containing the original main implementation.
"""

from __future__ import annotations

import argparse
import gc
import importlib
import importlib.util
import json
import os
import statistics
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import cupy as cp
import numpy as np
import scipy

from calculations.compute_backend import optional_cupy_backend
from calculations.topology import (SegmentTopology, PreparedTopology, SegmentRingSettings,
    PreparedSegment, prepare_segments, interpolate_segment_masks, rotate_segment_masks,
    resample_rotate_segment)

cs = importlib.import_module(
    "calculations.blood_flow_velocity.cross_section.generate_cross_section_signals"
)


def movie(frames, side, varying):
    y, x = np.mgrid[:side, :side].astype(np.float32)
    center = (side - 1) / 2
    mask = (np.abs(x - center) <= side / 6) & (np.abs(y - center) <= side / 3)
    spatial = 2 + 10 * np.maximum(0, 1 - ((x - center) / (side / 5)) ** 2)
    temporal = 1 + 0.3 * np.sin(np.arange(frames, dtype=np.float32) * 2 * np.pi / 32)
    values = temporal[:, None, None] * spatial[None]
    values[:, :2] = np.nan
    if varying:
        values[::13, side // 2, side // 2] = np.nan
    return values.astype(np.float32), mask


def select_backend(mode):
    os.environ["EYEFLOW_COMPUTE_BACKEND"] = mode
    optional_cupy_backend.cache_clear()


def topology_measurement(values, mask, *, windowed=False):
    side = values.shape[-1]
    angle = -59.0
    geometry = SegmentTopology(
        spatial_shape=mask.shape, optic_disc_center_xy=(0., 0.),
        labels=mask.astype(np.int32), centerline=mask,
        branch_ids=np.array([1], np.int32), annulus_masks=np.ones((1, *mask.shape), bool),
        segment_masks=mask[None, None], segment_centers_xy=np.array([[[(side - 1) / 2] * 2]]),
        window_bounds_xyxy=np.array([[[0, side, 0, side]]]), window_side_pixels=side,
        branch_identity=type("Branches", (), {
            "branch_ids": np.array([1], np.int32), "labels": mask.astype(np.int32),
            "stages": type("Stages", (), {"vessel": mask})(),
        })(),
    )
    angles = np.array([[angle]], np.float32)
    resized = interpolate_segment_masks(mask[None, None])
    prepared = PreparedTopology(geometry, angles, resized, rotate_segment_masks(resized, angles))
    device = optional_cupy_backend() is not None
    segments = prepare_segments(values, prepared, keep_on_device=device) if windowed else iter([
        PreparedSegment(0, 0, resample_rotate_segment(values, angle, return_device=device))
    ])
    return cs._generate_cross_section_signals_from_prepared(
        values, prepared, segments, SegmentRingSettings(0., 1., 1., 1),
        cs.CrossSectionSignalSettings(.01), retain_velocity_maps=False,
    )


def measure(module, mode, values, mask, repeats, windowed=False):
    select_backend(mode)
    cp.get_default_memory_pool().free_all_blocks()

    def call():
        if module is cs:
            return topology_measurement(values, mask, windowed=windowed)
        return module._cross_section_velocity_from_substack(
            values,
            mask,
            (values.shape[-1] // 2,) * 2,
            (0, 0),
            -59.0,
            module.CrossSectionSignalSettings(True, 0.5, True, 0.01),
            values.shape[-1],
        )

    # Warm every shape and mode, including compilation and allocator initialization.
    warm = call()
    cp.cuda.Stream.null.synchronize()
    del warm
    times = []
    for _ in range(repeats):
        gc.collect()
        cp.cuda.Stream.null.synchronize()
        start = time.perf_counter()
        result = call()
        cp.cuda.Stream.null.synchronize()
        times.append(time.perf_counter() - start)
        del result
    return {
        "seconds": times,
        "median_seconds": statistics.median(times),
        "cuda_pool_reserved_mib": cp.get_default_memory_pool().total_bytes() / 1024**2,
    }


def parity(values, mask):
    results = []
    for mode in ("cpu", "cupy"):
        select_backend(mode)
        results.append(topology_measurement(values, mask))
    errors = {}
    for name in ("velocity", "safe_velocity", "velocity_profiles",
                 "transverse_velocity_profiles_masked",
                 "longitudinal_velocity_profiles_unmasked",
                 "longitudinal_velocity_profiles_masked"):
        cpu, gpu = (getattr(result, name) for result in results)
        np.testing.assert_allclose(gpu, cpu, rtol=1e-5, atol=1e-5, equal_nan=True)
        errors[name] = float(np.nanmax(np.abs(cpu - gpu)))
    return errors


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--baseline-file", type=Path)
    parser.add_argument("--repeats", type=int, default=3)
    args = parser.parse_args()
    baseline = None
    if args.baseline_file:
        name = "calculations.blood_flow_velocity.cross_section._benchmark_baseline"
        spec = importlib.util.spec_from_file_location(name, args.baseline_file)
        baseline = importlib.util.module_from_spec(spec)
        sys.modules[name] = baseline
        spec.loader.exec_module(baseline)
    props = cp.cuda.runtime.getDeviceProperties(0)
    report = {
        "gpu": props["name"].decode(),
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "cupy": cp.__version__,
        "cuda_runtime": cp.cuda.runtime.runtimeGetVersion(),
        "cuda_driver": cp.cuda.runtime.driverGetVersion(),
        "cases": [],
        "scope": "One segment: fused resize/rotate, profiles and host summaries. "
        "Synchronized wall time; compilation excluded. Synthetic data; "
        "excludes branch labeling, profile postprocessing, and export.",
    }
    for frames, side, varying in ((64, 33, False), (256, 65, False), (256, 65, True)):
        values, mask = movie(frames, side, varying)
        row = {"frames": frames, "native_side": side, "varying_validity": varying}
        row["parity_max_abs_error"] = parity(values, mask)
        if baseline:
            row["main_cpu"] = measure(baseline, "cpu", values, mask, args.repeats)
        row["new_cpu"] = measure(cs, "cpu", values, mask, args.repeats)
        row["new_gpu"] = measure(cs, "cupy", values, mask, args.repeats)
        row["windowed_cpu"] = measure(cs, "cpu", values, mask, args.repeats, windowed=True)
        row["windowed_gpu"] = measure(cs, "cupy", values, mask, args.repeats, windowed=True)
        row["windowed_gpu_speedup"] = (
            row["windowed_cpu"]["median_seconds"] / row["windowed_gpu"]["median_seconds"]
        )
        row["gpu_speedup_vs_new_cpu"] = (
            row["new_cpu"]["median_seconds"] / row["new_gpu"]["median_seconds"]
        )
        report["cases"].append(row)
        args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
        print(json.dumps(row), flush=True)


if __name__ == "__main__":
    main()
