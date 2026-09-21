"""Benchmark bounded fused and staged cross-section processing.

Run a CPU-only report from the repository with::

    python benchmarks/benchmark_cross_section_gpu.py --cpu-only --output report.json

When a usable CuPy/CUDA installation is present, the default also verifies
CPU/GPU parity and records CUDA allocator peak usage.
"""

from __future__ import annotations

import argparse
import gc
import importlib
import json
import os
import statistics
import sys
import time
import tracemalloc
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import scipy

from calculations.compute_backend import optional_cupy_backend
from calculations.topology import (
    PreparedTopology,
    SegmentRingSettings,
    SegmentTopology,
    interpolate_segment_masks,
    prepare_segment_chunks,
    rotate_segment_masks,
)
from pipelines.waveform_velocity.spatial_gradient_profiles import (
    _spatial_gradient_chain,
)

cs = importlib.import_module(
    "calculations.blood_flow_velocity.cross_section.generate_cross_section_signals"
)


def movie(frames: int, side: int, varying: bool):
    y, x = np.mgrid[:side, :side].astype(np.float32)
    center = (side - 1) / 2
    mask = (np.abs(x - center) <= side / 6) & (np.abs(y - center) <= side / 3)
    spatial = 2 + 10 * np.maximum(0, 1 - ((x - center) / (side / 5)) ** 2)
    temporal = 1 + 0.3 * np.sin(
        np.arange(frames, dtype=np.float32) * 2 * np.pi / 32
    )
    values = temporal[:, None, None] * spatial[None]
    values[:, :2] = np.nan
    if varying:
        values[::13, side // 2, side // 2] = np.nan
    return values.astype(np.float32), mask


def select_backend(mode: str) -> None:
    os.environ["EYEFLOW_COMPUTE_BACKEND"] = mode
    optional_cupy_backend.cache_clear()


def _prepared_topology(mask: np.ndarray) -> PreparedTopology:
    side = mask.shape[-1]
    branch_ids = np.asarray([1], dtype=np.int32)
    branch_stages = type("Stages", (), {"vessel": mask})()
    branch_identity = type(
        "Branches",
        (),
        {
            "branch_ids": branch_ids,
            "labels": mask.astype(np.int32),
            "stages": branch_stages,
        },
    )()
    geometry = SegmentTopology(
        spatial_shape=mask.shape,
        optic_disc_center_xy=(0.0, 0.0),
        labels=mask.astype(np.int32),
        centerline=mask,
        branch_ids=branch_ids,
        annulus_masks=np.ones((1, *mask.shape), dtype=bool),
        segment_masks=mask[None, None],
        segment_centers_xy=np.asarray(
            [[[(side - 1) / 2, (side - 1) / 2]]], dtype=np.float32
        ),
        window_bounds_xyxy=np.asarray([[[0, side, 0, side]]], dtype=np.int32),
        window_side_pixels=side,
        optic_disc_mask=np.zeros(mask.shape, dtype=bool),
        ring_settings=SegmentRingSettings(0.0, 1.0, 1.0, 1),
        branch_identity=branch_identity,
    )
    angles = np.asarray([[-59.0]], dtype=np.float32)
    interpolated = interpolate_segment_masks(mask[None, None])
    empty = np.zeros_like(interpolated)
    return PreparedTopology(
        topology=geometry,
        rotation_degrees=angles,
        interpolated_masks=interpolated,
        rotated_masks=rotate_segment_masks(interpolated, angles),
        interpolated_competing_masks=empty,
        rotated_competing_masks=rotate_segment_masks(empty, angles),
    )


def topology_measurement(
    values: np.ndarray,
    mask: np.ndarray,
    *,
    transform_mode: str,
    working_memory_mb: float,
):
    prepared = _prepared_topology(mask)
    backend = optional_cupy_backend()
    staged = transform_mode == "staged"
    chunks = prepare_segment_chunks(
        values,
        prepared,
        working_memory_mb=working_memory_mb,
        keep_on_device=backend is not None,
        transform_mode=transform_mode,
        post_interpolation=_spatial_gradient_chain if staged else None,
        temporal_halo=6 if staged else 0,
        scratch_array_count=7 if staged else None,
    )
    settings = cs.CrossSectionSignalSettings(
        pixel_size_mm=0.01,
        working_memory_mb=working_memory_mb,
    )
    return cs._generate_cross_section_signals_from_prepared(
        values,
        prepared,
        chunks,
        prepared.topology.ring_settings,
        settings,
        retain_velocity_maps=False,
    )


def _synchronize() -> None:
    backend = optional_cupy_backend()
    if backend is not None:
        backend.cupy.cuda.get_current_stream().synchronize()


def measure(
    mode: str,
    values: np.ndarray,
    mask: np.ndarray,
    repeats: int,
    *,
    transform_mode: str,
    working_memory_mb: float,
):
    select_backend(mode)
    backend = optional_cupy_backend()
    if backend is not None:
        backend.cupy.get_default_memory_pool().free_all_blocks()

    def call():
        return topology_measurement(
            values,
            mask,
            transform_mode=transform_mode,
            working_memory_mb=working_memory_mb,
        )

    warm = call()
    _synchronize()
    del warm
    times = []
    peak_python_bytes = 0
    for _ in range(repeats):
        gc.collect()
        _synchronize()
        tracemalloc.start()
        start = time.perf_counter()
        result = call()
        _synchronize()
        times.append(time.perf_counter() - start)
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        peak_python_bytes = max(peak_python_bytes, peak)
        del result
    return {
        "seconds": times,
        "median_seconds": statistics.median(times),
        "python_tracemalloc_peak_mib": peak_python_bytes / 1024**2,
        "cuda_pool_reserved_mib": (
            backend.cupy.get_default_memory_pool().total_bytes() / 1024**2
            if backend is not None
            else None
        ),
    }


def parity(
    values: np.ndarray,
    mask: np.ndarray,
    *,
    transform_mode: str,
    working_memory_mb: float,
):
    results = []
    for mode in ("cpu", "cupy"):
        select_backend(mode)
        results.append(
            topology_measurement(
                values,
                mask,
                transform_mode=transform_mode,
                working_memory_mb=working_memory_mb,
            )
        )
    errors = {}
    for name in (
        "velocity",
        "safe_velocity",
        "velocity_profiles",
        "transverse_velocity_profiles_masked",
        "longitudinal_velocity_profiles_unmasked",
        "longitudinal_velocity_profiles_masked",
    ):
        cpu, gpu = (getattr(result, name) for result in results)
        np.testing.assert_allclose(gpu, cpu, rtol=1e-5, atol=1e-5, equal_nan=True)
        difference = np.abs(cpu - gpu)
        errors[name] = (
            float(np.nanmax(difference)) if np.any(np.isfinite(difference)) else 0.0
        )
    return errors


def _cuda_metadata() -> dict[str, object] | None:
    select_backend("auto")
    backend = optional_cupy_backend()
    if backend is None:
        return None
    cupy = backend.cupy
    properties = cupy.cuda.runtime.getDeviceProperties(0)
    name = properties["name"]
    return {
        "gpu": name.decode() if isinstance(name, bytes) else str(name),
        "cupy": cupy.__version__,
        "cuda_runtime": cupy.cuda.runtime.runtimeGetVersion(),
        "cuda_driver": cupy.cuda.runtime.driverGetVersion(),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--working-memory-mb", type=float, default=64.0)
    parser.add_argument("--cpu-only", action="store_true")
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args()
    cuda = None if args.cpu_only else _cuda_metadata()
    report = {
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "working_memory_mb": args.working_memory_mb,
        "cuda": cuda,
        "cases": [],
        "scope": (
            "One segment through bounded fused or staged transform, profile "
            "reduction, and host summaries. Staged mode uses the production "
            "two-moving-average spatial-gradient chain with a six-frame halo."
        ),
    }
    cases = (
        ((12, 17, True),)
        if args.quick
        else ((64, 33, False), (128, 65, False), (128, 65, True))
    )
    for frames, side, varying in cases:
        values, mask = movie(frames, side, varying)
        row: dict[str, object] = {
            "frames": frames,
            "native_side": side,
            "varying_validity": varying,
            "paths": {},
        }
        for transform_mode in ("fused", "staged"):
            measurements = {
                "cpu": measure(
                    "cpu",
                    values,
                    mask,
                    args.repeats,
                    transform_mode=transform_mode,
                    working_memory_mb=args.working_memory_mb,
                )
            }
            if cuda is not None:
                measurements["gpu"] = measure(
                    "cupy",
                    values,
                    mask,
                    args.repeats,
                    transform_mode=transform_mode,
                    working_memory_mb=args.working_memory_mb,
                )
                measurements["parity_max_abs_error"] = parity(
                    values,
                    mask,
                    transform_mode=transform_mode,
                    working_memory_mb=args.working_memory_mb,
                )
            row["paths"][transform_mode] = measurements
        report["cases"].append(row)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
        print(json.dumps(row), flush=True)


if __name__ == "__main__":
    main()
