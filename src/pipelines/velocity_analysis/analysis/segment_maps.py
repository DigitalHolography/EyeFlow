"""Per-beat rotated velocity maps and segment-mask output packing."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from itertools import pairwise

import numpy as np

from calculations.blood_flow_velocity.signal_analysis.per_beat._signal_utils import (
    normalize_cycle_boundaries,
)
from calculations.compute_backend import optional_cupy_backend
from calculations.math import interpft_axis0, next_power_of_two
from runtime_limits import cap_parallel_jobs
from utils.logger import Logger

_MAX_PARALLEL_SEGMENT_INTERPOLATIONS = 8
_GPU_INTERPOLATION_SCRATCH_BYTES = 64 * 1024**2
# Short temporal FFTs are faster on the threaded CPU once host transfers and
# serialization are included.
_MIN_GPU_MAP_TIME_COUNT = 64


def prepare_segment_velocity_maps_per_beat(
    artery_segments,
    vein_segments,
    cycle_boundary_indexes,
    *,
    index_base: int = 0,
) -> tuple[np.ndarray | None, np.ndarray | None]:
    """Interpolate each vessel's maps once for reuse by output products."""

    return (
        _prepare_vessel_velocity_maps_per_beat(
            artery_segments,
            cycle_boundary_indexes,
            index_base=index_base,
        ),
        _prepare_vessel_velocity_maps_per_beat(
            vein_segments,
            cycle_boundary_indexes,
            index_base=index_base,
        ),
    )


def _prepare_vessel_velocity_maps_per_beat(
    segments,
    cycle_boundary_indexes,
    *,
    index_base: int,
) -> np.ndarray | None:
    if segments is None:
        return None
    profile = segments.profile
    retained = profile.require_maps()
    maps = np.asarray(retained.values)
    compact_arguments = (
        {
            "segment_indexes": retained.indexes,
            "radius_count": profile.segment_shape[0],
            "branch_count": profile.segment_shape[1],
        }
        if maps.ndim == 4
        else {}
    )
    return interpolate_velocity_maps_per_beat(
        maps,
        cycle_boundary_indexes,
        index_base=index_base,
        **compact_arguments,
    )


def interpolate_velocity_maps_per_beat(
    velocity_maps: np.ndarray,
    cycle_boundary_indexes,
    *,
    index_base: int = 0,
    segment_indexes: np.ndarray | None = None,
    radius_count: int | None = None,
    branch_count: int | None = None,
) -> np.ndarray:
    """Interpolate dense or compact maps to the serialized output layout.

    Compact input has shape ``(segment, frame, y, x)`` and is paired with a
    ``(segment, 2)`` radius/branch index table. Dense legacy input with shape
    ``(radius, branch, frame, y, x)`` remains supported.
    """
    maps = np.asarray(velocity_maps, dtype=np.float32)
    if maps.ndim == 5:
        (
            resolved_radius_count,
            resolved_branch_count,
            frame_count,
            y_count,
            x_count,
        ) = maps.shape
        indexed_rows = [
            (radius_index, branch_index, maps[radius_index, branch_index])
            for radius_index in range(resolved_radius_count)
            for branch_index in range(resolved_branch_count)
        ]
    elif maps.ndim == 4:
        if segment_indexes is None or radius_count is None or branch_count is None:
            raise ValueError(
                "Compact velocity maps require segment_indexes, radius_count, "
                "and branch_count."
            )
        resolved_radius_count = int(radius_count)
        resolved_branch_count = int(branch_count)
        if resolved_radius_count < 0 or resolved_branch_count < 0:
            raise ValueError("radius_count and branch_count must be non-negative.")
        indexes = np.asarray(segment_indexes, dtype=np.int32)
        if indexes.shape != (maps.shape[0], 2):
            raise ValueError(
                "segment_indexes must have shape (segment, 2) matching maps."
            )
        if indexes.size:
            if (
                np.any(indexes[:, 0] < 0)
                or np.any(indexes[:, 0] >= resolved_radius_count)
                or np.any(indexes[:, 1] < 0)
                or np.any(indexes[:, 1] >= resolved_branch_count)
            ):
                raise ValueError("segment_indexes contain an out-of-range index.")
            if np.unique(indexes, axis=0).shape[0] != indexes.shape[0]:
                raise ValueError("segment_indexes must not contain duplicates.")
        _, frame_count, y_count, x_count = maps.shape
        indexed_rows = [
            (int(radius_index), int(branch_index), maps[row_index])
            for row_index, (radius_index, branch_index) in enumerate(indexes)
        ]
    else:
        raise ValueError(
            "velocity maps must have dense (radius, branch, frame, y, x) or "
            "compact (segment, frame, y, x) shape."
        )

    boundaries = normalize_cycle_boundaries(
        cycle_boundary_indexes,
        frame_count,
        index_base=index_base,
    )
    beat_count = boundaries.size - 1
    time_count = next_power_of_two(int(np.max(np.diff(boundaries))))
    output = np.full(
        (
            x_count,
            y_count,
            time_count,
            beat_count,
            resolved_branch_count,
            resolved_radius_count,
        ),
        np.nan,
        dtype=np.float32,
    )
    backend = optional_cupy_backend()
    if time_count < _MIN_GPU_MAP_TIME_COUNT:
        backend = None
    if backend is not None:
        Logger.log("Per-beat segment velocity-map interpolation backend: CuPy/CUDA.")

    def interpolate_segment(indexed_row) -> None:
        radius_index, branch_index, segment_maps = indexed_row
        for beat_index in range(beat_count):
            start = int(boundaries[beat_index])
            stop = int(boundaries[beat_index + 1]) + 1
            interpolated = interpft_axis0(
                segment_maps[start:stop],
                time_count + 1,
            )[:-1]
            output[:, :, :, beat_index, branch_index, radius_index] = interpolated.transpose(
                2, 1, 0
            )

    if backend is not None:
        # One GPU job at a time bounds scratch across segments. Batch pixels
        # across beats instead of allocating a dense six-dimensional device
        # output; serialized arrays stay float32 NumPy for HDF5 and movie writers.
        for radius_index, branch_index, segment_maps in indexed_rows:
            try:
                _interpolate_segment_gpu(
                    segment_maps,
                    boundaries,
                    time_count,
                    output,
                    radius_index,
                    branch_index,
                    backend.cupy,
                )
            except Exception as exc:  # noqa: BLE001 - optional CUDA kernels can fail in many ways
                Logger.log_debug(
                    "CuPy segment-map interpolation failed; using CPU fallback: "
                    f"{type(exc).__name__}: {exc}"
                )
                interpolate_segment((radius_index, branch_index, segment_maps))
        return output

    worker_count = _segment_map_worker_count(len(indexed_rows))
    if worker_count == 1:
        for indexed_row in indexed_rows:
            interpolate_segment(indexed_row)
    else:
        with ThreadPoolExecutor(
            max_workers=worker_count,
            thread_name_prefix="segment-map",
        ) as executor:
            for _ in executor.map(interpolate_segment, indexed_rows):
                pass
    return output


def _interpolate_segment_gpu(
    maps, boundaries, time_count, output, radius_index, branch_index, cupy,
) -> None:
    """Fourier-resample bounded pixel batches, including shared beat endpoints."""

    y_count, x_count = maps.shape[-2:]
    lengths = np.diff(boundaries) + 1
    if np.all(lengths == time_count + 1):
        # Avoid device transfers when every beat requires only a copy.
        for beat_index, start in enumerate(boundaries[:-1]):
            output[:, :, :, beat_index, branch_index, radius_index] = maps[
                int(start):int(start) + time_count
            ].transpose(2, 1, 0)
        return

    beat_count = lengths.size
    # Upload each spatial batch once, retaining all of its resampled beats until
    # one download. Budget the resident source/results plus Fourier temporaries;
    # the persistent six-dimensional host output is outside this scratch budget.
    bytes_per_pixel = (
        4 * (maps.shape[0] + time_count * beat_count)
        + 32 * (int(lengths.max()) + time_count + 1)
    )
    batch_width = _GPU_INTERPOLATION_SCRATCH_BYTES // (bytes_per_pixel * max(y_count, 1))
    if batch_width < 1:
        raise MemoryError("One segment-map column exceeds the GPU interpolation scratch budget.")
    for x_start in range(0, x_count, batch_width):
        x_stop = min(x_count, x_start + batch_width)
        source = cupy.asarray(maps[:, :, x_start:x_stop], dtype=cupy.float32)
        results = cupy.empty((time_count, beat_count, y_count, x_stop - x_start), dtype=cupy.float32)
        for beat_index, (start, stop) in enumerate(pairwise(boundaries)):
            results[:, beat_index] = interpft_axis0(
                source[int(start):int(stop) + 1],
                time_count + 1,
                array_module=cupy,
            )[:-1]
        output[x_start:x_stop, :, :, :, branch_index, radius_index] = (
            cupy.asnumpy(results).transpose(3, 2, 0, 1)
        )
        del source, results


def _segment_map_worker_count(segment_count: int) -> int:
    if segment_count <= 1:
        return 1
    return min(
        segment_count,
        cap_parallel_jobs(_MAX_PARALLEL_SEGMENT_INTERPOLATIONS),
    )
