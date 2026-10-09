"""Execute bounded, ordered temporal sampling of scalar or vector segment maps."""

from __future__ import annotations

from collections import deque
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor
from time import perf_counter

import numpy as np

from calculations.compute_backend import optional_cupy_backend
from utils.logger import Logger

from .extraction import extract_segment
from .models import (
    SampledSegment,
    SampledSegmentChunk,
    SampledSegmentChunks,
    SampledSegments,
    SegmentSamplingPlan,
)
from .transforms import interpolate_segments, resample_rotate_segment


def prepare_segment_chunks(
    data_map,
    prepared_topology: SegmentSamplingPlan,
    *,
    spatial_axes: tuple[int, int] = (-2, -1),
    working_memory_mb: float = 512.0,
    worker_count: int | None = None,
    keep_on_device: bool = False,
    transform_mode: str = "fused",
    post_interpolation: Callable[[np.ndarray], np.ndarray] | None = None,
    temporal_halo: int = 0,
    scratch_array_count: int | None = None,
    include_masked_before_rotation: bool = False,
) -> SampledSegmentChunks:
    """Stream bounded temporal chunks of every valid prepared segment.

    Fused mode performs resize and rotation in one affine operation. Staged
    mode interpolates first, filters periodic halo context, trims, then rotates.
    When ``include_masked_before_rotation`` is enabled, each chunk also carries
    a companion following the legacy scientific order ``interpolate/filter ->
    mask -> rotate``. Retained result arrays are outside this scratch budget.
    """

    if transform_mode not in {"fused", "sequential", "staged"}:
        raise ValueError(
            "transform_mode must be 'fused', 'sequential', or 'staged'."
        )
    if transform_mode == "staged" and post_interpolation is None:
        raise ValueError("staged transforms require post_interpolation.")
    if transform_mode != "staged" and post_interpolation is not None:
        raise ValueError("post_interpolation is only valid for staged transforms.")
    if temporal_halo < 0:
        raise ValueError("temporal_halo must be non-negative.")
    if worker_count is not None and worker_count < 1:
        raise ValueError("worker_count must be positive.")

    topology = prepared_topology.native
    frame_count = int(data_map.shape[0])
    valid_indexes = np.argwhere(
        topology.valid_segments & np.isfinite(prepared_topology.rotation_degrees)
    )
    if len(valid_indexes) == 0:
        return
    component_count = _component_count(data_map.shape, spatial_axes)
    workers, chunk_frames = _bounded_chunk_plan(
        work_count=len(valid_indexes),
        frame_count=frame_count,
        native_side=topology.window_side_pixels,
        output_side=prepared_topology.interpolated_masks.shape[-1],
        component_count=component_count,
        working_memory_mb=working_memory_mb,
        temporal_halo=temporal_halo if transform_mode == "staged" else 0,
        scratch_array_count=(
            scratch_array_count
            if scratch_array_count is not None
            else (7 if transform_mode == "staged" else 0)
        ),
        interpolated_array_count=(
            1
            if transform_mode != "fused" or include_masked_before_rotation
            else 0
        ),
        rotated_array_count=2 if include_masked_before_rotation else 1,
        requested_workers=worker_count,
        keep_on_device=keep_on_device,
    )
    jobs = [
        (int(ring), int(branch), start, min(start + chunk_frames, frame_count))
        for ring, branch in valid_indexes
        for start in range(0, frame_count, chunk_frames)
    ]
    Logger.log(
        "Starting bounded topology segment preparation: "
        f"segments={len(valid_indexes)}, workers={workers}, "
        f"chunk_frames={chunk_frames}, halo={temporal_halo}, "
        f"transform_mode={transform_mode}, budget={working_memory_mb:.1f} MiB."
    )

    def prepare_job(job) -> SampledSegmentChunk:
        ring, branch, output_start, output_stop = job
        halo = temporal_halo if transform_mode == "staged" else 0
        if halo and output_stop - output_start + 2 * halo < frame_count:
            context_start = output_start - halo
            context_stop = output_stop + halo
            source = _periodic_frame_slice(data_map, context_start, context_stop)
        elif halo:
            context_start = 0
            context_stop = frame_count
            source = _frame_slice(data_map, context_start, context_stop)
        else:
            context_start = output_start
            context_stop = output_stop
            source = _frame_slice(data_map, context_start, context_stop)
        extracted = extract_segment(
            source,
            topology,
            ring,
            branch,
            spatial_axes=spatial_axes,
        )
        angle = float(prepared_topology.rotation_degrees[ring, branch])
        side = prepared_topology.interpolated_masks.shape[-1]
        interpolated_for_mask = None
        if transform_mode == "fused":
            rotated = resample_rotate_segment(
                extracted,
                angle,
                side,
                return_device=keep_on_device,
            )
            if include_masked_before_rotation:
                interpolated_for_mask = interpolate_segments(extracted, side)
        else:
            interpolated = interpolate_segments(extracted, side)
            if post_interpolation is not None:
                interpolated = np.asarray(
                    post_interpolation(interpolated),
                    dtype=np.float32,
                )
                expected_shape = (*extracted.shape[:-2], side, side)
                if interpolated.shape != expected_shape:
                    raise ValueError(
                        "post_interpolation must preserve all array dimensions."
                    )
            trim_start = output_start - context_start
            trim = slice(trim_start, trim_start + output_stop - output_start)
            interpolated = interpolated[trim]
            interpolated_for_mask = interpolated
            rotated = resample_rotate_segment(
                interpolated,
                angle,
                side,
                return_device=keep_on_device,
            )
        rotated_masked = None
        if include_masked_before_rotation:
            if interpolated_for_mask is None:
                raise RuntimeError("masked transforms require interpolated data.")
            masked = np.asarray(interpolated_for_mask, dtype=np.float32).copy()
            masked[..., ~prepared_topology.interpolated_masks[ring, branch]] = np.nan
            rotated_masked = resample_rotate_segment(
                masked,
                angle,
                side,
                return_device=keep_on_device,
            )
        backend = optional_cupy_backend() if keep_on_device else None
        if backend is not None:
            backend.cupy.cuda.get_current_stream().synchronize()
        return SampledSegmentChunk(
            ring,
            branch,
            slice(output_start, output_stop),
            rotated,
            rotated_masked,
        )

    started = perf_counter()
    if workers == 1:
        for job in jobs:
            yield prepare_job(job)
    else:
        with ThreadPoolExecutor(
            max_workers=workers,
            thread_name_prefix="topology-chunk",
        ) as executor:
            pending = deque()
            job_iter = iter(jobs)
            for _ in range(min(workers, len(jobs))):
                pending.append(executor.submit(prepare_job, next(job_iter)))
            while pending:
                future = pending.popleft()
                yield future.result()
                try:
                    pending.append(executor.submit(prepare_job, next(job_iter)))
                except StopIteration:
                    pass
    Logger.log(
        "Completed bounded topology segment preparation in "
        f"{perf_counter() - started:.2f}s."
    )


def prepare_segments(
    data_map,
    prepared_topology: SegmentSamplingPlan,
    *,
    spatial_axes: tuple[int, int] = (-2, -1),
    worker_count: int = 1,
    keep_on_device: bool = False,
    transform_mode: str = "fused",
) -> SampledSegments:
    """Compatibility helper returning whole segments assembled from chunks."""

    current_index: tuple[int, int] | None = None
    parts: list[object] = []
    backend = optional_cupy_backend() if keep_on_device else None

    def joined():
        if backend is not None and parts and isinstance(parts[0], backend.cupy.ndarray):
            return backend.cupy.concatenate(parts, axis=0)
        return np.concatenate(parts, axis=0)

    for chunk in prepare_segment_chunks(
        data_map,
        prepared_topology,
        spatial_axes=spatial_axes,
        working_memory_mb=float("inf"),
        worker_count=worker_count,
        keep_on_device=keep_on_device,
        transform_mode=transform_mode,
    ):
        index = (chunk.ring_index, chunk.branch_index)
        if current_index is not None and index != current_index:
            raise RuntimeError("segment chunks must be grouped by segment index.")
        current_index = index
        parts.append(chunk.rotated)
        if int(chunk.frame_slice.stop) == int(data_map.shape[0]):
            yield SampledSegment(*current_index, joined())
            current_index = None
            parts = []


def _bounded_chunk_plan(
    *,
    work_count: int,
    frame_count: int,
    native_side: int,
    output_side: int,
    component_count: int,
    working_memory_mb: float,
    temporal_halo: int,
    scratch_array_count: int,
    interpolated_array_count: int,
    rotated_array_count: int,
    requested_workers: int | None,
    keep_on_device: bool,
) -> tuple[int, int]:
    """Choose worker count and output frames under one hard scratch bound."""

    if not np.isfinite(working_memory_mb):
        workers = 1 if keep_on_device else max(1, requested_workers or 1)
        return min(max(work_count, 1), workers), max(frame_count, 1)
    if working_memory_mb <= 0:
        raise ValueError("working_memory_mb must be finite and positive.")
    if frame_count < 1:
        raise ValueError("segment data must contain at least one frame.")
    budget = int(float(working_memory_mb) * 1024**2)
    native_bytes = max(native_side, 1) ** 2 * component_count * 4
    interpolated_bytes = output_side**2 * component_count * 4
    canvas_side = int(output_side * np.sqrt(2.0))
    rotated_bytes = canvas_side**2 * component_count * 4
    context_arrays = max(0, int(scratch_array_count)) + max(
        0, int(interpolated_array_count)
    )
    context_per_frame = native_bytes + context_arrays * interpolated_bytes
    required_context = min(frame_count, 1 + 2 * temporal_halo)
    maximum_extra_context = required_context - 1
    output_bytes = max(1, int(rotated_array_count)) * rotated_bytes
    minimum_worker_bytes = required_context * context_per_frame + output_bytes
    if minimum_worker_bytes > budget:
        raise MemoryError(
            "working_memory_mb is too small for one output frame and its "
            "required temporal filtering context."
        )

    maximum_workers = 1 if keep_on_device else min(
        max(work_count, 1),
        max(1, requested_workers or work_count),
        8,
    )
    workers = min(maximum_workers, max(1, budget // minimum_worker_bytes))
    per_worker_budget = budget // workers
    numerator = per_worker_budget - maximum_extra_context * context_per_frame
    chunk_frames = numerator // max(context_per_frame + output_bytes, 1)
    if chunk_frames < 1:
        workers = 1
        numerator = budget - maximum_extra_context * context_per_frame
        chunk_frames = numerator // max(context_per_frame + output_bytes, 1)
    if chunk_frames < 1:
        raise MemoryError(
            "working_memory_mb is too small for one output frame and its "
            "required temporal filtering context."
        )
    return workers, min(frame_count, int(chunk_frames))


def _component_count(
    shape: tuple[int, ...],
    spatial_axes: tuple[int, int],
) -> int:
    dimension_count = len(shape)
    y_axis, x_axis = (int(axis) % dimension_count for axis in spatial_axes)
    count = 1
    for axis, size in enumerate(shape):
        if axis not in (0, y_axis, x_axis):
            count *= int(size)
    return max(count, 1)


def _frame_slice(data_map, start: int, stop: int):
    slices = [slice(None)] * len(data_map.shape)
    slices[0] = slice(int(start), int(stop))
    return data_map[tuple(slices)]


def _periodic_frame_slice(data_map, start: int, stop: int):
    """Read a possibly wrapped frame interval without fancy dataset indexing."""

    frame_count = int(data_map.shape[0])
    if frame_count < 1:
        raise ValueError("periodic frame slicing requires at least one frame.")
    parts = []
    cursor = int(start)
    stop = int(stop)
    while cursor < stop:
        wrapped_start = cursor % frame_count
        part_length = min(stop - cursor, frame_count - wrapped_start)
        parts.append(
            np.asarray(
                _frame_slice(
                    data_map,
                    wrapped_start,
                    wrapped_start + part_length,
                )
            )
        )
        cursor += part_length
    if len(parts) == 1:
        return parts[0]
    return np.concatenate(parts, axis=0)
