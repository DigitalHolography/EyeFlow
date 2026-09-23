"""Logical execution order for reusable retinal segment topology."""

from __future__ import annotations

from collections.abc import Callable, Iterator, Mapping, MutableMapping
from collections import deque
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, replace
from time import perf_counter
from typing import NamedTuple

import numpy as np
from scipy import ndimage as ndi

from calculations.compute_backend import optional_cupy_backend
from utils.logger import Logger

from .cache import TopologyCacheKey, topology_cache_key
from .geometry import AnnulusGeometry
from .optic_disc import OpticDisc
from .segments import (
    SegmentTopology,
    build_segment_topology,
    extract_segment,
    resize_segment_topology_windows,
)
from .transforms import (
    determine_segment_rotations,
    interpolate_segment_masks,
    interpolate_segments,
    resample_rotate_segment,
    rotate_segment_masks,
    rotate_segments,
)


@dataclass(frozen=True)
class PreparedTopology:
    """Map-independent segment geometry ready for map transformations."""

    topology: SegmentTopology
    rotation_degrees: np.ndarray
    interpolated_masks: np.ndarray
    rotated_masks: np.ndarray


class PreparedSegment(NamedTuple):
    """One uniformly resized and upright segment from a retinal map."""

    ring_index: int
    branch_index: int
    rotated: object


PreparedSegments = Iterator[PreparedSegment]


class PreparedSegmentChunk(NamedTuple):
    """One bounded temporal slice of a uniformly transformed segment."""

    ring_index: int
    branch_index: int
    frame_slice: slice
    rotated: object
    rotated_masked: object | None = None


PreparedSegmentChunks = Iterator[PreparedSegmentChunk]


def prepare_topology(
    vessel_mask,
    optic_disc: OpticDisc,
    settings: AnnulusGeometry,
    *,
    output_side_pixels: int = 128,
    window_size_percentile_kept: float = 0.95,
    window_side_pixels: int | None = None,
) -> PreparedTopology:
    """Build segment geometry, orientations, and uniform masks once."""

    topology = build_segment_topology(
        vessel_mask,
        optic_disc,
        settings,
        window_size_percentile_kept=window_size_percentile_kept,
        window_side_pixels=window_side_pixels,
    )
    rotation_degrees = determine_segment_rotations(topology)
    interpolated_masks = interpolate_segment_masks(
        topology.segment_masks,
        output_side_pixels,
    )
    return PreparedTopology(
        topology=topology,
        rotation_degrees=rotation_degrees,
        interpolated_masks=interpolated_masks,
        rotated_masks=rotate_segment_masks(
            interpolated_masks,
            rotation_degrees,
        ),
    )


def prepare_topologies(
    vessel_masks: Mapping[str, object],
    optic_disc: OpticDisc,
    settings: AnnulusGeometry,
    *,
    source_id: str,
    cache: MutableMapping[TopologyCacheKey, object] | None = None,
    output_side_pixels: int = 128,
    window_size_percentile_kept: float = 0.95,
    window_side_pixels: int | None = None,
) -> dict[str, PreparedTopology]:
    """Prepare consistently sized topology for every named vessel mask."""

    masks = {
        str(name): np.asarray(mask, dtype=bool)
        for name, mask in vessel_masks.items()
    }
    if not masks:
        return {}
    image_shape = next(iter(masks.values())).shape
    if any(mask.shape != image_shape for mask in masks.values()):
        raise ValueError("All vessel masks must have the same spatial shape.")
    disc = optic_disc.mask_for(image_shape)

    if window_side_pixels is not None:
        return {
            name: _cached_topology(
                name,
                mask,
                optic_disc,
                settings,
                source_id=source_id,
                cache=cache,
                output_side_pixels=output_side_pixels,
                window_size_percentile_kept=window_size_percentile_kept,
                window_side_pixels=window_side_pixels,
            )
            for name, mask in masks.items()
        }

    initial = {
        name: _cached_topology(
            name,
            mask,
            optic_disc,
            settings,
            source_id=source_id,
            cache=cache,
            output_side_pixels=output_side_pixels,
            window_size_percentile_kept=window_size_percentile_kept,
            window_side_pixels=None,
        )
        for name, mask in masks.items()
    }
    shared_side = _shared_window_side(initial, window_size_percentile_kept)
    prepared_topologies = {
        name: (
            topology
            if topology.topology.window_side_pixels == shared_side
            else _resize_prepared_topology(
                topology,
                shared_side,
                output_side_pixels,
            )
        )
        for name, topology in initial.items()
    }
    if cache is not None:
        for name, prepared in prepared_topologies.items():
            key = topology_cache_key(
                source_id, name, masks[name], disc, settings,
                optic_disc_center=optic_disc.center,
                output_side_pixels=output_side_pixels,
                window_size_percentile_kept=window_size_percentile_kept,
                window_side_pixels=None,
            )
            # Cache the final joint-window geometry, not just the initial
            # per-vessel geometry, so every pipeline reuses the same objects.
            cache[key] = prepared
    Logger.log(f"Prepared topology: vessels={tuple(masks)}, shared_window={shared_side}px.")
    return prepared_topologies


def _resize_prepared_topology(
    prepared: PreparedTopology,
    window_side_pixels: int,
    output_side_pixels: int,
) -> PreparedTopology:
    topology = resize_segment_topology_windows(
        prepared.topology,
        window_side_pixels,
    )
    interpolated_masks = interpolate_segment_masks(
        topology.segment_masks,
        output_side_pixels,
    )
    return PreparedTopology(
        topology=topology,
        rotation_degrees=prepared.rotation_degrees,
        interpolated_masks=interpolated_masks,
        rotated_masks=rotate_segment_masks(
            interpolated_masks,
            prepared.rotation_degrees,
        ),
    )


def resolve_segment_rotations(
    prepared: PreparedTopology,
    reference_map,
    *,
    working_memory_mb: float = 512.0,
) -> PreparedTopology:
    """Resolve indeterminate centerline angles once from a mean reference map.

    The projection-score fallback is intentionally limited to segments whose
    topology centerline cannot determine an orientation. Segments remain
    unresolved when both methods fail and are skipped by segment iterators.
    """

    rotations = prepared.rotation_degrees.copy()
    unresolved = np.argwhere(
        prepared.topology.valid_segments & ~np.isfinite(rotations)
    )
    if len(unresolved) == 0:
        return prepared
    frame_count = int(reference_map.shape[0])
    if frame_count < 1:
        return prepared
    side = max(prepared.topology.window_side_pixels, 1)
    budget = int(float(working_memory_mb) * 1024**2)
    if budget <= 0:
        raise ValueError("working_memory_mb must be finite and positive.")
    frames_per_block = max(1, min(frame_count, budget // max(side * side * 8, 1)))
    for ring, branch in unresolved:
        total = np.zeros((side, side), dtype=np.float64)
        count = np.zeros((side, side), dtype=np.int64)
        for start in range(0, frame_count, frames_per_block):
            extracted = extract_segment(
                _frame_slice(reference_map, start, min(start + frames_per_block, frame_count)),
                prepared.topology,
                int(ring),
                int(branch),
            )
            finite = np.isfinite(extracted)
            total += np.sum(np.where(finite, extracted, 0.0), axis=0, dtype=np.float64)
            count += np.sum(finite, axis=0, dtype=np.int64)
        mean_native = np.full((side, side), np.nan, dtype=np.float32)
        np.divide(total, count, out=mean_native, where=count > 0)
        mean_image = interpolate_segments(
            mean_native,
            prepared.interpolated_masks.shape[-1],
        )
        mask = prepared.interpolated_masks[int(ring), int(branch)]
        mean_image[~mask] = np.nan
        angle = _projection_rotation_angle(
            mean_image,
            prepared.topology.segment_centers_xy[int(ring), int(branch)],
            prepared.topology.optic_disc_center_xy,
        )
        if np.isfinite(angle):
            rotations[int(ring), int(branch)] = np.float32(angle)

    return replace(
        prepared,
        rotation_degrees=rotations,
        rotated_masks=rotate_segment_masks(prepared.interpolated_masks, rotations),
    )


def _projection_rotation_angle(
    mean_image: np.ndarray,
    segment_center_xy,
    optic_disc_center_xy,
) -> float:
    image = np.asarray(mean_image, dtype=np.float32).copy()
    yy, xx = np.indices(image.shape, dtype=np.float32)
    cy = np.float32((image.shape[0] - 1) / 2.0)
    cx = np.float32((image.shape[1] - 1) / 2.0)
    image[(yy - cy) ** 2 + (xx - cx) ** 2 > min(cx, cy) ** 2] = np.nan
    image = np.nan_to_num(image, nan=0.0)
    image[image < 0] = 0
    if not np.any(image > 0):
        return float("nan")
    loc_x, loc_y = (float(value) for value in segment_center_xy)
    disc_x, disc_y = (float(value) for value in optic_disc_center_xy)
    alpha = np.degrees(np.arctan2(loc_y - disc_y, loc_x - disc_x))
    beta = int(np.mod(90 + np.floor(np.mod(alpha, 360.0) + 0.5), 180))
    angles = np.mod(np.arange(beta - 90, beta + 91), 180)
    scores = np.asarray([_projection_score(image, angle) for angle in angles])
    if not np.any(np.isfinite(scores)):
        return float("nan")
    return float(angles[int(np.nanargmax(scores))])


def _projection_score(image: np.ndarray, angle: float) -> float:
    rotated = ndi.rotate(
        image,
        angle,
        reshape=False,
        order=0,
        mode="constant",
        cval=0.0,
    )
    start = max(int(np.floor(rotated.shape[0] / 3)) - 1, 0)
    stop = int(np.ceil(2 * rotated.shape[0] / 3))
    projection = np.sum(rotated[start:stop], axis=0, dtype=np.float32)
    total = np.sum(projection, dtype=np.float32)
    return float("nan") if total <= 0 else float(np.max(projection) / total)


def prepare_segment_chunks(
    data_map,
    prepared_topology: PreparedTopology,
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
) -> PreparedSegmentChunks:
    """Stream bounded temporal chunks of every valid prepared segment.

    Fused mode performs resize and rotation in one affine operation. Staged
    mode interpolates first, filters the halo context, trims, then rotates.
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

    topology = prepared_topology.topology
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

    def prepare_job(job) -> PreparedSegmentChunk:
        ring, branch, output_start, output_stop = job
        halo = temporal_halo if transform_mode == "staged" else 0
        context_start = max(0, output_start - halo)
        context_stop = min(frame_count, output_stop + halo)
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
            trim = slice(
                output_start - context_start,
                output_stop - context_start,
            )
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
        return PreparedSegmentChunk(
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
    prepared_topology: PreparedTopology,
    *,
    spatial_axes: tuple[int, int] = (-2, -1),
    worker_count: int = 1,
    keep_on_device: bool = False,
    transform_mode: str = "fused",
) -> PreparedSegments:
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
            yield PreparedSegment(*current_index, joined())
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


def _cached_topology(
    vessel_name: str,
    vessel_mask: np.ndarray,
    optic_disc: OpticDisc,
    settings: AnnulusGeometry,
    *,
    source_id: str,
    cache: MutableMapping[TopologyCacheKey, object] | None,
    output_side_pixels: int,
    window_size_percentile_kept: float,
    window_side_pixels: int | None,
) -> PreparedTopology:
    optic_disc_mask = optic_disc.mask_for(vessel_mask.shape)
    key = topology_cache_key(
        source_id,
        vessel_name,
        vessel_mask,
        optic_disc_mask,
        settings,
        optic_disc_center=optic_disc.center,
        output_side_pixels=output_side_pixels,
        window_size_percentile_kept=window_size_percentile_kept,
        window_side_pixels=window_side_pixels,
    )
    if cache is not None:
        found = cache.get(key)
        if found is not None:
            if not isinstance(found, PreparedTopology):
                raise TypeError("Topology cache values must be PreparedTopology.")
            Logger.log(f"Topology cache hit: {vessel_name}; reusing prepared topology.")
            return found

    cache_status = "miss" if cache is not None else "disabled"
    Logger.log(
        f"Topology cache {cache_status}: {vessel_name}; preparing topology."
    )
    prepared = prepare_topology(
        vessel_mask,
        optic_disc,
        settings,
        output_side_pixels=output_side_pixels,
        window_size_percentile_kept=window_size_percentile_kept,
        window_side_pixels=window_side_pixels,
    )
    if cache is not None:
        cache[key] = prepared
    return prepared


def _shared_window_side(
    topologies: Mapping[str, PreparedTopology],
    percentile_kept: float,
) -> int:
    percentile = float(percentile_kept)
    if not 0.0 < percentile <= 1.0:
        raise ValueError("window_size_percentile_kept must be in (0, 1].")

    widths: list[int] = []
    heights: list[int] = []
    for prepared in topologies.values():
        topology = prepared.topology
        for annulus in topology.annulus_masks:
            for branch_id in topology.branch_ids:
                y, x = np.nonzero(annulus & (topology.labels == int(branch_id)))
                if x.size == 0:
                    continue
                widths.append(int(x.max() - x.min() + 1))
                heights.append(int(y.max() - y.min() + 1))

    if not widths:
        return 0

    width = int(np.quantile(widths, percentile, method="higher"))
    height = int(np.quantile(heights, percentile, method="higher"))
    side = max(width, height)
    return side if side % 2 == 1 else side + 1


def _array_summary(values) -> str:
    return (
        f"shape={values.shape}, dtype={values.dtype}, "
        f"allocated={values.nbytes / (1024 ** 3):.2f} GiB"
    )
