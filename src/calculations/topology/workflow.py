"""Logical execution order for reusable retinal segment topology."""

from __future__ import annotations

from collections.abc import Iterator, Mapping, MutableMapping
from collections import deque
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from time import perf_counter
from typing import NamedTuple

import numpy as np

from calculations.compute_backend import optional_cupy_backend
from utils.logger import Logger

from .cache import TopologyCacheKey, topology_cache_key
from .geometry import SegmentRingSettings
from .segments import (
    SegmentTopology,
    build_segment_topology,
    competing_segment_masks,
    extract_segment,
    resize_segment_topology_windows,
)
from .transforms import (
    determine_segment_rotations,
    interpolate_segment_masks,
    resample_rotate_segment,
    rotate_segment_masks,
)


@dataclass(frozen=True)
class PreparedTopology:
    """Map-independent segment geometry ready for map transformations."""

    topology: SegmentTopology
    rotation_degrees: np.ndarray
    interpolated_masks: np.ndarray
    rotated_masks: np.ndarray
    rotated_competing_masks: np.ndarray | None = None


class PreparedSegment(NamedTuple):
    """One uniformly resized and upright segment from a retinal map."""

    ring_index: int
    branch_index: int
    rotated: object


PreparedSegments = Iterator[PreparedSegment]


def prepare_topology(
    vessel_mask,
    optic_disc_mask,
    settings: SegmentRingSettings,
    *,
    output_side_pixels: int = 128,
    window_size_percentile_kept: float = 0.95,
    window_side_pixels: int | None = None,
    competing_vessel_mask=None,
) -> PreparedTopology:
    """Build segment geometry, orientations, and uniform masks once."""

    topology = build_segment_topology(
        vessel_mask,
        optic_disc_mask,
        settings,
        window_size_percentile_kept=window_size_percentile_kept,
        window_side_pixels=window_side_pixels,
    )
    rotation_degrees = determine_segment_rotations(topology)
    interpolated_masks = interpolate_segment_masks(
        topology.segment_masks,
        output_side_pixels,
    )
    competing_masks = competing_segment_masks(
        topology,
        vessel_mask,
        competing_vessel_mask,
    )
    interpolated_competing_masks = interpolate_segment_masks(
        competing_masks,
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
        rotated_competing_masks=rotate_segment_masks(
            interpolated_competing_masks,
            rotation_degrees,
        ),
    )


def prepare_topologies(
    vessel_masks: Mapping[str, object],
    optic_disc_mask,
    settings: SegmentRingSettings,
    *,
    source_id: str,
    cache: MutableMapping[TopologyCacheKey, object] | None = None,
    optic_disc_center=None,
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
    disc = _resolved_optic_disc_mask(
        optic_disc_mask,
        optic_disc_center,
        image_shape,
    )
    competing_masks = {}
    for name in masks:
        competing = np.zeros(image_shape, dtype=bool)
        for other_name, mask in masks.items():
            if other_name != name:
                np.logical_or(competing, mask, out=competing)
        competing_masks[name] = competing

    if window_side_pixels is not None:
        return {
            name: _cached_topology(
                name,
                mask,
                disc,
                settings,
                source_id=source_id,
                cache=cache,
                output_side_pixels=output_side_pixels,
                window_size_percentile_kept=window_size_percentile_kept,
                window_side_pixels=window_side_pixels,
                competing_vessel_mask=competing_masks[name],
            )
            for name, mask in masks.items()
        }

    initial = {
        name: _cached_topology(
            name,
            mask,
            disc,
            settings,
            source_id=source_id,
            cache=cache,
            output_side_pixels=output_side_pixels,
            window_size_percentile_kept=window_size_percentile_kept,
            window_side_pixels=None,
            competing_vessel_mask=competing_masks[name],
        )
        for name, mask in masks.items()
    }
    shared_side = _shared_window_side(initial, window_size_percentile_kept)
    return {
        name: (
            topology
            if topology.topology.window_side_pixels == shared_side
            else _resize_prepared_topology(
                topology,
                shared_side,
                output_side_pixels,
                masks[name],
                competing_masks[name],
            )
        )
        for name, topology in initial.items()
    }


def _resize_prepared_topology(
    prepared: PreparedTopology,
    window_side_pixels: int,
    output_side_pixels: int,
    vessel_mask: np.ndarray,
    competing_vessel_mask: np.ndarray,
) -> PreparedTopology:
    topology = resize_segment_topology_windows(
        prepared.topology,
        window_side_pixels,
    )
    interpolated_masks = interpolate_segment_masks(
        topology.segment_masks,
        output_side_pixels,
    )
    interpolated_competing_masks = interpolate_segment_masks(
        competing_segment_masks(topology, vessel_mask, competing_vessel_mask),
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
        rotated_competing_masks=rotate_segment_masks(
            interpolated_competing_masks,
            prepared.rotation_degrees,
        ),
    )


def prepare_segments(
    data_map,
    prepared_topology: PreparedTopology,
    *,
    spatial_axes: tuple[int, int] = (-2, -1),
    worker_count: int = 1,
    keep_on_device: bool = False,
) -> Iterator[PreparedSegment]:
    """Yield prepared stacks for valid segments without dense materialization."""

    total_started = perf_counter()
    topology = prepared_topology.topology
    valid_indexes = np.argwhere(
        topology.valid_segments & np.isfinite(prepared_topology.rotation_degrees)
    )
    progress_step = max(1, len(valid_indexes) // 10)
    Logger.log(
        "Starting streamed topology segment preparation: "
        f"source_shape={tuple(int(size) for size in data_map.shape)}, "
        f"source_type={type(data_map).__name__}, "
        f"valid_segments={len(valid_indexes)}."
    )
    if worker_count < 1:
        raise ValueError("worker_count must be positive.")

    def prepare_index(index) -> tuple[PreparedSegment, float, float]:
        ring_index, branch_index = index
        ring = int(ring_index)
        branch = int(branch_index)
        started = perf_counter()
        extracted = extract_segment(
            data_map,
            topology,
            ring,
            branch,
            spatial_axes=spatial_axes,
        )
        extraction_seconds = perf_counter() - started
        started = perf_counter()
        rotated = resample_rotate_segment(
            extracted,
            float(prepared_topology.rotation_degrees[ring, branch]),
            prepared_topology.interpolated_masks.shape[-1],
            return_device=keep_on_device,
        )
        backend = optional_cupy_backend() if keep_on_device else None
        if backend is not None:
            backend.cupy.cuda.get_current_stream().synchronize()
        transform_seconds = perf_counter() - started
        return (
            PreparedSegment(
                ring_index=ring,
                branch_index=branch,
                rotated=rotated,
            ),
            extraction_seconds,
            transform_seconds,
        )

    def prepared_results():
        if worker_count == 1:
            for index in valid_indexes:
                yield prepare_index(index)
            return
        with ThreadPoolExecutor(
            max_workers=min(worker_count, len(valid_indexes)),
            thread_name_prefix="topology-segment",
        ) as executor:
            pending = deque()
            indexes = iter(valid_indexes)
            for _ in range(min(worker_count, len(valid_indexes))):
                pending.append(executor.submit(prepare_index, next(indexes)))
            while pending:
                future = pending.popleft()
                yield future.result()
                try:
                    index = next(indexes)
                except StopIteration:
                    continue
                pending.append(executor.submit(prepare_index, index))

    for work_index, result in enumerate(prepared_results(), start=1):
        prepared, extraction_seconds, transform_seconds = result
        if work_index % progress_step == 0 or work_index == len(valid_indexes):
            Logger.log(
                f"Streamed segment progress: {work_index}/{len(valid_indexes)}; "
                f"index=({prepared.ring_index}, {prepared.branch_index}), "
                f"extraction={extraction_seconds:.2f}s, "
                f"fused_transform={transform_seconds:.2f}s; "
                f"rotated={_array_summary(prepared.rotated)}."
            )
        yield prepared

    Logger.log(
        "Completed streamed topology segment preparation in "
        f"{perf_counter() - total_started:.2f}s."
    )


def _cached_topology(
    vessel_name: str,
    vessel_mask: np.ndarray,
    optic_disc_mask: np.ndarray,
    settings: SegmentRingSettings,
    *,
    source_id: str,
    cache: MutableMapping[TopologyCacheKey, object] | None,
    output_side_pixels: int,
    window_size_percentile_kept: float,
    window_side_pixels: int | None,
    competing_vessel_mask: np.ndarray,
) -> PreparedTopology:
    key = topology_cache_key(
        source_id,
        vessel_name,
        vessel_mask,
        optic_disc_mask,
        settings,
        output_side_pixels=output_side_pixels,
        window_size_percentile_kept=window_size_percentile_kept,
        window_side_pixels=window_side_pixels,
        competing_vessel_mask=competing_vessel_mask,
    )
    if cache is not None:
        found = cache.get(key)
        if found is not None:
            if not isinstance(found, PreparedTopology):
                raise TypeError("Topology cache values must be PreparedTopology.")
            return found

    prepared = prepare_topology(
        vessel_mask,
        optic_disc_mask,
        settings,
        output_side_pixels=output_side_pixels,
        window_size_percentile_kept=window_size_percentile_kept,
        window_side_pixels=window_side_pixels,
        competing_vessel_mask=competing_vessel_mask,
    )
    if cache is not None:
        cache[key] = prepared
    return prepared


def _resolved_optic_disc_mask(
    optic_disc_mask,
    optic_disc_center,
    image_shape: tuple[int, int],
) -> np.ndarray:
    if optic_disc_mask is not None:
        disc = np.asarray(optic_disc_mask, dtype=bool)
        if disc.shape != image_shape:
            raise ValueError(
                f"optic_disc_mask must have shape {image_shape}, got {disc.shape}."
            )
        return disc

    disc = np.zeros(image_shape, dtype=bool)
    if optic_disc_center is None:
        y = image_shape[0] // 2
        x = image_shape[1] // 2
    else:
        center = np.asarray(optic_disc_center, dtype=np.float32).reshape(-1)
        if center.size < 2 or not np.all(np.isfinite(center[:2])):
            y = image_shape[0] // 2
            x = image_shape[1] // 2
        else:
            x = int(np.clip(np.floor(center[0] + 0.5), 0, image_shape[1] - 1))
            y = int(np.clip(np.floor(center[1] + 0.5), 0, image_shape[0] - 1))
    disc[y, x] = True
    return disc


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
