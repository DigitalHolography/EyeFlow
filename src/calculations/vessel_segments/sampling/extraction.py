"""Extract padded scalar or component-valued patches using native topology."""

from __future__ import annotations

from time import perf_counter

import numpy as np

from calculations.topology.segments import SegmentTopology
from calculations.topology.windows import window_target_slices
from utils.logger import Logger


def extract_segments(
    data_map,
    topology: SegmentTopology,
    *,
    spatial_axes: tuple[int, int] = (-2, -1),
) -> np.ndarray:
    """Extract one padded square per annulus and branch.

    Spatial axes are moved to the last two output axes. Other axes retain their
    relative order, so a ``(frame, y, x, component)`` vector map produces
    ``(annulus, branch, frame, component, local_y, local_x)`` when called with
    ``spatial_axes=(1, 2)``.
    """

    shape = tuple(int(size) for size in data_map.shape)
    y_axis, x_axis = _normalized_spatial_axes(len(shape), spatial_axes)
    assert (shape[y_axis], shape[x_axis]) == topology.spatial_shape
    nonspatial_shape = tuple(
        size for axis, size in enumerate(shape) if axis not in (y_axis, x_axis)
    )
    ring_count, branch_count = topology.segment_centers_xy.shape[:2]
    side = topology.window_side_pixels
    output_shape = (ring_count, branch_count, *nonspatial_shape, side, side)
    output_gib = np.prod(output_shape, dtype=np.int64) * 4 / (1024 ** 3)
    Logger.log(
        f"Allocating extracted segment array: shape={output_shape}, "
        f"allocated={output_gib:.2f} GiB."
    )
    extracted = np.full(
        output_shape,
        np.nan,
        dtype=np.float32,
    )
    if side == 0:
        return extracted

    target_prefix = (slice(None),) * len(nonspatial_shape)
    valid_indexes = np.argwhere(topology.valid_segments)
    progress_step = max(1, len(valid_indexes) // 10)
    read_seconds = 0.0
    extraction_started = perf_counter()
    for work_index, (ring_index, branch_index) in enumerate(valid_indexes, start=1):
        bounds = topology.window_bounds_xyxy[ring_index, branch_index]
        center = topology.segment_centers_xy[ring_index, branch_index]
        source_slices = [slice(None)] * len(shape)
        source_slices[x_axis] = slice(int(bounds[0]), int(bounds[1]))
        source_slices[y_axis] = slice(int(bounds[2]), int(bounds[3]))
        read_started = perf_counter()
        source = np.asarray(data_map[tuple(source_slices)], dtype=np.float32)
        read_seconds += perf_counter() - read_started
        source = np.moveaxis(source, (y_axis, x_axis), (-2, -1))
        target_y, target_x = window_target_slices(bounds, center, side)
        extracted[(
            int(ring_index),
            int(branch_index),
            *target_prefix,
            target_y,
            target_x,
        )] = source
        if work_index % progress_step == 0 or work_index == len(valid_indexes):
            Logger.log(
                f"Segment extraction progress: {work_index}/{len(valid_indexes)} "
                f"in {perf_counter() - extraction_started:.2f}s."
            )
    Logger.log(
        f"Segment source slicing accounted for {read_seconds:.2f}s across "
        f"{len(valid_indexes)} segment reads."
    )
    return extracted


def extract_segment(
    data_map,
    topology: SegmentTopology,
    ring_index: int,
    branch_index: int,
    *,
    spatial_axes: tuple[int, int] = (-2, -1),
) -> np.ndarray:
    """Extract one padded segment while retaining all non-spatial axes."""

    shape = tuple(int(size) for size in data_map.shape)
    y_axis, x_axis = _normalized_spatial_axes(len(shape), spatial_axes)
    assert (shape[y_axis], shape[x_axis]) == topology.spatial_shape
    nonspatial_shape = tuple(
        size for axis, size in enumerate(shape) if axis not in (y_axis, x_axis)
    )
    side = topology.window_side_pixels
    extracted = np.full((*nonspatial_shape, side, side), np.nan, dtype=np.float32)
    index = (int(ring_index), int(branch_index))
    if side == 0 or not topology.valid_segments[index]:
        return extracted

    bounds = topology.window_bounds_xyxy[index]
    center = topology.segment_centers_xy[index]
    source_slices = [slice(None)] * len(shape)
    source_slices[x_axis] = slice(int(bounds[0]), int(bounds[1]))
    source_slices[y_axis] = slice(int(bounds[2]), int(bounds[3]))
    source = np.asarray(data_map[tuple(source_slices)], dtype=np.float32)
    source = np.moveaxis(source, (y_axis, x_axis), (-2, -1))
    target_y, target_x = window_target_slices(bounds, center, side)
    extracted[..., target_y, target_x] = source
    return extracted


def _normalized_spatial_axes(
    dimension_count: int,
    spatial_axes: tuple[int, int],
) -> tuple[int, int]:
    if dimension_count < 2:
        raise ValueError("data_map must contain two spatial axes.")
    y_axis, x_axis = (int(axis) % dimension_count for axis in spatial_axes)
    if y_axis == x_axis:
        raise ValueError("spatial_axes must identify two different axes.")
    return y_axis, x_axis
