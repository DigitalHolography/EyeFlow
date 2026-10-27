"""Extract padded scalar or component-valued patches using native topology."""

from __future__ import annotations

import numpy as np

from calculations.topology.segments import SegmentTopology
from calculations.topology.windows import window_target_slices



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
