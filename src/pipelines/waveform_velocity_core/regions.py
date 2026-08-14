"""Shared retinal quadrant geometry for velocity and metric outputs."""

from __future__ import annotations

import numpy as np

REGION_NAMES = (
    "north_west",
    "north_east",
    "south_west",
    "south_east",
)

QUADRANTS_GROUP_NAME = "Quadrants"

_TRIGONOMETRIC_QUADRANT_ORDER = (
    "north_east",
    "north_west",
    "south_west",
    "south_east",
)
_TRIGONOMETRIC_QUADRANT_INDICES = tuple(
    REGION_NAMES.index(name) for name in _TRIGONOMETRIC_QUADRANT_ORDER
)
_EYEFLOW_SPATIAL_Y_INVERTED = True


def region_membership(
    branch_ids: np.ndarray,
    branch_label_map: np.ndarray,
    segment_center_xy: np.ndarray,
    optic_disc_center: np.ndarray,
) -> np.ndarray:
    """Assign every branch/radius to its majority quadrant."""
    if branch_label_map.ndim != 2:
        raise ValueError(
            f"BranchLabelMap must have shape (y, x), got {branch_label_map.shape}."
        )
    if segment_center_xy.ndim != 3 or segment_center_xy.shape[2] != 2:
        raise ValueError(
            "SegmentCenterXY must have shape (branch, radius, 2), got "
            f"{segment_center_xy.shape}."
        )

    n_branches, n_radii = segment_center_xy.shape[:2]
    if branch_ids.size != n_branches:
        raise ValueError("BranchIds and SegmentCenterXY must have the same branch count.")

    height, width = branch_label_map.shape
    pixel_y, pixel_x = np.indices((height, width), dtype=float)
    west = pixel_x < optic_disc_center[0]
    east = ~west
    south = pixel_y < optic_disc_center[1]
    north = ~south
    quadrant_masks = np.asarray(
        (north & west, north & east, south & west, south & east),
        dtype=bool,
    )

    area_by_quadrant = np.zeros((n_branches, 4), dtype=np.float32)
    for branch_index, branch_id in enumerate(branch_ids):
        branch_pixels = branch_label_map == branch_id
        area_by_quadrant[branch_index] = np.asarray(
            [np.count_nonzero(branch_pixels & mask) for mask in quadrant_masks],
            dtype=np.float32,
        )

    if not np.all(np.any(area_by_quadrant, axis=1)):
        raise ValueError(
            "BranchLabelMap must contain at least one pixel for every branch."
        )

    chosen_quadrants = np.full(n_branches, -1, dtype=np.int32)
    priority_order = np.asarray(_TRIGONOMETRIC_QUADRANT_INDICES, dtype=int)
    for branch_index, area in enumerate(area_by_quadrant):
        chosen_quadrants[branch_index] = priority_order[
            int(np.argmax(area[priority_order]))
        ]

    assigned_quadrants = np.asarray(
        [chosen_quadrants == index for index in range(4)],
        dtype=bool,
    )
    assigned_quadrants = np.broadcast_to(
        assigned_quadrants[:, :, np.newaxis],
        (4, n_branches, n_radii),
    )
    return assigned_quadrants.copy()


def optic_disc_center_xy(source_data, image_shape: tuple[int, int]) -> np.ndarray:
    center = np.asarray(source_data.optic_disc_center, dtype=float).reshape(-1)
    if center.size < 2 or not np.all(np.isfinite(center[:2])):
        return np.asarray(
            [image_shape[1] / 2.0, image_shape[0] / 2.0],
            dtype=float,
        )
    return center[:2].copy()


def normalize_spatial_frame(
    branch_label_map: np.ndarray,
    optic_disc_center: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    if not _EYEFLOW_SPATIAL_Y_INVERTED:
        return branch_label_map, optic_disc_center

    normalized_center = optic_disc_center.copy()
    normalized_center[1] = branch_label_map.shape[0] - 1 - normalized_center[1]
    return np.flip(branch_label_map, axis=0).copy(), normalized_center
