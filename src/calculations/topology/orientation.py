"""Determine vessel direction from native centerline geometry."""

from __future__ import annotations

import numpy as np

from .geometry import image_half_diagonal
from .segments import SegmentTopology


def determine_segment_rotations(topology: SegmentTopology) -> np.ndarray:
    """Return the rotation that makes every valid vessel segment upright.

    The returned array has ``(annulus, branch)`` shape. Directions are measured
    along the existing centerline and point away from the optic disc. This
    preserves a consistent sign for later vector-basis correction.
    """

    rotations = np.full(topology.valid_segments.shape, np.nan, dtype=np.float32)
    for annulus_index, branch_index in np.argwhere(topology.valid_segments):
        segment = topology.annulus_masks[annulus_index] & (
            topology.labels == int(topology.branch_ids[branch_index])
        )
        segment_centerline = topology.centerline & segment
        tilt = _segment_tilt(
            segment_centerline,
            topology.optic_disc_center_xy,
        )
        if np.isfinite(tilt):
            rotations[annulus_index, branch_index] = np.float32(tilt + 90.0)
    return rotations


def _segment_tilt(
    segment_centerline: np.ndarray,
    optic_disc_center_xy: tuple[float, float],
) -> float:
    point_y, point_x = np.nonzero(segment_centerline)
    if point_x.size < 2:
        return float("nan")

    points = np.column_stack((point_x, point_y)).astype(np.float64)
    center = np.median(points, axis=0)
    centered = points - center
    moment = centered.T @ centered
    eigenvalues, eigenvectors = np.linalg.eigh(moment)
    if eigenvalues[-1] <= eigenvalues[0]:
        return float("nan")

    axis = eigenvectors[:, -1]
    projections = centered @ axis
    if np.ptp(projections) < 1.0:
        return float("nan")

    radii = _radius_grid(segment_centerline.shape, optic_disc_center_xy)[
        point_y,
        point_x,
    ]
    lower_projection, upper_projection = np.quantile(
        projections,
        (0.25, 0.75),
    )
    lower_radius = float(
        np.median(radii[projections <= lower_projection])
    )
    upper_radius = float(
        np.median(radii[projections >= upper_projection])
    )
    if upper_radius < lower_radius:
        axis = -axis
    elif np.isclose(upper_radius, lower_radius):
        optic_center_x, optic_center_y = optic_disc_center_xy
        radial_direction = np.asarray(
            (center[0] - optic_center_x, center[1] - optic_center_y),
            dtype=np.float64,
        )
        radial_alignment = float(np.dot(axis, radial_direction))
        if np.isclose(radial_alignment, 0.0):
            return float("nan")
        if radial_alignment < 0.0:
            axis = -axis
    return float(np.degrees(np.arctan2(axis[1], axis[0])))


def _radius_grid(
    image_shape: tuple[int, int],
    optic_disc_center_xy: tuple[float, float],
) -> np.ndarray:
    ny, nx = image_shape
    center_x, center_y = optic_disc_center_xy
    scale = np.float32(1.0 / max(image_half_diagonal(ny, nx), 1.0))
    y = (np.arange(ny, dtype=np.float32)[:, None] - np.float32(center_y)) * scale
    x = (np.arange(nx, dtype=np.float32)[None, :] - np.float32(center_x)) * scale
    return np.sqrt(x**2 + y**2)
