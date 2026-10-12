"""Mask helpers used by retinal velocity estimation."""

from __future__ import annotations

import numpy as np

from calculations.topology import downsample_gap_zones


def gap_velocity_mask(
    analysis_zones: np.ndarray,
    spatial_shape: tuple[int, int],
) -> np.ndarray:
    """Prepare native gap support for a video without changing the persisted mask.

    A 512x512 mask becomes 256x256 using conservative 2x2 reduction. If its
    dimensions already match the video's spatial shape, return a copy.
    """
    zones = np.asarray(analysis_zones)
    if zones.ndim != 2 or zones.dtype != np.bool_:
        raise TypeError("analysis_zones must be a 2-D boolean array.")
    if zones.shape == tuple(spatial_shape):
        return zones.copy()
    if zones.shape != tuple(2 * size for size in spatial_shape):
        raise ValueError(
            "Gap-mask dimensions must match the video or be exactly twice its "
            f"spatial dimensions; got mask {zones.shape} and video {spatial_shape}."
        )
    return downsample_gap_zones(zones)


def elliptical_mask(ny: int, nx: int, radius_frac: float, center=None) -> np.ndarray:
    radius = max(0.0, min(1.0, float(radius_frac)))
    semi_x = max((nx / 2.0) * radius, np.finfo(float).eps)
    semi_y = max((ny / 2.0) * radius, np.finfo(float).eps)
    y_grid, x_grid = np.ogrid[:ny, :nx]
    cy, cx = (ny / 2.0, nx / 2.0) if center is None else center
    normalized = ((x_grid - cx) / semi_x) ** 2 + ((y_grid - cy) / semi_y) ** 2
    return normalized <= 1.0


def elliptical_annulus_mask(
    ny: int,
    nx: int,
    *,
    outer_radius_frac: float = 0.5,
    inner_radius_frac: float = 0.2,
) -> np.ndarray:
    outer = elliptical_mask(ny, nx, outer_radius_frac)
    inner = elliptical_mask(ny, nx, inner_radius_frac)
    return outer & ~inner
