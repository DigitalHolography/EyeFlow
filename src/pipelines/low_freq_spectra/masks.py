"""Adapt native gap support to the HD spectral video's spatial grid."""

from __future__ import annotations

import numpy as np

from calculations.topology import downsample_gap_zones


def gap_velocity_mask(
    analysis_zones: np.ndarray,
    spatial_shape: tuple[int, int],
) -> np.ndarray:
    """Interpolate Boolean support with nearest neighbors; preserve native masks."""
    return downsample_gap_zones(analysis_zones, target_shape=spatial_shape)
