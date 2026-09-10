"""Per-beat interpolation for cross-section profiles."""

from __future__ import annotations

import numpy as np
from scipy.signal import resample

from calculations.blood_flow_velocity.signal_analysis.per_beat._signal_utils import (
    normalize_cycle_boundaries,
)
from calculations.math import next_power_of_two


def interpolate_velocity_profiles_per_beat(
    velocity_profiles: np.ndarray,
    cycle_boundary_indexes,
    *,
    index_base: int = 0,
    valid_segments: np.ndarray | None = None,
) -> np.ndarray:
    """Interpolate profiles to ``(x, time, beat, branch, radius)``.

    Transform all spatial samples in one valid segment/beat together. Invalid
    dense ring/branch slots stay NaN without entering an FFT.
    """

    profiles = np.asarray(velocity_profiles, dtype=np.float32)
    if profiles.ndim != 4:
        raise ValueError(
            "velocity_profiles must have shape "
            "(radius, branch, frame, transverse_sample)."
        )

    radius_count, branch_count, frame_count, x_count = profiles.shape
    if valid_segments is None:
        active = np.any(np.isfinite(profiles), axis=(2, 3))
    else:
        active = np.asarray(valid_segments, dtype=bool)
        if active.shape != (radius_count, branch_count):
            raise ValueError(
                "valid_segments must match the profile radius and branch axes."
            )
    boundaries = normalize_cycle_boundaries(
        cycle_boundary_indexes,
        frame_count,
        index_base=index_base,
    )
    beat_count = boundaries.size - 1
    time_count = next_power_of_two(int(np.max(np.diff(boundaries))))
    output = np.full(
        (x_count, time_count, beat_count, branch_count, radius_count),
        np.nan,
        dtype=np.float32,
    )

    for radius_index, branch_index in np.argwhere(active):
        for beat_index in range(beat_count):
            start = int(boundaries[beat_index])
            stop = int(boundaries[beat_index + 1]) + 1
            beat_profiles = profiles[
                radius_index,
                branch_index,
                start:stop,
                :,
            ]
            output[
                :,
                :,
                beat_index,
                branch_index,
                radius_index,
            ] = resample(
                beat_profiles,
                time_count + 1,
                axis=0,
            )[:-1].T.astype(np.float32, copy=False)

    return output


__all__ = ["interpolate_velocity_profiles_per_beat"]
