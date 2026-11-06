"""Per-beat interpolation for sampled vessel profiles."""

from __future__ import annotations

import numpy as np
from scipy.signal import resample

from calculations.math.cycle_boundaries import (
    normalize_cycle_boundaries,
)
from calculations.math import next_power_of_two


def interpolate_profiles_per_beat(
    profiles: np.ndarray,
    cycle_boundary_indexes,
    *,
    index_base: int = 0,
    valid_segments: np.ndarray | None = None,
) -> np.ndarray:
    """Interpolate profiles to ``(x, time, beat, branch, radius)``.

    Transform all spatial samples in one valid segment/beat together. Invalid
    dense ring/branch slots stay NaN without entering an FFT.
    """

    profiles = np.asarray(profiles, dtype=np.float32)
    if profiles.ndim != 4:
        raise ValueError(
            "profiles must have shape (radius, branch, frame, transverse_sample)."
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


def interpolate_segment_maps_per_beat(
    segment_maps: np.ndarray,
    cycle_boundary_indexes,
    *,
    index_base: int = 0,
    valid_segments: np.ndarray | None = None,
    array_allocator=None,
) -> np.ndarray:
    """Interpolate dense segment maps into the serialized per-beat layout.

    Scalar input is ``(radius, branch, frame, y, x)``. A final component axis
    is also supported. Output is ``(x, y, time, beat, branch, radius, ...)``.
    """

    maps = segment_maps if hasattr(segment_maps, "shape") else np.asarray(segment_maps, dtype=np.float32)
    if maps.ndim not in (5, 6):
        raise ValueError(
            "segment_maps must have shape (radius, branch, frame, y, x) "
            "with an optional final component axis."
        )
    radius_count, branch_count, frame_count, y_count, x_count = maps.shape[:5]
    if valid_segments is None:
        active = np.zeros((radius_count, branch_count), dtype=bool)
        for radius in range(radius_count):
            for branch in range(branch_count):
                for frame in range(frame_count):
                    if np.any(np.isfinite(maps[radius, branch, frame])):
                        active[radius, branch] = True
                        break
    else:
        active = np.asarray(valid_segments, dtype=bool)
        if active.shape != (radius_count, branch_count):
            raise ValueError(
                "valid_segments must match the map radius and branch axes."
            )

    boundaries = normalize_cycle_boundaries(
        cycle_boundary_indexes,
        frame_count,
        index_base=index_base,
    )
    beat_count = boundaries.size - 1
    time_count = next_power_of_two(int(np.max(np.diff(boundaries))))
    component_shape = maps.shape[5:]
    working_bytes = x_count * y_count * max(1, int(np.prod(component_shape))) * 4 * (
        2 * (int(np.max(np.diff(boundaries))) + 1) + 4 * (time_count + 1)
    )
    if working_bytes > 512 * 2**20:
        raise MemoryError("One displacement segment/beat FFT exceeds the 512 MiB working-memory allowance.")
    shape = (
            x_count,
            y_count,
            time_count,
            beat_count,
            branch_count,
            radius_count,
            *component_shape,
    )
    output = np.empty(shape, dtype=np.float32) if array_allocator is None else array_allocator(shape)
    output.fill(np.nan)

    for radius_index, branch_index in np.argwhere(active):
        for beat_index in range(beat_count):
            start = int(boundaries[beat_index])
            stop = int(boundaries[beat_index + 1]) + 1
            interpolated = _resample_map_axis0(
                maps[radius_index, branch_index, start:stop],
                time_count + 1,
            )[:-1]
            axes = (2, 1, 0, *range(3, interpolated.ndim))
            output[
                :,
                :,
                :,
                beat_index,
                branch_index,
                radius_index,
                ...,
            ] = interpolated.transpose(axes)

    return output


def _resample_map_axis0(values: np.ndarray, target_length: int) -> np.ndarray:
    source = np.asarray(values, dtype=np.float32)
    source_length = int(source.shape[0])
    flattened = source.reshape(source_length, -1)
    active = np.any(np.isfinite(flattened), axis=0)
    interpolated = np.full(
        (target_length, flattened.shape[1]),
        np.nan,
        dtype=np.float32,
    )
    if np.any(active):
        interpolated[:, active] = resample(
            flattened[:, active],
            target_length,
            axis=0,
        ).astype(np.float32, copy=False)
    return interpolated.reshape(target_length, *source.shape[1:])


__all__ = [
    "interpolate_profiles_per_beat",
    "interpolate_segment_maps_per_beat",
]
