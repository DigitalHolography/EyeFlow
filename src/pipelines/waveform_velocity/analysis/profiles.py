"""Numerical analysis of waveform-velocity cross-section profiles."""

from __future__ import annotations

import numpy as np

from calculations.math import nanmean_float32
from calculations.topology import dilate_segment_masks

DEFAULT_PROFILE_MASK_DILATION_ITERATIONS = 10


def velocity_fft_transverse_profiles(
    velocity_maps_per_beat: np.ndarray,
    segment_masks: np.ndarray,
    *,
    mask_dilation_pixels: int = DEFAULT_PROFILE_MASK_DILATION_ITERATIONS,
) -> tuple[np.ndarray, np.ndarray]:
    """Return FFT profiles shaped ``(x, frequency, beat, branch, radius)``.

    ``velocity_maps_per_beat`` must already have shape
    ``(x, y, time, beat, branch, radius)``. The FFT is applied along its time
    axis independently for every pixel, beat, branch, and radius.
    """

    maps = np.asarray(velocity_maps_per_beat, dtype=np.float32)
    masks = np.asarray(segment_masks, dtype=bool)
    if maps.ndim != 6:
        raise ValueError(
            "velocity_maps_per_beat must have shape "
            "(x, y, time, beat, branch, radius)."
        )
    expected_mask_shape = (
        maps.shape[5],
        maps.shape[4],
        maps.shape[1],
        maps.shape[0],
    )
    if masks.shape != expected_mask_shape:
        raise ValueError(
            "segment_masks must have shape (radius, branch, y, x) matching "
            "velocity_maps_per_beat."
        )

    dilated_masks = dilate_segment_masks(
        masks,
        iterations=mask_dilation_pixels,
        horizontal_only=True,
    )
    output_shape = (
        maps.shape[0],
        maps.shape[2],
        maps.shape[3],
        maps.shape[4],
        maps.shape[5],
    )
    unmasked = np.full(output_shape, np.nan, dtype=np.float32)
    masked = np.full(output_shape, np.nan, dtype=np.float32)
    for radius_index in range(maps.shape[5]):
        for branch_index in range(maps.shape[4]):
            magnitude = np.abs(
                np.fft.fft(
                    maps[..., branch_index, radius_index],
                    axis=2,
                )
            ).astype(np.float32, copy=False)
            unmasked[..., branch_index, radius_index] = nanmean_float32(
                magnitude,
                axis=1,
            )
            xy_mask = dilated_masks[radius_index, branch_index].T
            masked[..., branch_index, radius_index] = nanmean_float32(
                np.where(
                    xy_mask[:, :, None, None],
                    magnitude,
                    np.float32(np.nan),
                ),
                axis=1,
            )
    return unmasked, masked


__all__ = [
    "DEFAULT_PROFILE_MASK_DILATION_ITERATIONS",
    "velocity_fft_transverse_profiles",
]
