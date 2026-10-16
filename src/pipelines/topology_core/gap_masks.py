"""Intensity-derived gap support in the native HD (y, x) image frame."""

from __future__ import annotations

import warnings

import numpy as np
from scipy import ndimage as ndi
from skimage.exposure import equalize_adapthist
from skimage.filters import threshold_otsu

from calculations.topology.capillary import (
    GAPS_VESSEL_DILATION_RADIUS_PIXELS as DILATION_RADIUS_PIXELS,
)

CONTRAST_PERCENTILES = (1.0, 99.0)
CLAHE_CLIP_LIMIT = 0.01
CLAHE_NBINS = 256
SMOOTHING_SIGMA_PIXELS = 1.0
_MEDIAN_TILE_BYTES = 8 * 1024 * 1024


def temporal_median_image(volume) -> np.ndarray:
    """Reduce frames exactly, reading spatial tiles rather than the whole movie.

    Nonfinite samples are ignored; pixels without finite samples remain NaN.
    Spatial axes are retained verbatim, including for square images.
    """
    frames, height, width = volume.shape
    if np.dtype(volume.dtype).kind not in "iuf":
        raise TypeError("moment0ff must contain real numeric intensities.")
    if min(frames, height, width) < 1:
        raise ValueError("moment0ff must have nonempty (frame, y, x) axes.")
    tile_side = max(1, int(np.sqrt(_MEDIAN_TILE_BYTES / (frames * 8))))
    result = np.empty((height, width), dtype=np.float32)
    for y in range(0, height, tile_side):
        for x in range(0, width, tile_side):
            tile = np.array(
                volume[:, y:y + tile_side, x:x + tile_side], dtype=np.float64,
                copy=True,
            )
            tile[~np.isfinite(tile)] = np.nan
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                result[y:y + tile_side, x:x + tile_side] = np.nanmedian(tile, axis=0)
    return result


def gaps_from_m0ff(volume) -> tuple[np.ndarray, dict[str, object]]:
    """Enhance and threshold vessels, dilate with a 7x7 square, then invert."""
    image = temporal_median_image(volume)
    valid = np.isfinite(image)
    attrs = {
        "temporal_aggregation": "median",
        "nonfinite_policy": "ignore temporal nonfinite samples; invalid pixels excluded",
        "contrast_percentiles": CONTRAST_PERCENTILES,
        "contrast_method": "percentile_stretch_then_CLAHE",
        "clahe_clip_limit": CLAHE_CLIP_LIMIT,
        "clahe_nbins": CLAHE_NBINS,
        "clahe_kernel_size": [max(1, size // 8) for size in image.shape],
        "smoothing_method": "Gaussian, finite-support normalized",
        "smoothing_sigma_pixels": SMOOTHING_SIGMA_PIXELS,
        "smoothing_border_mode": "reflect",
        "binarization": "valid AND NOT dilate(smoothed_intensity >= global_Otsu_threshold)",
        "gap_polarity": "dark",
        "dilation_radius_pixels": DILATION_RADIUS_PIXELS,
        "dilation_footprint": "7x7_square",
        "dilation_distance_metric": "chebyshev",
        "dilation_iterations": 1,
        "dilation_applies_to": "bright_vessels",
        "dilation_border_value": 0,
        "anatomical_exclusions": "none",
        "valid_pixel_count": int(valid.sum()),
        "otsu_threshold": float("nan"),
        "contrast_low": float("nan"),
        "contrast_high": float("nan"),
    }
    empty = np.zeros(image.shape, dtype=bool)
    if not valid.any():
        return empty, attrs
    low, high = np.percentile(image[valid], CONTRAST_PERCENTILES)
    # Sparse dark/bright features may lie outside both percentiles.
    if high <= low:
        low, high = float(image[valid].min()), float(image[valid].max())
    attrs.update(contrast_low=float(low), contrast_high=float(high))
    if high <= low:
        return empty, attrs  # A uniform field contains no distinguishable gaps.
    normalized = np.clip((image.astype(np.float64) - low) / (high - low), 0, 1)
    normalized[~valid] = np.median(normalized[valid])
    enhanced = equalize_adapthist(
        normalized, kernel_size=attrs["clahe_kernel_size"],
        clip_limit=CLAHE_CLIP_LIMIT, nbins=CLAHE_NBINS,
    )
    weights = ndi.gaussian_filter(valid.astype(float), SMOOTHING_SIGMA_PIXELS)
    smoothed = ndi.gaussian_filter(np.where(valid, enhanced, 0), SMOOTHING_SIGMA_PIXELS)
    np.divide(smoothed, weights, out=smoothed, where=weights > 0)
    if np.ptp(smoothed[valid]) == 0:
        return empty, attrs
    threshold = float(threshold_otsu(smoothed[valid]))
    attrs["otsu_threshold"] = threshold
    vessels = valid & (smoothed >= threshold)
    footprint_width = 2 * DILATION_RADIUS_PIXELS + 1
    vessels = ndi.binary_dilation(
        vessels, structure=np.ones((footprint_width, footprint_width), dtype=bool),
        iterations=1, border_value=0,
    )
    return valid & ~vessels, attrs
