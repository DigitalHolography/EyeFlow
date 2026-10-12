"""Native analysis support outside retinal vessels, optic disc, and choroid."""

from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi

GAPS_CUTOFF_DIAMETER_WIDTH_FRACTION = 0.8
GAPS_VESSEL_DILATION_RADIUS_PIXELS = 3


def capillary_analysis_zones(
    retinal_vessel_mask: np.ndarray,
    optic_disc_mask: np.ndarray,
    choroid_vessel_mask: np.ndarray,
) -> np.ndarray:
    """Return central gaps outside dilated vessels and the optic disc.

    Inputs and output use native ``(y, x)`` coordinates. A true output pixel
    belongs to the capillary analysis support. Retinal and choroid vessels
    are dilated by three pixels using one pass of a 7x7 square footprint.
    The optic disc is excluded without dilation. Support is restricted to an image-centered
    circle with diameter 80% of the image width. Inputs are never mutated.
    """
    masks = []
    for name, mask in (
        ("retinal_vessel_mask", retinal_vessel_mask),
        ("optic_disc_mask", optic_disc_mask),
        ("choroid_vessel_mask", choroid_vessel_mask),
    ):
        value = np.asarray(mask)
        if value.ndim != 2 or value.dtype != np.bool_:
            raise TypeError(f"{name} must be a 2-D boolean array.")
        if masks and value.shape != masks[0].shape:
            raise ValueError(
                f"{name} must share the aligned mask shape {masks[0].shape}, "
                f"got {value.shape}."
            )
        masks.append(value)
    ny, nx = masks[0].shape
    if ny < 1 or nx < 1:
        raise ValueError("Gap masks must have positive image dimensions.")
    footprint_width = 2 * GAPS_VESSEL_DILATION_RADIUS_PIXELS + 1
    vessels = ndi.binary_dilation(
        masks[0] | masks[2],
        structure=np.ones((footprint_width, footprint_width), dtype=bool),
        iterations=1,
        border_value=0,
    )
    y, x = np.ogrid[:ny, :nx]
    radius = GAPS_CUTOFF_DIAMETER_WIDTH_FRACTION * nx / 2.0
    cutoff_disc = (x - nx / 2.0) ** 2 + (y - ny / 2.0) ** 2 <= radius**2
    return cutoff_disc & ~(vessels | masks[1])


def downsample_gap_zones(analysis_zones: np.ndarray) -> np.ndarray:
    """Reduce native gap support by two, retaining only wholly included 2x2 blocks.

    One excluded pixel excludes the entire reduced pixel. The input is never
    mutated; this reduction is for velocity-mask application, not persistence.
    Axes stay in native ``(y, x)`` order.
    """
    zones = np.asarray(analysis_zones)
    if zones.ndim != 2 or zones.dtype != np.bool_:
        raise TypeError("analysis_zones must be a 2-D boolean array.")
    ny, nx = zones.shape
    if ny < 1 or nx < 1 or ny % 2 or nx % 2:
        raise ValueError("Gap-mask downsampling requires positive, even spatial dimensions.")
    return zones.reshape(ny // 2, 2, nx // 2, 2).all(axis=(1, 3))
