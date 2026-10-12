"""Native analysis support outside retinal vessels, optic disc, and choroid."""

from __future__ import annotations

import numpy as np


def capillary_analysis_zones(
    retinal_vessel_mask: np.ndarray,
    optic_disc_mask: np.ndarray,
    choroid_vessel_mask: np.ndarray,
) -> np.ndarray:
    """Return the complement of three aligned, binarized exclusion masks.

    Inputs and output use native ``(y, x)`` coordinates. A true output pixel
    belongs to the capillary analysis support. Inputs are never mutated.
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
    return ~(masks[0] | masks[1] | masks[2])
