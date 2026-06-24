"""Image and array transformation helpers shared by scientific calculations."""

from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi


def rotate_image_with_nan(image: np.ndarray, angle: float) -> np.ndarray:
    valid = np.isfinite(image)
    fill = np.float32(np.nanmin(image[valid]) - 1.0) if np.any(valid) else np.float32(0.0)
    filled = np.where(valid, image, fill)
    rotated = ndi.rotate(filled, angle, reshape=False, order=1, mode="constant", cval=float(fill))
    mask = ndi.rotate(
        valid.astype(np.float32),
        angle,
        reshape=False,
        order=1,
        mode="constant",
        cval=0.0,
    )
    rotated[mask < 0.5] = np.nan
    return rotated.astype(np.float32, copy=False)


def rotate_array_threshold(array: np.ndarray, angle: float, threshold: float = 0.5) -> np.ndarray:
    rotated = ndi.rotate(
        np.asarray(array, dtype=np.float32),
        angle,
        reshape=False,
        order=1,
        mode="constant",
        cval=0.0,
    )
    return rotated >= np.float32(threshold)
