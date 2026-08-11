"""Image and array transformation helpers shared by scientific calculations."""

from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi


def rotate_image_with_nan(image: np.ndarray, angle: float) -> np.ndarray:
    """Rotate an image without mixing missing pixels into finite values.

    Linear interpolation is applied separately to zero-filled values and to
    their validity weights.  Dividing the two rotated arrays performs a
    normalized convolution: a finite output depends only on finite input
    pixels, so an arbitrary NaN fill value cannot create out-of-range values.
    """

    valid = np.isfinite(image)
    values = ndi.rotate(
        np.where(valid, image, np.float32(0.0)),
        angle,
        reshape=False,
        order=1,
        mode="constant",
        cval=0.0,
    )
    weights = ndi.rotate(
        valid.astype(np.float32),
        angle,
        reshape=False,
        order=1,
        mode="constant",
        cval=0.0,
    )
    rotated = np.full(values.shape, np.nan, dtype=np.float32)
    np.divide(
        values,
        weights,
        out=rotated,
        where=weights >= np.float32(0.5),
    )
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
