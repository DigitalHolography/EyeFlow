"""Write PNG artifacts for EyeFlow runs."""

from pathlib import Path

import numpy as np
from skimage.io import imsave


def write_png_file(path: str | Path, image) -> Path:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    imsave(target, _uint8_image(image), check_contrast=False)
    return target


def _uint8_image(image) -> np.ndarray:
    array = np.asarray(image)
    if array.dtype == np.uint8:
        return array
    if array.dtype == bool:
        return (array.astype(np.uint8) * 255)
    if np.issubdtype(array.dtype, np.floating):
        return _normalize_float(array)
    clipped = np.clip(array, 0, 255)
    return clipped.astype(np.uint8, copy=False)


def _normalize_float(array: np.ndarray) -> np.ndarray:
    finite = np.isfinite(array)
    if not np.any(finite):
        return np.zeros(array.shape, dtype=np.uint8)
    values = array[finite]
    min_value = np.min(values)
    span = np.max(values) - min_value
    if span <= 0:
        return np.where(finite, 255, 0).astype(np.uint8)
    scaled = np.zeros(array.shape, dtype=np.float32)
    scaled[finite] = (array[finite] - min_value) / span
    return np.rint(scaled * 255).astype(np.uint8)
