"""Array cleanup and normalization helpers shared by scientific calculations."""

from __future__ import annotations

import numpy as np


def as_float32_vector(values) -> np.ndarray:
    if values is None:
        return np.asarray([], dtype=np.float32)
    return np.asarray(values, dtype=np.float32).reshape(-1)


def as_nonnegative_int_indexes(values) -> np.ndarray:
    if values is None:
        return np.asarray([], dtype=np.int32)
    indexes = np.asarray(values, dtype=np.int32).reshape(-1)
    return indexes[indexes >= 0]


def finite_image(image: np.ndarray) -> np.ndarray:
    data = np.asarray(image, dtype=np.float32)
    finite = np.isfinite(data)
    if not np.any(finite):
        return np.zeros_like(data)
    clean = data.copy()
    clean[~finite] = float(np.nanmin(data[finite]))
    return clean


def nan_to_mean(values: np.ndarray) -> np.ndarray:
    clean = np.asarray(values, dtype=np.float32).copy()
    finite = np.isfinite(clean)
    if not np.any(finite):
        return np.zeros_like(clean)
    clean[~finite] = float(np.nanmean(clean[finite]))
    return clean


def standardize(values: np.ndarray) -> np.ndarray:
    clean = nan_to_mean(values)
    std = float(np.nanstd(clean))
    if std <= 0:
        return np.zeros_like(clean)
    return ((clean - float(np.nanmean(clean))) / std).astype(np.float32)


def rescale(values: np.ndarray) -> np.ndarray:
    data = finite_image(values).astype(np.float32)
    vmin = float(np.nanmin(data))
    vmax = float(np.nanmax(data))
    if vmax <= vmin:
        return np.zeros_like(data)
    return (data - vmin) / (vmax - vmin)
