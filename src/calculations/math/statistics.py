"""Statistical array helpers shared by scientific calculations."""

from __future__ import annotations

import numpy as np


def nanmean_float32(values, axis=None):
    array = np.asarray(values, dtype=np.float32)
    finite = np.isfinite(array)
    count = np.sum(finite, axis=axis)
    total = np.sum(np.where(finite, array, 0.0), axis=axis, dtype=np.float32)
    result = np.divide(
        total,
        count,
        out=np.full_like(total, np.nan, dtype=np.float32),
        where=count > 0,
    )
    if isinstance(result, np.ndarray):
        return result.astype(np.float32, copy=False)
    return np.float32(result)
