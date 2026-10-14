"""Mean-preserving SVD filtering of sample-by-channel matrices."""

from __future__ import annotations

from dataclasses import dataclass
import operator

import numpy as np


@dataclass(frozen=True)
class SvdFilterResult:
    filtered: np.ndarray
    singular_values: np.ndarray
    temporal_mean: np.ndarray
    valid_rows: np.ndarray
    valid_columns: np.ndarray
    modes_removed: int
    applied: bool


def remove_leading_svd_modes(matrix, *, mode_count: int = 2) -> SvdFilterResult:
    """Remove dominant centered modes, then restore each channel's mean.

    Entirely invalid columns are excluded. Fit only rows finite across all
    remaining columns; excluded rows and existing NaNs are left unchanged.
    Retain at least one available mode. Insufficient fit dimensions return
    unchanged values with applied=False instead of erasing all variation.
    """
    if isinstance(mode_count, (bool, np.bool_)):
        raise ValueError("mode_count must be a non-negative integer.")
    try:
        mode_count = operator.index(mode_count)
    except TypeError as error:
        raise ValueError("mode_count must be a non-negative integer.") from error
    if mode_count < 0:
        raise ValueError("mode_count must be a non-negative integer.")
    values = np.asarray(matrix, dtype=np.float64).copy()
    if values.ndim != 2 or any(size == 0 for size in values.shape):
        raise ValueError("SVD filtering requires a nonempty (sample, channel) matrix.")
    finite = np.isfinite(values)
    values[~finite] = np.nan
    columns = np.any(finite, axis=0)
    rows = np.all(finite[:, columns], axis=1) if columns.any() else np.zeros(values.shape[0], bool)
    means = np.full(values.shape[1], np.nan, dtype=np.float64)
    fit = values[np.ix_(rows, columns)]
    if fit.shape[0]:
        means[columns] = np.mean(fit, axis=0)
    if min(fit.shape) <= mode_count:
        return SvdFilterResult(values.astype(np.float32), np.empty(0), means, rows, columns, 0, False)
    centered = fit - means[columns]
    temporal_modes, singular_values, channel_modes = np.linalg.svd(centered, full_matrices=False)
    kept = singular_values.copy()
    kept[:mode_count] = 0.0
    reconstructed = (temporal_modes * kept) @ channel_modes
    values[np.ix_(rows, columns)] = reconstructed + means[columns]
    return SvdFilterResult(
        values.astype(np.float32), singular_values, means, rows, columns, mode_count, True,
    )
