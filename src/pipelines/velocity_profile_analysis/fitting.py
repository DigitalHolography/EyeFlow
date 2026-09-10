"""Independent weighted quadratic fits to (x, time, beat, branch, radius)."""

from __future__ import annotations

import numpy as np

FLOAT_OUTPUTS = (
    "a",
    "b",
    "c",
    "fit_r_squared",
    "fit_rmse",
    "fit_rss",
    "fit_weighted_r_squared",
    "fit_weighted_rmse",
    "fit_weighted_rss",
    "index_center",
    "index_left_zero",
    "index_right_zero",
    "Qv_fit",
    "Qv",
)
COUNT_OUTPUTS = ("n_fit_samples", "n_area_samples")
DEFAULT_TIME_BLOCK_SIZE = 256


def border_weights(sample_count: int) -> np.ndarray:
    """Half weight in the outer quarters of the original index domain."""
    weights = np.ones(sample_count, dtype=np.float64)
    if sample_count > 1:
        u = np.arange(sample_count, dtype=np.float64) / (sample_count - 1)
        weights[(u < 0.25) | (u > 0.75)] = 0.5
    return weights


def _allocate(shape, *, dtype=np.float64):
    outputs = {name: np.full(shape, np.nan, dtype=dtype) for name in FLOAT_OUTPUTS}
    outputs.update({name: np.zeros(shape, dtype=np.int32) for name in COUNT_OUTPUTS})
    return outputs


def analyze_velocity_profiles(v, *, time_block_size=DEFAULT_TIME_BLOCK_SIZE):
    """Read bounded time slabs from a NumPy array or HDF5 dataset.

    Return float32 measurements and int32 counts with (time, beat, branch,
    radius) axes. Neither the entire input nor a full fitted movie is copied.
    """
    shape = getattr(v, "shape", None)
    if shape is None or len(shape) != 5:
        raise ValueError("v must have shape (x, time, beat, branch, radius).")
    if np.dtype(v.dtype).kind not in "biuf":
        raise ValueError("v must contain real numeric velocity samples.")
    if not isinstance(time_block_size, (int, np.integer)) or time_block_size < 1:
        raise ValueError("time_block_size must be a positive integer.")
    nx, nt, nb, nk, nr = shape
    outputs = _allocate((nt, nb, nk, nr), dtype=np.float32)
    x = np.arange(nx, dtype=np.float64)
    weights = border_weights(nx)
    midpoint = (nx - 1) / 2.0
    scale = max(midpoint, 1.0)
    z = (x - midpoint) / scale
    design = np.column_stack((z * z, z, np.ones(nx)))
    for beat, branch, radius in np.ndindex(nb, nk, nr):
        for start in range(0, nt, time_block_size):
            stop = min(start + time_block_size, nt)
            values = np.asarray(v[:, start:stop, beat, branch, radius], dtype=np.float64)
            block = _fit_block(values, x, weights, design, midpoint, scale)
            for name, result in block.items():
                outputs[name][start:stop, beat, branch, radius] = result
    return outputs


def _fit_block(values, x, weights, design, midpoint, scale):
    result = _allocate((values.shape[1],))
    valid = np.isfinite(values)
    result["n_fit_samples"][:] = valid.sum(axis=0)
    if values.shape[0] < 3:
        return result
    # Each distinct finite mask needs one factorization for all its time samples.
    _, groups = np.unique(np.packbits(valid.T, axis=1), axis=0, return_inverse=True)
    for group in range(int(groups.max()) + 1):
        columns = np.flatnonzero(groups == group)
        finite = valid[:, columns[0]]
        n = int(finite.sum())
        if n < 3:
            continue
        y = values[np.ix_(finite, columns)]
        w = weights[finite]
        sqrt_w = np.sqrt(w)
        matrix = design[finite]
        try:
            coefficients, _, rank, _ = np.linalg.lstsq(
                matrix * sqrt_w[:, None],
                y * sqrt_w[:, None],
                rcond=None,
            )
        except np.linalg.LinAlgError:
            continue
        if rank < 3:
            continue
        alpha, beta, gamma = coefficients
        converted = np.array(
            (
                alpha / scale**2,
                beta / scale - 2 * alpha * midpoint / scale**2,
                gamma - beta * midpoint / scale + alpha * (midpoint / scale) ** 2,
            )
        )
        good = np.all(np.isfinite(coefficients), axis=0) & np.all(np.isfinite(converted), axis=0)
        if not np.any(good):
            continue
        columns, y = columns[good], y[:, good]
        coefficients, converted = coefficients[:, good], converted[:, good]
        for name, row in zip(("a", "b", "c"), converted, strict=True):
            result[name][columns] = row
        fitted = matrix @ coefficients
        residual_sq = (y - fitted) ** 2
        rss = residual_sq.sum(axis=0)
        wrss = (w[:, None] * residual_sq).sum(axis=0)
        sst = ((y - y.mean(axis=0)) ** 2).sum(axis=0)
        wmean = (w[:, None] * y).sum(axis=0) / w.sum()
        wsst = (w[:, None] * (y - wmean) ** 2).sum(axis=0)
        result["fit_rss"][columns] = rss
        result["fit_rmse"][columns] = np.sqrt(rss / n)
        result["fit_weighted_rss"][columns] = wrss
        result["fit_weighted_rmse"][columns] = np.sqrt(wrss / w.sum())
        # An exactly constant observed profile has no defined R-squared, even
        # when reduction roundoff would produce a tiny nonzero SST.
        nonconstant = np.any(y != y[0], axis=0)
        ok = (sst > 0) & nonconstant
        result["fit_r_squared"][columns[ok]] = 1 - rss[ok] / sst[ok]
        ok = (wsst > 0) & nonconstant
        result["fit_weighted_r_squared"][columns[ok]] = 1 - wrss[ok] / wsst[ok]
        _geometry_and_areas(result, columns, coefficients, x[finite], y, fitted, midpoint, scale)
    return result


def _geometry_and_areas(result, columns, coefficients, observed_x, y, fitted, midpoint, scale):
    alpha, beta, gamma = coefficients
    tolerance = 64 * np.finfo(np.float64).eps * np.maximum(1.0, np.max(np.abs(y), axis=0))
    for j in np.flatnonzero(alpha < -tolerance):
        column = columns[j]
        a, b, c = alpha[j], beta[j], gamma[j]
        center = midpoint - scale * b / (2 * a)
        if not np.isfinite(center):
            continue
        result["index_center"][column] = center
        disc = b * b - 4 * a * c
        # Roundoff around a repeated root must not create two artificial zeros.
        disc_tolerance = 64 * np.finfo(np.float64).eps * max(1.0, a * a, b * b, abs(4 * a * c))
        if not np.isfinite(disc) or disc <= disc_tolerance:
            continue
        q = -0.5 * (b + np.copysign(np.sqrt(disc), b))
        roots = np.sort(midpoint + scale * np.array((q / a, c / q)))
        if not np.all(np.isfinite(roots)) or roots[0] >= roots[1]:
            continue
        left, right = roots
        result["index_left_zero"][column] = left
        result["index_right_zero"][column] = right
        support = (observed_x >= left) & (observed_x <= right)
        count = np.count_nonzero(support)
        result["n_area_samples"][column] = count
        if count:
            result["Qv_fit"][column] = fitted[support, j].sum()
            result["Qv"][column] = y[support, j].sum()
