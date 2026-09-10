"""Independent weighted quadratic fits to per-beat velocity profiles."""

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
    """Give half weight to samples in the outer quarters of the index domain."""

    weights = np.ones(sample_count, dtype=np.float64)
    if sample_count > 1:
        normalized = np.arange(sample_count, dtype=np.float64) / (sample_count - 1)
        weights[(normalized < 0.25) | (normalized > 0.75)] = 0.5
    return weights


def _allocate(shape, *, dtype=np.float64):
    outputs = {name: np.full(shape, np.nan, dtype=dtype) for name in FLOAT_OUTPUTS}
    outputs.update({name: np.zeros(shape, dtype=np.int32) for name in COUNT_OUTPUTS})
    return outputs


def analyze_velocity_profiles(values, *, time_block_size=DEFAULT_TIME_BLOCK_SIZE):
    """Analyze ``(x, time, beat, branch, radius)`` data in bounded time slabs."""

    shape = getattr(values, "shape", None)
    if shape is None or len(shape) != 5:
        raise ValueError("values must have shape (x, time, beat, branch, radius).")
    if np.dtype(values.dtype).kind not in "biuf":
        raise ValueError("values must contain real numeric velocity samples.")
    if not isinstance(time_block_size, (int, np.integer)) or time_block_size < 1:
        raise ValueError("time_block_size must be a positive integer.")

    sample_count, time_count, beat_count, branch_count, radius_count = shape
    outputs = _allocate(
        (time_count, beat_count, branch_count, radius_count),
        dtype=np.float32,
    )
    x = np.arange(sample_count, dtype=np.float64)
    weights = border_weights(sample_count)
    midpoint = (sample_count - 1) / 2.0
    scale = max(midpoint, 1.0)
    normalized_x = (x - midpoint) / scale
    design = np.column_stack(
        (normalized_x * normalized_x, normalized_x, np.ones(sample_count))
    )

    for beat, branch, radius in np.ndindex(
        beat_count,
        branch_count,
        radius_count,
    ):
        for start in range(0, time_count, time_block_size):
            stop = min(start + time_block_size, time_count)
            block_values = np.asarray(
                values[:, start:stop, beat, branch, radius],
                dtype=np.float64,
            )
            block = _fit_block(
                block_values,
                x,
                weights,
                design,
                midpoint,
                scale,
            )
            for name, result in block.items():
                outputs[name][start:stop, beat, branch, radius] = result
    return outputs


def _fit_block(values, x, weights, design, midpoint, scale):
    result = _allocate((values.shape[1],))
    valid = np.isfinite(values)
    result["n_fit_samples"][:] = valid.sum(axis=0)
    if values.shape[0] < 3:
        return result

    # Reuse one factorization for all time samples sharing the same finite mask.
    _, groups = np.unique(np.packbits(valid.T, axis=1), axis=0, return_inverse=True)
    for group in range(int(groups.max()) + 1):
        columns = np.flatnonzero(groups == group)
        finite = valid[:, columns[0]]
        sample_count = int(finite.sum())
        if sample_count < 3:
            continue
        observed = values[np.ix_(finite, columns)]
        finite_weights = weights[finite]
        square_root_weights = np.sqrt(finite_weights)
        matrix = design[finite]
        try:
            coefficients, _, rank, _ = np.linalg.lstsq(
                matrix * square_root_weights[:, None],
                observed * square_root_weights[:, None],
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
                gamma
                - beta * midpoint / scale
                + alpha * (midpoint / scale) ** 2,
            )
        )
        good = np.all(np.isfinite(coefficients), axis=0) & np.all(
            np.isfinite(converted),
            axis=0,
        )
        if not np.any(good):
            continue
        columns = columns[good]
        observed = observed[:, good]
        coefficients = coefficients[:, good]
        converted = converted[:, good]
        for name, row in zip(("a", "b", "c"), converted, strict=True):
            result[name][columns] = row

        fitted = matrix @ coefficients
        residual_squared = (observed - fitted) ** 2
        rss = residual_squared.sum(axis=0)
        weighted_rss = (finite_weights[:, None] * residual_squared).sum(axis=0)
        sst = ((observed - observed.mean(axis=0)) ** 2).sum(axis=0)
        weighted_mean = (
            (finite_weights[:, None] * observed).sum(axis=0)
            / finite_weights.sum()
        )
        weighted_sst = (
            finite_weights[:, None] * (observed - weighted_mean) ** 2
        ).sum(axis=0)
        result["fit_rss"][columns] = rss
        result["fit_rmse"][columns] = np.sqrt(rss / sample_count)
        result["fit_weighted_rss"][columns] = weighted_rss
        result["fit_weighted_rmse"][columns] = np.sqrt(
            weighted_rss / finite_weights.sum()
        )
        nonconstant = np.any(observed != observed[0], axis=0)
        valid_r_squared = (sst > 0) & nonconstant
        result["fit_r_squared"][columns[valid_r_squared]] = (
            1 - rss[valid_r_squared] / sst[valid_r_squared]
        )
        valid_weighted = (weighted_sst > 0) & nonconstant
        result["fit_weighted_r_squared"][columns[valid_weighted]] = (
            1 - weighted_rss[valid_weighted] / weighted_sst[valid_weighted]
        )
        _geometry_and_areas(
            result,
            columns,
            coefficients,
            x[finite],
            observed,
            fitted,
            midpoint,
            scale,
        )
    return result


def _geometry_and_areas(
    result,
    columns,
    coefficients,
    observed_x,
    observed,
    fitted,
    midpoint,
    scale,
):
    alpha, beta, gamma = coefficients
    tolerance = 64 * np.finfo(np.float64).eps * np.maximum(
        1.0,
        np.max(np.abs(observed), axis=0),
    )
    for column_index in np.flatnonzero(alpha < -tolerance):
        column = columns[column_index]
        a, b, c = alpha[column_index], beta[column_index], gamma[column_index]
        center = midpoint - scale * b / (2 * a)
        if not np.isfinite(center):
            continue
        result["index_center"][column] = center
        discriminant = b * b - 4 * a * c
        discriminant_tolerance = 64 * np.finfo(np.float64).eps * max(
            1.0,
            a * a,
            b * b,
            abs(4 * a * c),
        )
        if not np.isfinite(discriminant) or discriminant <= discriminant_tolerance:
            continue
        q = -0.5 * (b + np.copysign(np.sqrt(discriminant), b))
        roots = np.sort(midpoint + scale * np.array((q / a, c / q)))
        if not np.all(np.isfinite(roots)) or roots[0] >= roots[1]:
            continue
        left, right = roots
        result["index_left_zero"][column] = left
        result["index_right_zero"][column] = right
        support = (observed_x >= left) & (observed_x <= right)
        area_sample_count = np.count_nonzero(support)
        result["n_area_samples"][column] = area_sample_count
        if area_sample_count:
            result["Qv_fit"][column] = fitted[support, column_index].sum()
            result["Qv"][column] = observed[support, column_index].sum()


__all__ = [
    "COUNT_OUTPUTS",
    "DEFAULT_TIME_BLOCK_SIZE",
    "FLOAT_OUTPUTS",
    "analyze_velocity_profiles",
    "border_weights",
]
