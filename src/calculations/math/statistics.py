"""Statistical array helpers shared by scientific calculations."""

from __future__ import annotations

import warnings

import numpy as np


DEFAULT_FLOAT_DTYPE = np.float32


def _as_float_array(values, dtype=DEFAULT_FLOAT_DTYPE) -> np.ndarray:
    return np.asarray(values, dtype=dtype)


def _nan_result_for_reduction(
    values: np.ndarray,
    axis: int | tuple[int, ...] | None,
    keepdims: bool,
) -> np.ndarray:
    if axis is None:
        return np.asarray(np.nan, dtype=DEFAULT_FLOAT_DTYPE)

    axes = (axis,) if isinstance(axis, int) else tuple(axis)
    normalized_axes = {
        item if item >= 0 else values.ndim + item
        for item in axes
    }
    shape = list(values.shape)
    if keepdims:
        for item in normalized_axes:
            shape[item] = 1
    else:
        shape = [
            size
            for index, size in enumerate(shape)
            if index not in normalized_axes
        ]
    return np.full(tuple(shape), np.nan, dtype=DEFAULT_FLOAT_DTYPE)


def nanmean(
    values,
    axis: int | tuple[int, ...] | None = None,
    keepdims: bool = False,
    dtype=DEFAULT_FLOAT_DTYPE,
) -> np.ndarray:
    """Return a warning-free NaN-aware mean using the shared float policy."""
    array = _as_float_array(values, dtype)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmean(array, axis=axis, keepdims=keepdims)


def nanmedian(
    values,
    axis: int | tuple[int, ...] | None = None,
    keepdims: bool = False,
    dtype=DEFAULT_FLOAT_DTYPE,
) -> np.ndarray:
    """Return a warning-free NaN-aware median using the shared float policy."""
    array = _as_float_array(values, dtype)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmedian(array, axis=axis, keepdims=keepdims)


def nansum(
    values,
    axis: int | tuple[int, ...] | None = None,
    keepdims: bool = False,
    dtype=DEFAULT_FLOAT_DTYPE,
) -> np.ndarray:
    """Return a NaN-aware sum using the shared float policy."""
    array = _as_float_array(values, dtype)
    return np.nansum(array, axis=axis, keepdims=keepdims)


def nanmin(
    values,
    axis: int | tuple[int, ...] | None = None,
    keepdims: bool = False,
    dtype=DEFAULT_FLOAT_DTYPE,
) -> np.ndarray:
    """Return a warning-free NaN-aware minimum, including empty inputs."""
    array = _as_float_array(values, dtype)
    if array.size == 0:
        return _nan_result_for_reduction(array, axis, keepdims)

    finite = np.isfinite(array)
    if axis is None and not np.any(finite):
        return np.asarray(np.nan, dtype=DEFAULT_FLOAT_DTYPE)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        result = np.nanmin(array, axis=axis, keepdims=keepdims)

    if axis is not None:
        finite_count = np.sum(finite, axis=axis, keepdims=keepdims)
        result = np.where(finite_count > 0, result, np.nan)
    return result


def nanmax(
    values,
    axis: int | tuple[int, ...] | None = None,
    keepdims: bool = False,
    dtype=DEFAULT_FLOAT_DTYPE,
) -> np.ndarray:
    """Return a warning-free NaN-aware maximum, including empty inputs."""
    array = _as_float_array(values, dtype)
    if array.size == 0:
        return _nan_result_for_reduction(array, axis, keepdims)

    finite = np.isfinite(array)
    if axis is None and not np.any(finite):
        return np.asarray(np.nan, dtype=DEFAULT_FLOAT_DTYPE)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        result = np.nanmax(array, axis=axis, keepdims=keepdims)

    if axis is not None:
        finite_count = np.sum(finite, axis=axis, keepdims=keepdims)
        result = np.where(finite_count > 0, result, np.nan)
    return result


def _nanarg_reduction(
    values: np.ndarray,
    axis: int | tuple[int, ...] | None,
    keepdims: bool,
    operation,
) -> np.ndarray:
    if values.size == 0:
        result_shape = _nan_result_for_reduction(values, axis, keepdims).shape
        return np.full(result_shape, -1, dtype=int)

    finite = np.isfinite(values)
    if isinstance(axis, tuple):
        axes = tuple(
            sorted(item if item >= 0 else values.ndim + item for item in axis)
        )
        if any(item < 0 or item >= values.ndim for item in axes):
            raise np.AxisError(axis, values.ndim)
        remaining = tuple(item for item in range(values.ndim) if item not in axes)
        permutation = remaining + axes
        reduced_values = np.transpose(values, permutation)
        reduced_finite = np.transpose(finite, permutation)
        reduced_size = int(np.prod([values.shape[item] for item in axes], dtype=int))
        leading_shape = tuple(values.shape[item] for item in remaining)
        reduced_values = reduced_values.reshape((*leading_shape, reduced_size))
        reduced_finite = reduced_finite.reshape((*leading_shape, reduced_size))
        finite_count = np.sum(reduced_finite, axis=-1)
        fill_value = np.inf if operation is np.nanargmin else -np.inf
        reduced_values = np.where(reduced_finite, reduced_values, fill_value)
        result = operation(reduced_values, axis=-1)
        result = np.where(finite_count > 0, result, -1)
        if keepdims:
            result_shape = tuple(
                1 if item in axes else values.shape[item]
                for item in range(values.ndim)
            )
            return result.reshape(result_shape)
        return result

    if axis is None:
        if not np.any(finite):
            return np.asarray(-1, dtype=int)
        fill_value = np.inf if operation is np.nanargmin else -np.inf
        reduced_values = np.where(finite, values, fill_value)
        return operation(reduced_values, axis=axis, keepdims=keepdims)

    finite_count = np.sum(finite, axis=axis, keepdims=keepdims)
    fill_value = np.inf if operation is np.nanargmin else -np.inf
    reduced_values = np.where(finite, values, fill_value)
    result = operation(reduced_values, axis=axis, keepdims=keepdims)
    return np.where(finite_count > 0, result, -1)


def nanargmin(
    values,
    axis: int | tuple[int, ...] | None = None,
    keepdims: bool = False,
    dtype=DEFAULT_FLOAT_DTYPE,
) -> np.ndarray:
    """Return a NaN-aware argmin, using ``-1`` for empty/all-NaN slices."""
    return _nanarg_reduction(
        _as_float_array(values, dtype), axis, keepdims, np.nanargmin
    )


def nanargmax(
    values,
    axis: int | tuple[int, ...] | None = None,
    keepdims: bool = False,
    dtype=DEFAULT_FLOAT_DTYPE,
) -> np.ndarray:
    """Return a NaN-aware argmax, using ``-1`` for empty/all-NaN slices."""
    return _nanarg_reduction(
        _as_float_array(values, dtype), axis, keepdims, np.nanargmax
    )


safe_nanmean = nanmean
safe_nanmedian = nanmedian


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
