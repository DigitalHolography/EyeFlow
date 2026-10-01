"""Centered sliding-window reductions with periodic boundary handling."""

from __future__ import annotations

import operator
from enum import Enum

import numpy as np

try:
    from numpy.exceptions import AxisError
except ImportError:  # NumPy 1.24 compatibility
    from numpy import AxisError


class SlidingWindowMethod(str, Enum):
    """Reduction methods supported by :func:`centered_sliding_window`."""

    MEDIAN = "median"
    AVERAGE = "average"


def centered_sliding_window(
    signal,
    window_size: int,
    method: SlidingWindowMethod = SlidingWindowMethod.MEDIAN,
    *,
    window_stride: int = 1,
    axis: int = 0,
    array_module=np,
):
    """Reduce centered windows after making ``signal`` periodic along ``axis``.

    A result is produced for centers ``0, window_stride, 2 * window_stride, ...``.
    Consequently, a stride of one preserves the input shape while larger strides
    downsample the selected axis. NaNs follow the selected array module's regular
    median or mean propagation behavior.
    """

    size = _positive_integer(window_size, "window_size")
    if size % 2 == 0:
        raise ValueError("window_size must be odd for a centered sliding window.")
    stride = _positive_integer(window_stride, "window_stride")
    try:
        reduction = SlidingWindowMethod(method)
    except ValueError as error:
        choices = ", ".join(item.value for item in SlidingWindowMethod)
        raise ValueError(
            f"Unsupported sliding-window method {method!r}; expected one of: {choices}."
        ) from error

    xp = array_module
    values = xp.asarray(signal)
    normalized_axis = _normalize_axis(axis, values.ndim)
    axis_length = int(values.shape[normalized_axis])
    if axis_length == 0:
        return values.copy()

    radius = size // 2
    centers = xp.arange(0, axis_length, stride)
    offsets = xp.arange(-radius, radius + 1)
    periodic_indices = (centers[:, None] + offsets[None, :]) % axis_length
    windows = xp.take(values, periodic_indices, axis=normalized_axis)
    window_axis = normalized_axis + 1

    if reduction is SlidingWindowMethod.MEDIAN:
        return xp.median(windows, axis=window_axis)
    return xp.mean(windows, axis=window_axis)


def _positive_integer(value: int, name: str) -> int:
    try:
        normalized = operator.index(value)
    except TypeError as error:
        raise TypeError(f"{name} must be an integer.") from error
    if normalized < 1:
        raise ValueError(f"{name} must be positive.")
    return normalized


def _normalize_axis(axis: int, ndim: int) -> int:
    try:
        normalized = operator.index(axis)
    except TypeError as error:
        raise TypeError("axis must be an integer.") from error
    if normalized < -ndim or normalized >= ndim:
        raise AxisError(normalized, ndim=ndim)
    return normalized % ndim


__all__ = ["SlidingWindowMethod", "centered_sliding_window"]
