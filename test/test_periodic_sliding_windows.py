"""Tests for centered periodic sliding-window reductions."""

import numpy as np
import pytest

from calculations.math.periodic import (
    SlidingWindowMethod,
    centered_sliding_window,
)


def test_average_wraps_centered_windows_around_both_edges():
    stack = np.zeros((9, 1, 2), dtype=np.float32)
    stack[3, 0, 0] = 70.0
    stack[:, 0, 1] = 3.0

    actual = centered_sliding_window(
        stack,
        7,
        SlidingWindowMethod.AVERAGE,
    )

    np.testing.assert_allclose(
        actual[:, 0, 0],
        [10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 0.0, 0.0],
    )
    np.testing.assert_array_equal(actual[:, 0, 1], 3.0)
    assert actual.dtype == np.float32
    assert stack[3, 0, 0] == 70.0


def test_median_filters_each_signal_column_independently():
    signal = np.asarray(
        [[8.0, 1.0], [1.0, 10.0], [7.0, 3.0], [2.0, 20.0]],
        dtype=np.float32,
    )

    actual = centered_sliding_window(
        signal,
        3,
        SlidingWindowMethod.MEDIAN,
    )

    np.testing.assert_array_equal(
        actual,
        [[2.0, 10.0], [7.0, 3.0], [2.0, 10.0], [7.0, 3.0]],
    )


def test_nan_propagates_only_to_periodic_windows_containing_it():
    stack = np.broadcast_to(
        np.arange(11, dtype=np.float32)[:, None, None],
        (11, 1, 2),
    ).copy()
    stack[4, 0, 0] = np.nan

    actual = centered_sliding_window(
        stack,
        7,
        SlidingWindowMethod.AVERAGE,
    )

    np.testing.assert_array_equal(
        np.isnan(actual[:, 0, 0]),
        [False, True, True, True, True, True, True, True, False, False, False],
    )
    assert np.all(np.isfinite(actual[:, 0, 1]))


def test_axis_and_stride_downsample_window_centers():
    signal = np.arange(12, dtype=np.float32).reshape(2, 6)

    actual = centered_sliding_window(
        signal,
        3,
        SlidingWindowMethod.AVERAGE,
        window_stride=2,
        axis=1,
    )

    np.testing.assert_allclose(actual, [[2.0, 2.0, 4.0], [8.0, 8.0, 10.0]])
    assert actual.shape == (2, 3)


def test_empty_signal_and_single_sample_window_preserve_shape():
    empty = np.empty((0, 2, 3), dtype=np.float32)
    assert centered_sliding_window(empty, 9).shape == empty.shape
    signal = np.asarray([1.0, np.nan, 3.0], dtype=np.float32)[:, None, None]
    np.testing.assert_array_equal(
        centered_sliding_window(signal, 1, SlidingWindowMethod.AVERAGE),
        signal,
    )


@pytest.mark.parametrize("window_size", [0, -1, 2, 4])
def test_window_size_must_be_positive_and_odd(window_size):
    with pytest.raises(ValueError):
        centered_sliding_window(np.arange(3), window_size)


@pytest.mark.parametrize("window_stride", [0, -1])
def test_window_stride_must_be_positive(window_stride):
    with pytest.raises(ValueError, match="window_stride"):
        centered_sliding_window(
            np.arange(3),
            3,
            window_stride=window_stride,
        )


def test_method_and_axis_are_validated_for_empty_signals():
    with pytest.raises(ValueError, match="Unsupported"):
        centered_sliding_window(np.empty((0, 2)), 3, "mode")
    with pytest.raises(IndexError):
        centered_sliding_window(np.empty((0, 2)), 3, axis=2)
