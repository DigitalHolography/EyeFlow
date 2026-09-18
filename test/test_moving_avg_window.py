"""Tests for truncated, NaN-propagating temporal moving averages."""

import numpy as np
import pytest

from calculations.math.spatial_gradient import moving_avg_window


def test_seven_frame_average_truncates_edges_and_filters_pixels_independently():
    stack = np.zeros((9, 1, 2), dtype=np.float32)
    stack[3, 0, 0] = 70.0
    stack[:, 0, 1] = 3.0

    actual = moving_avg_window(stack)

    # The impulse is divided by the actual window sizes: 4, 5, 6, 7, ...
    np.testing.assert_allclose(
        actual[:, 0, 0], [17.5, 14.0, 70.0 / 6.0, 10.0, 10.0, 10.0, 70.0 / 6.0, 0.0, 0.0]
    )
    np.testing.assert_array_equal(actual[:, 0, 1], 3.0)
    assert actual.dtype == np.float32
    assert stack[3, 0, 0] == 70.0


def test_nan_propagates_only_to_windows_containing_it():
    stack = np.broadcast_to(np.arange(11, dtype=np.float32)[:, None, None], (11, 1, 2)).copy()
    stack[4, 0, 0] = np.nan

    actual = moving_avg_window(stack)

    np.testing.assert_array_equal(
        np.isnan(actual[:, 0, 0]),
        [False, True, True, True, True, True, True, True, False, False, False],
    )
    np.testing.assert_allclose(actual[[0, 8, 9, 10], 0, 0], [1.5, 7.5, 8.0, 8.5])
    assert np.all(np.isfinite(actual[:, 0, 1]))


@pytest.mark.parametrize("values", [[8.0], [0.0, 1.0, 4.0, 8.0]])
def test_recordings_shorter_than_window_average_only_available_frames(values):
    stack = np.asarray(values, dtype=np.float32)[:, None, None]

    actual = moving_avg_window(stack)

    np.testing.assert_allclose(actual, np.mean(values))
    assert actual.shape == stack.shape


def test_empty_recording_and_single_frame_window():
    empty = np.empty((0, 2, 3), dtype=np.float32)
    assert moving_avg_window(empty).shape == empty.shape
    stack = np.asarray([1.0, np.nan, 3.0], dtype=np.float32)[:, None, None]
    np.testing.assert_array_equal(moving_avg_window(stack, window=1), stack)
