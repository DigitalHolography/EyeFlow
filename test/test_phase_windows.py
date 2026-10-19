"""Configurable closed phase windows, NaN-aware medians, and signed deltas."""

import warnings

import numpy as np
import pytest

from calculations.math.phase_windows import (
    early_late_phase_windows, phase_window_medians, late_minus_early,
)


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_default_windows_include_the_shared_third_boundary(dtype):
    windows = early_late_phase_windows(np.arange(6, dtype=dtype) / 6)
    np.testing.assert_array_equal(windows["early"][0], [True, True, True, False, False, False])
    np.testing.assert_array_equal(windows["late"][0], [False, False, True, True, True, True])
    assert windows["early"][1] == [0, 1 / 3]
    assert windows["late"][1] == [1 / 3, 1]


def test_custom_window_bounds_are_closed_and_preserved_as_fractions():
    windows = early_late_phase_windows(
        np.array([0, 0.25, 0.5, 0.75]), early_interval=(0.25, 0.5), late_interval=(0.5, 1),
    )
    np.testing.assert_array_equal(windows["early"][0], [False, True, True, False])
    np.testing.assert_array_equal(windows["late"][0], [False, False, True, True])
    assert windows["early"][1] == [0.25, 0.5]
    assert windows["late"][1] == [0.5, 1]


def test_generic_medians_support_arbitrary_phase_axis_empty_windows_and_missing_cells():
    values = np.array([[1, np.nan, 5, 9], [np.nan, np.nan, 2, 4]])
    phase = np.arange(4) / 4
    windows = early_late_phase_windows(phase, early_interval=(0, 0.25), late_interval=(0.5, 1))
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        medians = phase_window_medians(values, windows, axis=1)
        np.testing.assert_array_equal(medians["early"], [1, np.nan])
        np.testing.assert_array_equal(medians["late"], [7, 3])
        np.testing.assert_array_equal(late_minus_early(medians["early"], medians["late"]), [6, np.nan])
        empty = early_late_phase_windows(phase, early_interval=(0.1, 0.2))
        assert np.isnan(phase_window_medians(values, empty, axis=1)["early"]).all()
    assert late_minus_early(5, 2) == -3


@pytest.mark.parametrize("interval", [(-0.1, 0.3), (0.7, 0.2), (0.2, 1.1), (np.nan, 1), (0.2,), (0.2, 0.3, 0.4)])
def test_invalid_fractional_intervals_are_rejected(interval):
    with pytest.raises(ValueError, match="two ordered finite bounds"):
        early_late_phase_windows([0, 0.25, 0.5], early_interval=interval)


@pytest.mark.parametrize("phase", [[], [0, 1], [0.5, 0.25], [0, 0], [0, np.nan]])
def test_invalid_source_phase_coordinates_are_rejected(phase):
    with pytest.raises(ValueError, match="strictly increasing"):
        early_late_phase_windows(phase)


def test_phase_mask_shape_must_match_the_reduction_axis():
    windows = early_late_phase_windows(np.arange(4) / 4)
    with pytest.raises(ValueError, match="match the phase axis"):
        phase_window_medians(np.ones((3, 2)), windows)
