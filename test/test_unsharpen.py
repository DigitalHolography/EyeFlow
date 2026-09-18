"""Checks for ImageJ unsharp-mask strength, spatial support and missing data."""

import numpy as np
import pytest

from calculations.math.spatial_gradient import unsharpen


def test_imagej_mask_weight_and_gaussian_tail_reference_without_temporal_filtering():
    stack = np.full((2, 129, 129), 100.0, dtype=np.float32)
    stack[0, 64, 64] += 1000.0
    stack[1] = 7.0
    original = stack.copy()

    actual = unsharpen(stack)

    # GaussianBlur.makeGaussianKernel(8, 0.01, 129), including its tail correction.
    # ImageJ's mask weight of 0.6 corresponds to an ordinary sharpening amount of 1.5.
    center_coefficient = 0.050052348524332047
    reference_coefficients = {
        0: center_coefficient,
        1: 0.049662839621305466,
        8: 0.03035828471183777,
        16: 0.00677384901791811,
        20: 0.002191213658079505,
        24: 0.00024346818099729717,
        25: 0.00006086704524932429,
        26: 0.0,
    }
    for offset, coefficient in reference_coefficients.items():
        blurred = 100.0 + 1000.0 * center_coefficient * coefficient
        source = 1100.0 if offset == 0 else 100.0
        expected = (source - 0.6 * blurred) / 0.4
        np.testing.assert_allclose(actual[0, 64, 64 + offset], expected, rtol=1e-6, atol=1e-5)
    np.testing.assert_allclose(actual[1, 1:-1, 1:-1], 7.0, atol=1e-6)
    for edge in (actual[:, 0], actual[:, -1], actual[:, :, 0], actual[:, :, -1]):
        np.testing.assert_array_equal(edge, 0.0)
    np.testing.assert_array_equal(stack, original)
    assert actual.dtype == np.float32


@pytest.mark.parametrize("position", [(64, 64), (0, 0)])
def test_nan_propagates_through_spatial_kernel_and_stays_missing_at_border(position):
    stack = np.ones((2, 129, 129), dtype=np.float32)
    row, column = position
    stack[0, row, column] = np.nan

    actual = unsharpen(stack)

    expected_missing = np.zeros((129, 129), dtype=bool)
    expected_missing[max(0, row - 25) : row + 26, max(0, column - 25) : column + 26] = True
    np.testing.assert_array_equal(np.isnan(actual[0]), expected_missing)
    assert np.all(np.isfinite(actual[1]))
    assert np.isnan(actual[0, row, column])


def test_negative_halos_are_clipped_and_zero_signal_stays_zero():
    stack = np.zeros((1, 129, 129), dtype=np.float32)
    stack[0, 64, 64] = 1.0

    actual = unsharpen(stack)

    assert 2.0 < actual[0, 64, 64] < 2.5
    assert np.count_nonzero(actual) == 1
    assert np.all(actual >= 0.0)


def test_empty_recording():
    stack = np.empty((0, 129, 129), dtype=np.float32)
    actual = unsharpen(stack)
    assert actual.shape == stack.shape
    assert actual.dtype == np.float32
