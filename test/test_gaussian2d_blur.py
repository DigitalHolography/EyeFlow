"""Checks for ImageJ's radius-6 Gaussian blur of float images."""

import numpy as np
import pytest

from calculations.math.spatial_gradient import gaussian2d_blur


def test_imagej_gaussian_reference_and_frame_independence():
    stack = np.zeros((2, 129, 129), dtype=np.float32)
    stack[0, 64, 64] = 1000.0
    stack[1] = 7.0
    original = stack.copy()

    actual = gaussian2d_blur(stack)

    # GaussianBlur.makeGaussianKernel(6, 0.0002, 129), including tail correction.
    center_coefficient = 0.06649307161569595
    reference_coefficients = {
        0: center_coefficient,
        1: 0.06557594239711761,
        6: 0.04033008962869644,
        12: 0.008998858742415905,
        18: 0.0007386713405139744,
        24: 0.00001904302007460501,
        25: 0.000004760755018651253,
        26: 0.0,
    }
    for offset, coefficient in reference_coefficients.items():
        expected = 1000.0 * center_coefficient * coefficient
        np.testing.assert_allclose(actual[0, 64, 64 + offset], expected, rtol=1e-6, atol=1e-7)
        np.testing.assert_allclose(actual[0, 64 + offset, 64], expected, rtol=1e-6, atol=1e-7)
    np.testing.assert_allclose(actual[1, 1:-1, 1:-1], 7.0, atol=1e-6)
    for edge in (actual[:, 0], actual[:, -1], actual[:, :, 0], actual[:, :, -1]):
        np.testing.assert_array_equal(edge, 0.0)
    np.testing.assert_array_equal(stack, original)
    assert actual.dtype == np.float32


def test_imagej_nearest_extension_before_zeroing_output_border():
    stack = np.zeros((1, 129, 129), dtype=np.float32)
    stack[0, 0, 0] = 1.0

    actual = gaussian2d_blur(stack)

    # The replicated corner contributes the product of the two positive half-kernel sums.
    np.testing.assert_allclose(actual[0, 1, 1], 0.2178587943278677, rtol=1e-6)
    np.testing.assert_array_equal(actual[0, 0], 0.0)
    np.testing.assert_array_equal(actual[0, :, 0], 0.0)


@pytest.mark.parametrize("position", [(64, 64), (0, 0)])
def test_nan_propagates_spatially_and_remains_missing_on_border(position):
    stack = np.ones((2, 129, 129), dtype=np.float32)
    row, column = position
    stack[0, row, column] = np.nan

    actual = gaussian2d_blur(stack)

    expected_missing = np.zeros((129, 129), dtype=bool)
    expected_missing[max(0, row - 25) : row + 26, max(0, column - 25) : column + 26] = True
    np.testing.assert_array_equal(np.isnan(actual[0]), expected_missing)
    assert np.all(np.isfinite(actual[1]))


def test_empty_recording():
    stack = np.empty((0, 129, 129), dtype=np.float32)
    actual = gaussian2d_blur(stack)
    assert actual.shape == stack.shape
    assert actual.dtype == np.float32
