"""Tests for generic segment profile reductions."""

from __future__ import annotations

import unittest

import numpy as np

from calculations.topology.profiles import (
    fit_inverse_parabola_profiles_with_roots,
    longitudinal_profiles,
    mean_profiles,
    profile_deviation_power,
    transverse_profiles,
)


class TopologyProfileTests(unittest.TestCase):
    def test_inverse_parabola_fit_preserves_profile_support(self) -> None:
        x = np.linspace(-90.0, 90.0, 181, dtype=np.float32)
        profile = np.full(x.shape, np.nan, dtype=np.float32)
        support = np.zeros(x.shape, dtype=bool)
        support[80:107] = True
        profile[support] = 56.0 - 0.25 * np.square(x[support] - 3.0)

        fitted, roots = fit_inverse_parabola_profiles_with_roots(
            profile,
            x_values=x,
        )

        self.assertEqual(profile.shape, fitted.shape)
        np.testing.assert_array_equal(np.isfinite(fitted), support)
        np.testing.assert_allclose(fitted[support], profile[support], atol=1e-5)
        expected_root_offset = np.sqrt(56.0 / 0.25)
        np.testing.assert_allclose(
            roots,
            [3.0 - expected_root_offset, 3.0 + expected_root_offset],
            atol=1e-5,
        )

    def test_transverse_and_longitudinal_profiles_preserve_segment_axes(self) -> None:
        segments = np.asarray(
            [[[[[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]]]],
            dtype=np.float32,
        )
        masks = np.asarray([[[[True, True, True], [False, False, False]]]])

        transverse = transverse_profiles(segments, masks)
        longitudinal = longitudinal_profiles(segments, masks)

        self.assertEqual((1, 1, 1, 3), transverse.shape)
        self.assertEqual((1, 1, 1, 2), longitudinal.shape)
        np.testing.assert_allclose(transverse, [[[[1.0, 2.0, 3.0]]]])
        np.testing.assert_allclose(
            longitudinal,
            [[[[2.0, np.nan]]]],
            equal_nan=True,
        )

    def test_profile_mean_and_deviation_power_use_the_requested_axis(self) -> None:
        profiles = np.asarray(
            [[[[1.0, 2.0, 3.0], [3.0, 4.0, 5.0]]]],
            dtype=np.float32,
        )

        mean = mean_profiles(profiles, axis=2)
        power = profile_deviation_power(profiles, axis=2)
        reused_mean_power = profile_deviation_power(profiles, mean, axis=2)

        np.testing.assert_allclose(mean, [[[2.0, 3.0, 4.0]]])
        np.testing.assert_allclose(power, np.ones_like(profiles))
        np.testing.assert_array_equal(reused_mean_power, power)
        self.assertEqual(np.float32, power.dtype)


if __name__ == "__main__":
    unittest.main()
