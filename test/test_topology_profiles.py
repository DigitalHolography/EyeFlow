"""Tests for generic segment profile reductions."""

from __future__ import annotations

import unittest

import numpy as np

from calculations.topology.profiles import (
    longitudinal_profiles,
    mean_profiles,
    profile_deviation_power,
    transverse_profiles,
)


class TopologyProfileTests(unittest.TestCase):
    def test_spatial_mask_broadcasts_over_a_frame_stack(self) -> None:
        segments = np.arange(2 * 3 * 4, dtype=np.float32).reshape(2, 3, 4)
        mask = np.zeros((3, 4), dtype=bool)
        mask[1, 1:3] = True

        transverse = transverse_profiles(segments, mask)
        longitudinal = longitudinal_profiles(segments, mask)

        expected_transverse = np.full((2, 4), np.nan, dtype=np.float32)
        expected_transverse[:, 1:3] = segments[:, 1, 1:3]
        expected_longitudinal = np.full((2, 3), np.nan, dtype=np.float32)
        expected_longitudinal[:, 1] = np.mean(segments[:, 1, 1:3], axis=1)
        np.testing.assert_allclose(transverse, expected_transverse, equal_nan=True)
        np.testing.assert_allclose(longitudinal, expected_longitudinal, equal_nan=True)

    def test_cuda_profiles_stay_on_device_and_match_cpu(self) -> None:
        from calculations.compute_backend import optional_cupy_backend
        backend = optional_cupy_backend()
        if backend is None:
            self.skipTest("CuPy/CUDA unavailable.")
        cupy = backend.cupy
        segments = np.arange(2 * 3 * 4, dtype=np.float32).reshape(1, 1, 2, 3, 4)
        segments[:, :, 0] = np.nan
        masks = np.zeros((1, 1, 3, 4), bool)
        masks[..., 1, 1:3] = True
        for reduce in (transverse_profiles, longitudinal_profiles):
            expected = reduce(segments, masks)
            actual = reduce(cupy.asarray(segments), masks)
            self.assertIsInstance(actual, cupy.ndarray)
            self.assertEqual(cupy.float32, actual.dtype)
            np.testing.assert_allclose(cupy.asnumpy(actual), expected, equal_nan=True)

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
