"""Tests for shared numerical helpers used by waveform metrics."""

import unittest
import warnings

import numpy as np

from calculations.math import (
    harmonic_pack,
    nanargmax,
    nanargmin,
    nanmax,
    nanmean,
    nanmedian,
    nanmin,
    nansum,
    rfft_normalized,
)


class SharedStatisticsTests(unittest.TestCase):
    def test_nan_reductions_are_warning_free_and_float32(self) -> None:
        values = np.asarray([[1.0, np.nan], [np.nan, np.nan]], dtype=np.float64)

        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            mean = nanmean(values, axis=0)
            median = nanmedian(values, axis=0)
            minimum = nanmin(values, axis=0)
            maximum = nanmax(values, axis=0)
            argmin = nanargmin(values, axis=0)
            argmax = nanargmax(values, axis=0)

        np.testing.assert_allclose(mean, [1.0, np.nan], equal_nan=True)
        np.testing.assert_allclose(median, [1.0, np.nan], equal_nan=True)
        np.testing.assert_allclose(minimum, [1.0, np.nan], equal_nan=True)
        np.testing.assert_allclose(maximum, [1.0, np.nan], equal_nan=True)
        np.testing.assert_array_equal(argmin, [0, -1])
        np.testing.assert_array_equal(argmax, [0, -1])
        self.assertEqual(np.float32, mean.dtype)
        self.assertEqual(np.float32, median.dtype)
        self.assertEqual(np.float32, minimum.dtype)
        self.assertEqual(np.float32, maximum.dtype)

    def test_empty_reductions_return_safe_values(self) -> None:
        empty = np.asarray([], dtype=np.float64)

        np.testing.assert_allclose(nanmin(empty), np.nan, equal_nan=True)
        np.testing.assert_allclose(nanmax(empty), np.nan, equal_nan=True)
        self.assertEqual(-1, int(nanargmin(empty)))
        self.assertEqual(-1, int(nanargmax(empty)))
        self.assertEqual(0.0, float(nansum(empty)))


class HarmonicPackTests(unittest.TestCase):
    def test_harmonic_pack_matches_normalized_fft(self) -> None:
        signal = np.asarray(
            [1.0, 2.0, 0.5, 3.0, 1.5, 0.0, 2.5, 1.0],
            dtype=np.float64,
        )
        pack = harmonic_pack(signal, max_harmonic=3, axis=0)
        expected = np.fft.rfft(signal.astype(np.float32)) / signal.size

        np.testing.assert_allclose(pack["Vfull"], expected, rtol=1e-6, atol=1e-6)
        np.testing.assert_allclose(pack["V"], expected[:4], rtol=1e-6, atol=1e-6)
        np.testing.assert_allclose(
            pack["vb"],
            np.fft.irfft(expected[:4].tolist() + [0] * (expected.size - 4), n=signal.size) * signal.size,
            rtol=1e-6,
            atol=1e-6,
        )
        self.assertEqual(3, pack["H"])

    def test_harmonic_pack_supports_batched_waveforms(self) -> None:
        signals = np.arange(16, dtype=np.float32).reshape(8, 2)
        pack = harmonic_pack(signals, max_harmonic=3, axis=0)

        self.assertEqual((4, 2), pack["V"].shape)
        self.assertEqual((5, 2), pack["Vfull"].shape)
        self.assertEqual((8, 2), pack["vb"].shape)
        np.testing.assert_allclose(
            pack["Vfull"],
            rfft_normalized(signals, axis=0),
        )


if __name__ == "__main__":
    unittest.main()
