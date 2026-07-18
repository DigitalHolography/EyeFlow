import unittest

import numpy as np

from calculations.blood_flow_velocity.analysis_preparation.beat_detection import (
    find_systole_index,
)
from calculations.dopplerview_analysis import analyze_arterial_waveforms


def _pulse_train(amplitudes, *, centers=None, sample_count=800):
    x = np.arange(sample_count, dtype=np.float32)
    if centers is None:
        centers = np.arange(100, 100 * (len(amplitudes) + 1), 100)
    signal = np.zeros(sample_count, dtype=np.float32)
    for amplitude, center in zip(amplitudes, centers, strict=True):
        signal += np.float32(amplitude) * np.exp(
            np.float32(-0.5) * ((x - np.float32(center)) / np.float32(12.0)) ** 2
        )
    return signal.astype(np.float32)


class RobustBeatDetectionTests(unittest.TestCase):
    def test_clean_primary_boundaries_are_unchanged(self):
        result = find_systole_index(
            _pulse_train(np.ones(7)),
            dt=np.float32(0.01),
        )

        expected = np.arange(88, 689, 100, dtype=np.int32)
        np.testing.assert_array_equal(result.initial_systole_indexes, expected)
        np.testing.assert_array_equal(result.systole_indexes, expected)
        self.assertEqual(result.recovered_systole_indexes.size, 0)
        self.assertEqual(float(result.nominal_period_samples), 100.0)
        np.testing.assert_array_equal(
            result.interval_duration_valid,
            np.ones(6, dtype=np.bool_),
        )

    def test_attenuated_midpoint_upstroke_is_recovered_and_flagged(self):
        amplitudes = np.ones(7, dtype=np.float32)
        amplitudes[3] = np.float32(0.8)
        result = find_systole_index(
            _pulse_train(amplitudes),
            dt=np.float32(0.01),
        )

        self.assertNotIn(388, result.initial_systole_indexes)
        np.testing.assert_array_equal(
            result.recovered_systole_indexes,
            np.asarray([388], dtype=np.int32),
        )
        np.testing.assert_array_equal(
            result.systole_indexes,
            np.arange(88, 689, 100, dtype=np.int32),
        )
        self.assertTrue(np.all(result.interval_duration_valid))

    def test_true_pause_without_supporting_upstroke_is_not_filled(self):
        centers = np.asarray([100, 200, 300, 500, 600, 700], dtype=np.int32)
        result = find_systole_index(
            _pulse_train(np.ones(centers.size), centers=centers),
            dt=np.float32(0.01),
        )

        self.assertEqual(result.recovered_systole_indexes.size, 0)
        self.assertIn(200, np.diff(result.systole_indexes))
        self.assertFalse(np.all(result.interval_duration_valid))

    def test_weak_midpoint_fluctuation_is_not_promoted_to_a_boundary(self):
        amplitudes = np.ones(7, dtype=np.float32)
        amplitudes[3] = np.float32(0.55)
        result = find_systole_index(
            _pulse_train(amplitudes),
            dt=np.float32(0.01),
        )

        self.assertEqual(result.recovered_systole_indexes.size, 0)
        self.assertIn(200, np.diff(result.systole_indexes))
        self.assertFalse(np.all(result.interval_duration_valid))

    def test_moderately_irregular_rhythm_is_not_modified(self):
        centers = np.cumsum(
            np.asarray([100, 90, 112, 95, 125, 88, 108], dtype=np.int32)
        )
        result = find_systole_index(
            _pulse_train(np.ones(centers.size), centers=centers),
            dt=np.float32(0.01),
        )

        self.assertEqual(result.recovered_systole_indexes.size, 0)
        np.testing.assert_array_equal(
            result.systole_indexes,
            result.initial_systole_indexes,
        )

    def test_affine_signal_transform_preserves_boundaries(self):
        amplitudes = np.ones(7, dtype=np.float32)
        amplitudes[3] = np.float32(0.8)
        signal = _pulse_train(amplitudes)
        baseline = find_systole_index(signal, dt=np.float32(0.01))
        transformed = find_systole_index(
            signal * np.float32(7.0) + np.float32(23.0),
            dt=np.float32(0.01),
        )

        np.testing.assert_array_equal(
            transformed.systole_indexes,
            baseline.systole_indexes,
        )
        np.testing.assert_array_equal(
            transformed.recovered_systole_indexes,
            baseline.recovered_systole_indexes,
        )

    def test_waveform_refresh_exposes_aligned_diagnostics(self):
        amplitudes = np.ones(7, dtype=np.float32)
        amplitudes[3] = np.float32(0.8)
        artery = _pulse_train(amplitudes)
        outputs = analyze_arterial_waveforms(
            artery,
            artery * np.float32(0.7),
            dt_seconds=np.float32(0.01),
        )

        boundaries = np.asarray(outputs["beat_indices"])
        ratios = np.asarray(outputs["beat_detection_interval_period_ratio"])
        valid = np.asarray(outputs["beat_detection_interval_duration_valid"])
        self.assertEqual(ratios.size, boundaries.size - 1)
        self.assertEqual(valid.size, boundaries.size - 1)
        np.testing.assert_array_equal(
            outputs["beat_detection_recovered_indices"],
            np.asarray([388], dtype=np.int32),
        )


if __name__ == "__main__":
    unittest.main()
