"""Tests for MATLAB-grounded cardiac-cycle calculations."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np
from scipy.signal import butter, filtfilt

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from calculations.blood_flow_velocity.signal_analysis.cardiac_cycle import (  # noqa: E402
    SystoleDetectionError,
    analyze_cardiac_cycles,
    cardiac_cycles_from_available_vessel,
    find_systole_index,
)
from calculations.blood_flow_velocity.signal_analysis.per_beat.signal import (  # noqa: E402
    per_beat_signal_analysis,
)
from calculations.math import interpft_real  # noqa: E402


class SpectralCardiacCycleTests(unittest.TestCase):
    def test_edge_case_velocity_yields_four_complete_correctly_aligned_beats(self) -> None:
        with np.load(Path(__file__).parent / "fixtures" / "cardiac_velocity_edge_case.npz") as fixture:
            velocity = fixture["arterial_velocity_mm_s"]
            dt_seconds = float(fixture["dt_seconds"])
            expected = fixture["expected_boundaries"]
        # The fixture's filtered array used SciPy's default 15-sample padding.
        # Keep its boundary reference, but verify MATLAB's 12-sample reflection.
        cutoff = np.float32(15) / (np.float32(.5) / np.float32(dt_seconds))
        b, a = butter(4, cutoff)
        expected_filtered = filtfilt(b, a, velocity, padlen=12).astype(np.float32)
        original = velocity.copy()
        analysis = analyze_cardiac_cycles(velocity, dt_seconds=dt_seconds)
        detection = analysis.systole
        np.testing.assert_array_equal(detection.systole_indexes, expected)
        np.testing.assert_array_equal(detection.signal_filtered, expected_filtered)
        np.testing.assert_array_equal(velocity, original)
        self.assertAlmostEqual(
            detection.min_peak_height, np.percentile(detection.derivative_signal, 95), places=5,
        )
        np.testing.assert_allclose(
            np.diff(detection.systole_indexes) * dt_seconds,
            [0.7188487, 0.7188487, 0.7672328, 0.7810568], atol=1e-6,
        )
        per_beat = per_beat_signal_analysis(
            velocity, detection.systole_indexes, 12, index_base=0,
        )
        self.assertEqual((4, 128), per_beat.velocity_signal_per_beat.shape)
        for index, (start, stop) in enumerate(zip(expected[:-1], expected[1:])):
            np.testing.assert_array_equal(
                per_beat.velocity_signal_per_beat[index],
                interpft_real(velocity[start : stop + 1], 129)[:-1],
            )
        peaks = detection.systole_indexes
        self.assertTrue(np.all(detection.derivative_signal[peaks] >= detection.min_peak_height))


    def test_unusable_artery_uses_vein_then_full_record(self) -> None:
        time = np.arange(1000, dtype=np.float32) * .01
        vein = np.sin(2 * np.pi * time).astype(np.float32)
        unavailable = np.full(vein.size, np.nan, dtype=np.float32)
        analysis, source = cardiac_cycles_from_available_vessel(
            unavailable, vein, dt_seconds=.01,
        )
        self.assertEqual(source, "vein")
        np.testing.assert_array_equal(
            analysis.systole.systole_indexes,
            find_systole_index(vein, dt=np.float32(.01)).systole_indexes,
        )
        analysis, source = cardiac_cycles_from_available_vessel(
            unavailable, unavailable, dt_seconds=.01,
        )
        self.assertEqual(source, "none")
        np.testing.assert_array_equal(analysis.systole.systole_indexes, [0, vein.size - 1])

    def test_systole_derivative_is_velocity_per_second(self) -> None:
        dt_seconds = 0.01
        time = np.arange(1000, dtype=np.float32) * dt_seconds
        waveform = (10.0 + time + np.sin(2.0 * np.pi * time)).astype(np.float32)

        result = analyze_cardiac_cycles(waveform, dt_seconds=dt_seconds)

        expected = np.gradient(
            result.systole.signal_filtered,
            np.float32(dt_seconds),
        )
        np.testing.assert_allclose(result.systole.derivative_signal, expected)


    def test_variable_cadence_is_not_forced_to_the_median_period(self) -> None:
        dt_seconds = 0.01
        expected_intervals = np.asarray(
            [32, 40, 48, 40, 32, 48, 40, 40],
            dtype=np.int32,
        )
        centers = np.concatenate(
            (np.asarray([30], dtype=np.int32), 30 + np.cumsum(expected_intervals))
        )
        frame_indexes = np.arange(int(centers[-1]) + 40, dtype=np.float32)
        waveform = np.full(frame_indexes.size, 20.0, dtype=np.float32)
        for center in centers:
            waveform += 30.0 * np.exp(
                -0.5 * ((frame_indexes - center) / 5.0) ** 2
            )

        detection = find_systole_index(
            waveform,
            dt=np.float32(dt_seconds),
        )

        np.testing.assert_array_equal(
            np.diff(detection.systole_indexes),
            expected_intervals,
        )


    def test_one_detected_peak_uses_full_record_fallback(self) -> None:
        dt_seconds = 0.01
        time = np.arange(150, dtype=np.float32) * dt_seconds
        waveform = np.sin(2.0 * np.pi * time).astype(np.float32)

        with self.assertRaises(SystoleDetectionError):
            find_systole_index(waveform, dt=np.float32(dt_seconds))

        analysis, source = cardiac_cycles_from_available_vessel(
            waveform,
            waveform,
            dt_seconds=dt_seconds,
        )

        self.assertEqual(source, "none")
        np.testing.assert_array_equal(
            analysis.systole.systole_indexes,
            [0, waveform.size - 1],
        )


if __name__ == "__main__":
    unittest.main()
