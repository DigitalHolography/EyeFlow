"""Tests for MATLAB-grounded cardiac-cycle calculations."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from calculations.blood_flow_velocity.signal_analysis.cardiac_cycle import (  # noqa: E402
    MATLAB_PADDING_FACTOR,
    SystoleDetectionError,
    analyze_cardiac_cycles,
    cardiac_cycles_from_available_vessel,
    find_systole_index,
    spectral_cardiac_cycle_analysis,
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
            expected_filtered = fixture["filtered_velocity_mm_s"]
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
            velocity, detection.systole_indexes, 13, index_base=0,
        )
        self.assertEqual((4, 128), per_beat.velocity_signal_per_beat.shape)
        for index, (start, stop) in enumerate(zip(expected[:-1], expected[1:])):
            np.testing.assert_array_equal(
                per_beat.velocity_signal_per_beat[index],
                interpft_real(velocity[start : stop + 1], 129)[:-1],
            )
        peaks = detection.systole_indexes
        self.assertTrue(np.all(detection.derivative_signal[peaks] >= detection.min_peak_height))

    def test_edge_case_boundaries_are_invariant_to_physical_velocity_scale(self) -> None:
        with np.load(Path(__file__).parent / "fixtures" / "cardiac_velocity_edge_case.npz") as fixture:
            velocity = fixture["arterial_velocity_mm_s"]
            dt = np.float32(fixture["dt_seconds"])
            expected = fixture["expected_boundaries"]
        for scale in (1e-4, 1.0, 1e4):
            with self.subTest(scale=scale):
                detection = find_systole_index(velocity * np.float32(scale), dt=dt)
                np.testing.assert_array_equal(detection.systole_indexes, expected)

    def test_strict_height_gate_reports_a_missed_lower_amplitude_beat(self) -> None:
        # A strict global percentile can reject a genuine weaker
        # upstroke. Preserve and warn about the resulting two-period interval.
        frames = np.arange(700, dtype=np.float32)
        waveform = np.full(frames.size, -0.5, dtype=np.float32)
        centers = np.arange(50, 651, 100)
        for center, amplitude in zip(centers, (1.0, .8, 1.2, .9, 1.1, .85, 1.05)):
            waveform += amplitude * np.exp(-0.5 * ((frames - center) / 7.0) ** 2)
        detection = find_systole_index(waveform, dt=np.float32(.01))
        np.testing.assert_array_equal(detection.systole_indexes, [43, 243, 343, 443, 543, 643])
        self.assertEqual(len(detection.suspected_missed_beat_gaps), 1)
        gap = detection.suspected_missed_beat_gaps[0]
        self.assertEqual((gap.start_index, gap.stop_index, gap.estimated_multiple), (43, 243, 2))

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

    def test_matlab_padding_and_harmonic_fit(self) -> None:
        dt_seconds = 0.05
        time = np.arange(200, dtype=np.float32) * dt_seconds
        waveform = (
            1.0
            + np.sin(2.0 * np.pi * 1.0 * time)
            + 0.45 * np.sin(2.0 * np.pi * 2.0 * time)
            + 0.25 * np.sin(2.0 * np.pi * 3.0 * time)
        ).astype(np.float32)

        result = spectral_cardiac_cycle_analysis(
            waveform,
            dt_seconds,
            systole_count=10,
        )

        expected_fft_size = waveform.size * (1 + 2 * MATLAB_PADDING_FACTOR)
        self.assertEqual(result.fft_coefficients.size, expected_fft_size)
        self.assertEqual(result.frequencies_hz.size, expected_fft_size // 2 + 1)
        self.assertAlmostEqual(result.estimated_fundamental_hz, 1.0)
        self.assertAlmostEqual(result.fundamental_hz, 1.0, places=6)
        self.assertAlmostEqual(result.heart_rate_hz, 1.0, places=6)
        self.assertAlmostEqual(result.heart_rate_bpm, 60.0, places=5)
        self.assertAlmostEqual(result.period_seconds, 1.0, places=6)
        np.testing.assert_allclose(
            result.frequencies_hz[result.peak_indexes],
            [1.0, 2.0, 3.0],
            atol=1e-6,
        )

    def test_lowest_prominent_peak_is_the_matlab_fundamental(self) -> None:
        dt_seconds = 0.025
        time = np.arange(400, dtype=np.float32) * dt_seconds
        waveform = (
            0.3 * np.sin(2.0 * np.pi * 1.2 * time)
            + np.sin(2.0 * np.pi * 2.4 * time)
        ).astype(np.float32)

        result = spectral_cardiac_cycle_analysis(
            waveform,
            dt_seconds,
            systole_count=12,
        )

        # Reference values generated by MATLAB R2026a from the current
        # BloodFlowVelocity/SpectralWaveformAnalysis.m (N=2, prominence=10%).
        self.assertAlmostEqual(
            result.fundamental_hz,
            1.2,
            delta=1e-12,
        )
        self.assertAlmostEqual(
            result.heart_rate_hz,
            1.2,
            delta=1e-12,
        )
        self.assertAlmostEqual(
            result.heart_rate_ste_hz,
            0.0,
            delta=1e-12,
        )
        np.testing.assert_allclose(
            result.frequencies_hz[result.peak_indexes],
            [1.2, 2.4],
            atol=1e-7,
        )
        self.assertAlmostEqual(result.period_seconds, 1.0 / result.heart_rate_hz)

    def test_cardiac_cycle_analysis_combines_boundaries_and_spectral_rate(self) -> None:
        dt_seconds = 0.01
        time = np.arange(1000, dtype=np.float32) * dt_seconds
        waveform = (10.0 + 2.0 * np.sin(2.0 * np.pi * 1.2 * time)).astype(
            np.float32
        )

        result = analyze_cardiac_cycles(
            waveform,
            dt_seconds=dt_seconds,
        )

        self.assertGreater(result.systole.systole_indexes.size, 2)
        self.assertAlmostEqual(result.spectral.heart_rate_hz, 1.2, places=6)

    def test_exercise_rate_uses_median_adaptive_spacing_and_64_point_beats(
        self,
    ) -> None:
        dt_seconds = 0.01
        time = np.arange(1000, dtype=np.float32) * dt_seconds
        waveform = (
            10.0
            + 2.0 * np.sin(2.0 * np.pi * 2.5 * time)
            + 0.25 * np.sin(2.0 * np.pi * 5.0 * time)
        ).astype(np.float32)

        detection = find_systole_index(
            waveform,
            dt=np.float32(dt_seconds),
        )
        per_beat = per_beat_signal_analysis(
            waveform,
            detection.systole_indexes,
            band_limited_signal_harmonic_count=4,
        )

        np.testing.assert_allclose(
            np.diff(detection.systole_indexes),
            40,
            atol=1,
        )
        self.assertAlmostEqual(detection.estimated_period_samples, 40.0, delta=1.0)
        self.assertGreaterEqual(detection.min_peak_distance, 35)
        self.assertLess(detection.min_peak_distance, 40)
        self.assertEqual(per_beat.velocity_signal_per_beat.shape[1], 64)

    def test_stable_cadence_rejects_a_close_weaker_false_upstroke(self) -> None:
        dt_seconds = 0.006912
        frame_indexes = np.arange(512, dtype=np.float32)
        waveform = np.full(frame_indexes.size, 20.0, dtype=np.float32)
        for center, amplitude in (
            (24, 20.0),
            (103, 30.0),
            (225, 30.0),
            (348, 30.0),
            (469, 30.0),
        ):
            waveform += amplitude * np.exp(
                -0.5 * ((frame_indexes - center) / 12.0) ** 2
            )

        detection = find_systole_index(
            waveform,
            dt=np.float32(dt_seconds),
        )

        np.testing.assert_allclose(
            detection.systole_indexes,
            [91, 213, 336, 457],
            atol=1,
        )
        np.testing.assert_allclose(
            np.diff(detection.systole_indexes),
            [122, 123, 121],
            atol=1,
        )

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

    def test_two_or_three_period_gaps_are_retained_and_reported(self) -> None:
        dt_seconds = 0.01
        time = np.arange(1000, dtype=np.float32) * dt_seconds

        for omitted_beats, expected_multiple in (((8,), 2), ((8, 9), 3)):
            with self.subTest(expected_multiple=expected_multiple):
                waveform = np.zeros(time.size, dtype=np.float32)
                for beat_index, center in enumerate(np.arange(0.3, 9.8, 0.4)):
                    if beat_index in omitted_beats:
                        continue
                    waveform += np.exp(
                        -0.5 * ((time - center) / 0.04) ** 2
                    ).astype(np.float32)

                detection = find_systole_index(
                    waveform,
                    dt=np.float32(dt_seconds),
                )

                self.assertEqual(len(detection.suspected_missed_beat_gaps), 1)
                gap = detection.suspected_missed_beat_gaps[0]
                self.assertEqual(gap.estimated_multiple, expected_multiple)
                self.assertEqual(gap.interval_samples, expected_multiple * 40)

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
