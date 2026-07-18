import unittest

import numpy as np

from calculations.blood_flow_velocity import (
    BeatQualitySettings,
    PerBeatAnalysisInput,
    per_beat_signal_analysis,
    run_per_beat_analysis,
)
from input_output.schema import EyeFlowOutputPaths
from pipelines.waveform_shape_metrics.metrics.calculator import (
    WaveformShapeMetricsCalculator,
)
from pipelines.waveform_shape_metrics.metrics.runner import (
    run_waveform_shape_metric_calculations,
)
from pipelines.waveform_shape_metrics.velocity.constants import (
    LEGACY_BAND_LIMITED_SIGNAL_HARMONIC_COUNT,
)
from pipelines.waveform_shape_metrics.velocity.outputs import (
    pack_velocity_per_beat_outputs,
)


def _periodic_signal(sample_count=701, *, amplitude=1.0, phase=0.0):
    indexes = np.arange(sample_count, dtype=np.float32)
    signal = np.full(sample_count, 2.0, dtype=np.float32)
    for center in range(35, sample_count, 100):
        signal += np.float32(amplitude) * np.exp(
            np.float32(-0.5)
            * (
                (indexes - np.float32(center + phase))
                / np.float32(12.0)
            )
            ** 2
        )
    return signal.astype(np.float32)


def _quality_input(artery, vein, boundaries, *, period_ratio=None, segments=False):
    segment_kwargs = {}
    if segments:
        segment_kwargs = {
            "arterial_velocity_segments": np.stack((artery, artery * 1.1)).reshape(
                2, 1, -1
            ),
            "venous_velocity_segments": np.stack((vein, vein * 1.1)).reshape(
                2, 1, -1
            ),
        }
    return PerBeatAnalysisInput(
        arterial_velocity_signal=artery,
        venous_velocity_signal=vein,
        systolic_acceleration_peak_indexes=boundaries,
        band_limited_signal_harmonic_count=(
            LEGACY_BAND_LIMITED_SIGNAL_HARMONIC_COUNT
        ),
        dt_seconds=0.01,
        index_base=0,
        interval_period_ratio=period_ratio,
        quality_settings=BeatQualitySettings(),
        **segment_kwargs,
    )


class PerBeatQualityTests(unittest.TestCase):
    def setUp(self):
        self.artery = _periodic_signal()
        self.vein = _periodic_signal(amplitude=0.7, phase=3.0)
        self.boundaries = np.arange(0, 701, 100, dtype=np.int32)
        self.period_ratio = np.ones(self.boundaries.size - 1, dtype=np.float32)

    def test_clean_beats_are_all_accepted(self):
        result = run_per_beat_analysis(
            _quality_input(
                self.artery,
                self.vein,
                self.boundaries,
                period_ratio=self.period_ratio,
            )
        )

        np.testing.assert_array_equal(
            result.quality.accepted_mask,
            np.ones(self.boundaries.size - 1, dtype=np.bool_),
        )
        np.testing.assert_array_equal(
            result.quality.rejection_flags,
            np.zeros(self.boundaries.size - 1, dtype=np.uint16),
        )

    def test_spike_or_dip_rejects_the_paired_interval(self):
        for artifact_amplitude in (8.0, -8.0):
            with self.subTest(artifact_amplitude=artifact_amplitude):
                artery = self.artery.copy()
                artery[330:333] += np.float32(artifact_amplitude)
                result = run_per_beat_analysis(
                    _quality_input(
                        artery,
                        self.vein,
                        self.boundaries,
                        period_ratio=self.period_ratio,
                        segments=True,
                    )
                )

                expected = np.ones(self.boundaries.size - 1, dtype=np.bool_)
                expected[3] = False
                np.testing.assert_array_equal(result.quality.accepted_mask, expected)
                np.testing.assert_array_equal(
                    result.quality.accepted_mask,
                    result.quality.quality_score <= 1.0,
                )
                self.assertNotEqual(int(result.quality.rejection_flags[3]) & 4, 0)
                self.assertEqual(result.artery.signal.velocity_signal_per_beat.shape[0], 6)
                self.assertEqual(result.vein.signal.velocity_signal_per_beat.shape[0], 6)
                self.assertEqual(result.beat_period_seconds.size, 6)
                self.assertEqual(
                    result.artery.segments.velocity_signal_per_beat_per_segment.shape[1],
                    6,
                )
                self.assertEqual(
                    result.vein.segments.velocity_signal_per_beat_per_segment.shape[1],
                    6,
                )

    def test_rejected_long_interval_does_not_set_accepted_fft_length(self):
        boundaries = np.asarray([0, 100, 200, 400, 500], dtype=np.int32)
        ratios = np.asarray([1.0, 1.0, 2.0, 1.0], dtype=np.float32)
        artery = _periodic_signal(sample_count=501)
        vein = _periodic_signal(sample_count=501, amplitude=0.7, phase=3.0)

        result = run_per_beat_analysis(
            _quality_input(artery, vein, boundaries, period_ratio=ratios)
        )
        expected_mask = np.asarray([True, True, False, True])
        expected = per_beat_signal_analysis(
            artery,
            boundaries,
            LEGACY_BAND_LIMITED_SIGNAL_HARMONIC_COUNT,
            index_base=0,
            interval_mask=expected_mask,
        )

        np.testing.assert_array_equal(result.quality.accepted_mask, expected_mask)
        self.assertEqual(result.artery.signal.velocity_signal_per_beat.shape, (3, 128))
        np.testing.assert_allclose(
            result.artery.signal.velocity_signal_per_beat,
            expected.velocity_signal_per_beat,
        )

    def test_paired_broad_deformation_is_rejected_by_template_shape(self):
        indexes = np.arange(self.artery.size, dtype=np.float32)
        deformation = np.exp(
            np.float32(-0.5)
            * ((indexes - np.float32(370.0)) / np.float32(25.0)) ** 2
        ).astype(np.float32)
        artery = self.artery + deformation
        vein = self.vein + np.float32(0.7) * deformation

        result = run_per_beat_analysis(
            _quality_input(
                artery,
                vein,
                self.boundaries,
                period_ratio=self.period_ratio,
            )
        )

        expected = np.ones(self.boundaries.size - 1, dtype=np.bool_)
        expected[3] = False
        np.testing.assert_array_equal(result.quality.accepted_mask, expected)
        self.assertLess(
            float(result.quality.artery.raw_bandlimited_residual[3]),
            BeatQualitySettings().raw_bandlimited_residual_limit,
        )
        self.assertEqual(int(result.quality.rejection_flags[3]) & 48, 48)

    def test_interval_mask_length_is_validated(self):
        with self.assertRaisesRegex(ValueError, "one value per boundary interval"):
            per_beat_signal_analysis(
                self.artery,
                self.boundaries,
                LEGACY_BAND_LIMITED_SIGNAL_HARMONIC_COUNT,
                index_base=0,
                interval_mask=np.ones(2, dtype=np.bool_),
            )

    def test_all_rejected_intervals_fail_clearly(self):
        with self.assertRaisesRegex(
            ValueError,
            "No beat intervals passed waveform quality control",
        ):
            run_per_beat_analysis(
                _quality_input(
                    self.artery,
                    self.vein,
                    self.boundaries,
                    period_ratio=np.full(
                        self.boundaries.size - 1,
                        2.0,
                        dtype=np.float32,
                    ),
                )
            )

    def test_quality_outputs_are_aligned_and_auditable(self):
        artery = self.artery.copy()
        artery[330:333] -= np.float32(8.0)
        result = run_per_beat_analysis(
            _quality_input(
                artery,
                self.vein,
                self.boundaries,
                period_ratio=self.period_ratio,
            )
        )
        metrics = pack_velocity_per_beat_outputs(result)
        paths = EyeFlowOutputPaths.active().beat_quality

        self.assertIsNotNone(paths)
        np.testing.assert_array_equal(
            metrics[paths.accepted_mask][0].reshape(-1),
            result.quality.accepted_mask,
        )
        self.assertEqual(int(metrics[paths.candidate_count]), 7)
        self.assertEqual(int(metrics[paths.accepted_count]), 6)
        self.assertEqual(int(metrics[paths.rejected_count]), 1)
        _, attrs = metrics[paths.rejection_flags]
        self.assertIn("4=artery_raw_bandlimited_residual", attrs["flagDefinitions"])

        shape_metrics = run_waveform_shape_metric_calculations(metrics)
        metric_names = set(WaveformShapeMetricsCalculator._metric_names())
        scalar_metrics = {
            key: value
            for key, value in shape_metrics.items()
            if "/global/" in key and key.rsplit("/", maxsplit=1)[-1] in metric_names
        }
        self.assertEqual(len(scalar_metrics), 4 * len(metric_names))
        for key, value in scalar_metrics.items():
            self.assertEqual(np.asarray(value).size, 6, key)


if __name__ == "__main__":
    unittest.main()
