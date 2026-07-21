import unittest

import numpy as np

from calculations.blood_flow_velocity import PerBeatAnalysisInput, run_per_beat_analysis
from calculations.blood_flow_velocity.analysis_preparation.beat_detection import (
    find_systole_index,
)
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


def _velocity_signal(amplitudes):
    indexes = np.arange(800, dtype=np.float32)
    signal = np.full(800, 2.0, dtype=np.float32)
    for amplitude, center in zip(
        amplitudes,
        np.arange(100, 701, 100),
        strict=True,
    ):
        signal += np.float32(amplitude) * np.exp(
            np.float32(-0.5)
            * ((indexes - np.float32(center)) / np.float32(12.0)) ** 2
        )
    return signal.astype(np.float32)


def _metric_outputs(artery, vein, boundaries):
    result = run_per_beat_analysis(
        PerBeatAnalysisInput(
            arterial_velocity_signal=artery,
            venous_velocity_signal=vein,
            systolic_acceleration_peak_indexes=boundaries,
            band_limited_signal_harmonic_count=(
                LEGACY_BAND_LIMITED_SIGNAL_HARMONIC_COUNT
            ),
            dt_seconds=0.01,
            index_base=0,
        )
    )
    return run_waveform_shape_metric_calculations(
        pack_velocity_per_beat_outputs(result)
    )


class WaveformMetricBeatRecoveryTests(unittest.TestCase):
    def test_recovered_segmentation_matches_oracle_for_every_scalar_metric(self):
        artery_amplitudes = np.ones(7, dtype=np.float32)
        artery_amplitudes[3] = np.float32(0.8)
        artery = _velocity_signal(artery_amplitudes)
        vein = _velocity_signal(
            np.asarray([0.7, 0.8, 0.65, 0.9, 0.75, 0.85, 0.7], dtype=np.float32)
        )
        detection = find_systole_index(artery, dt=np.float32(0.01))
        oracle = np.arange(88, 689, 100, dtype=np.int32)

        np.testing.assert_array_equal(detection.systole_indexes, oracle)
        self.assertNotIn(388, detection.initial_systole_indexes)
        recovered_metrics = _metric_outputs(
            artery,
            vein,
            detection.systole_indexes,
        )
        oracle_metrics = _metric_outputs(artery, vein, oracle)
        legacy_metrics = _metric_outputs(
            artery,
            vein,
            detection.initial_systole_indexes,
        )

        metric_names = set(WaveformShapeMetricsCalculator._metric_names())
        scalar_keys = [
            key
            for key in oracle_metrics
            if key.rsplit("/", maxsplit=1)[-1] in metric_names
            and "/global/" in key
        ]
        self.assertEqual(len(scalar_keys), 4 * len(metric_names))
        for key in scalar_keys:
            np.testing.assert_array_equal(
                np.asarray(recovered_metrics[key]),
                np.asarray(oracle_metrics[key]),
                err_msg=key,
            )

        affected = [
            key
            for key in scalar_keys
            if np.asarray(legacy_metrics[key]).shape
            != np.asarray(oracle_metrics[key]).shape
            or not np.allclose(
                np.asarray(legacy_metrics[key]),
                np.asarray(oracle_metrics[key]),
                equal_nan=True,
            )
        ]
        self.assertGreater(len(affected), len(metric_names))
        for name in ("RI", "PI", "t_min_over_T", "t_rise_over_T", "CF"):
            self.assertTrue(any(key.endswith(f"/{name}") for key in affected), name)


if __name__ == "__main__":
    unittest.main()
