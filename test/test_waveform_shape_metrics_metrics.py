"""Tests for waveform-shape metric calculation over packed EyeFlow outputs."""

import unittest

import numpy as np

import pipelines  # noqa: F401
from input_output.schema import EyeFlowOutputPaths
from pipeline_engine import PIPELINE_REGISTRY, PipelineDAG
from pipelines.waveform_shape_metrics.metrics.runner import (
    run_waveform_shape_metric_calculations,
)


class WaveformShapeMetricsTests(unittest.TestCase):
    def test_waveform_shape_metrics_is_the_only_registered_pipeline(self):
        self.assertIn("waveform_shape_metrics", PIPELINE_REGISTRY)
        self.assertNotIn("waveform_shape_metrics_angioeye", PIPELINE_REGISTRY)

        plan = PipelineDAG(PIPELINE_REGISTRY.values()).resolve_targets(
            ["waveform_shape_metrics"]
        )

        self.assertEqual(plan.names[-1], "waveform_shape_metrics")
        self.assertNotIn("waveform_shape_metrics_angioeye", plan.names)

    def test_runner_reads_packed_metrics_and_prefixes_outputs(self):
        schema = EyeFlowOutputPaths.active()
        packed_metrics = self._global_artery_inputs(schema)

        metrics = run_waveform_shape_metric_calculations(packed_metrics)

        self.assertTrue(metrics)
        self.assertTrue(
            all(
                key.startswith(f"{schema.waveform_shape_metrics_root}/artery/")
                for key in metrics
            )
        )
        self.assertFalse(any("/vein/" in key for key in metrics))

    def test_runner_rejects_an_incomplete_waveform_pair(self):
        schema = EyeFlowOutputPaths.active()
        packed_metrics = {
            schema.beat_period_seconds: (
                np.asarray([[0.8]], dtype=np.float32),
                {"unit": "s"},
            ),
            schema.artery_per_beat.velocity_signal: np.ones(
                (16, 1),
                dtype=np.float32,
            ),
        }

        with self.assertRaisesRegex(ValueError, "Incomplete global waveform"):
            run_waveform_shape_metric_calculations(packed_metrics)

    @staticmethod
    def _global_artery_inputs(schema):
        sample_count = 32
        time = np.linspace(
            0,
            2 * np.pi,
            sample_count,
            endpoint=False,
            dtype=np.float32,
        )
        waveform = np.stack(
            (1 + 0.4 * np.sin(time), 1.2 + 0.3 * np.sin(time + 0.2)),
            axis=1,
        ).astype(np.float32)
        return {
            schema.beat_period_seconds: (
                np.asarray([[0.8, 0.9]], dtype=np.float32),
                {"unit": "s"},
            ),
            schema.artery_per_beat.velocity_signal: waveform,
            schema.artery_per_beat.velocity_signal_band_limited: (
                waveform,
                {"unit": "mm/s"},
            ),
        }


if __name__ == "__main__":
    unittest.main()
