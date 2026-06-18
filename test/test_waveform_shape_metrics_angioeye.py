"""Tests for the EyeFlow adapter around AngioEye waveform metrics."""

import unittest

import h5py
import numpy as np

import pipelines  # noqa: F401
from input_output.inputs import EyeFlowView
from input_output.schema import EyeFlowOutputPaths
from pipeline_engine import PIPELINE_REGISTRY, PipelineDAG
from pipelines.waveform_shape_metrics_angioeye.runner import (
    run_waveform_shape_metrics_angioeye,
)


class _Context:
    def __init__(self, h5file):
        self.ef = EyeFlowView(h5file)


class WaveformShapeMetricsAngioEyeTests(unittest.TestCase):
    def test_dag_runs_eyeflow_waveforms_first(self):
        plan = PipelineDAG(PIPELINE_REGISTRY.values()).resolve_targets(
            ["waveform_shape_metrics_angioeye"]
        )

        self.assertEqual(
            plan.names[-2:],
            ("waveform_shape_metrics", "waveform_shape_metrics_angioeye"),
        )

    def test_runner_reads_schema_paths_and_prefixes_outputs(self):
        schema = EyeFlowOutputPaths.active()
        with self._memory_h5() as h5file:
            self._write_global_artery_inputs(h5file, schema)

            metrics = run_waveform_shape_metrics_angioeye(_Context(h5file))

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
        with self._memory_h5() as h5file:
            h5file.create_dataset(
                schema.beat_period_seconds,
                data=np.asarray([[0.8]], dtype=np.float32),
            )
            h5file.create_dataset(
                schema.artery_per_beat.velocity_signal,
                data=np.ones((16, 1), dtype=np.float32),
            )

            with self.assertRaisesRegex(ValueError, "Incomplete global waveform"):
                run_waveform_shape_metrics_angioeye(_Context(h5file))

    @staticmethod
    def _write_global_artery_inputs(h5file, schema):
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
        h5file.create_dataset(
            schema.beat_period_seconds,
            data=np.asarray([[0.8, 0.9]], dtype=np.float32),
        )
        h5file.create_dataset(
            schema.artery_per_beat.velocity_signal,
            data=waveform,
        )
        h5file.create_dataset(
            schema.artery_per_beat.velocity_signal_band_limited,
            data=waveform,
        )

    @staticmethod
    def _memory_h5():
        return h5py.File(
            "angioeye_metrics_test.h5",
            "w",
            driver="core",
            backing_store=False,
        )


if __name__ == "__main__":
    unittest.main()
