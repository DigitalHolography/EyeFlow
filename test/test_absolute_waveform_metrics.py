"""Tests for the native absolute waveform metric pipeline."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

import pipelines  # noqa: E402
from input_output.schema import EyeFlowOutputPaths  # noqa: E402
from pipeline_engine import PIPELINE_REGISTRY, PipelineDAG  # noqa: E402
from pipelines.absolute_waveform_metrics.calculator import (  # noqa: E402
    AbsoluteWaveformMetricsCalculator,
)
from pipelines.absolute_waveform_metrics.outputs import (  # noqa: E402
    pack_absolute_waveform_outputs,
)
from pipelines.absolute_waveform_metrics.runner import (  # noqa: E402
    run_absolute_waveform_metrics,
)
from pipelines.waveform_velocity_core import runner as core_runner  # noqa: E402


class _State:
    def __init__(self, values=None):
        self.values = dict(values or {})

    def get(self, key, default=None):
        return self.values.get(key, default)

    def set(self, key, value):
        self.values[key] = value


class AbsoluteWaveformMetricsTests(unittest.TestCase):
    def test_pipeline_is_native_package_with_velocity_dependency(self) -> None:
        pipeline_root = SRC_DIR / "pipelines"
        self.assertFalse((pipeline_root / "absolute_waveform_metrics.py").exists())
        self.assertTrue(
            (pipeline_root / "absolute_waveform_metrics" / "calculator.py").exists()
        )

        pipelines.load_pipeline_catalog()
        descriptor = PIPELINE_REGISTRY["absolute_waveform_metrics"]
        self.assertEqual(("waveform_velocity",), descriptor.dag_requires)
        self.assertEqual(
            ("per_beat", "segments", "quadrants"),
            tuple(option.name for option in descriptor.options),
        )
        self.assertEqual(
            "Processing/Metrics/absolute_waveform_metrics",
            EyeFlowOutputPaths.active().absolute_waveform_metrics_root,
        )

        plan = PipelineDAG(PIPELINE_REGISTRY.values()).resolve_targets(
            ["absolute_waveform_metrics"]
        )
        self.assertEqual(
            (
                "waveform_velocity_core",
                "waveform_velocity",
                "absolute_waveform_metrics",
            ),
            plan.names,
        )

    def test_global_outputs_preserve_amplitude_and_prefix_schema_root(self) -> None:
        schema = EyeFlowOutputPaths.active()
        time = np.linspace(0.0, 2.0 * np.pi, 32, endpoint=False)
        waveform = np.stack(
            (
                2.0 + 0.5 * np.sin(time),
                3.0 + 0.25 * np.sin(time),
            ),
            axis=0,
        ).astype(np.float32)
        inputs = {
            schema.beat_period_seconds: np.asarray([[0.8, 0.9]], dtype=np.float32),
            schema.artery_per_beat.velocity_signal: waveform,
            schema.artery_per_beat.velocity_signal_band_limited: waveform,
        }

        outputs = pack_absolute_waveform_outputs(inputs)

        vmax = (
            f"{schema.absolute_waveform_metrics_root}/artery/global/raw/vmax"
        )
        self.assertIn(vmax, outputs)
        np.testing.assert_allclose(outputs[vmax].data, [2.5, 3.25], rtol=1e-5)
        self.assertIn(
            f"{schema.absolute_waveform_metrics_root}/artery/global/"
            "raw_vs_bandlimited/raw_band_corr",
            outputs,
        )

    def test_segment_outputs_keep_native_sample_beat_branch_radius_order(self) -> None:
        schema = EyeFlowOutputPaths.active()
        waveform = np.ones((2, 8), dtype=np.float32)
        segments = np.ones((8, 2, 3, 2), dtype=np.float32)
        inputs = {
            schema.beat_period_seconds: np.asarray([[0.8, 0.9]], dtype=np.float32),
            schema.artery_per_beat.velocity_signal: waveform,
            schema.artery_per_beat.velocity_signal_band_limited: waveform,
            schema.artery_per_beat.segment_velocity_signal: segments,
            schema.artery_per_beat.segment_velocity_signal_band_limited: segments,
        }

        outputs = pack_absolute_waveform_outputs(inputs, include_segments=True)

        key = (
            f"{schema.absolute_waveform_metrics_root}/artery/by_segment/"
            "raw_segment/vmax"
        )
        self.assertEqual((2, 3, 2), outputs[key].data.shape)

    def test_incomplete_waveform_pair_is_rejected(self) -> None:
        schema = EyeFlowOutputPaths.active()
        inputs = {
            schema.beat_period_seconds: np.asarray([[0.8]], dtype=np.float32),
            schema.artery_per_beat.velocity_signal: np.ones(
                (1, 8),
                dtype=np.float32,
            ),
        }

        with self.assertRaisesRegex(ValueError, "Incomplete global waveform"):
            pack_absolute_waveform_outputs(inputs)

    def test_quadrant_outputs_use_absolute_metric_root_and_regions(self) -> None:
        schema = EyeFlowOutputPaths.active()
        waveform = np.ones((2, 8), dtype=np.float32)
        segment_waveform = np.empty((8, 2, 2, 2), dtype=np.float32)
        segment_waveform[:, :, 0, :] = 1.0
        segment_waveform[:, :, 1, :] = 3.0
        labels = np.zeros((8, 8), dtype=np.int32)
        labels[1, 1] = 1
        labels[1, 6] = 2
        segments = SimpleNamespace(
            branch_ids=np.asarray([1, 2], dtype=np.int32),
            labels=labels,
            segment_center_xy=np.zeros((2, 2, 2), dtype=float),
            velocity=np.zeros((2, 2, 3), dtype=np.float32),
        )
        source_data = SimpleNamespace(optic_disc_center=np.asarray([3.0, 2.0]))
        inputs = {
            schema.beat_period_seconds: np.asarray([[0.8, 0.9]], dtype=np.float32),
            schema.artery_per_beat.velocity_signal: waveform,
            schema.artery_per_beat.velocity_signal_band_limited: waveform,
            schema.artery_per_beat.segment_velocity_signal: segment_waveform,
            schema.artery_per_beat.segment_velocity_signal_band_limited: (
                segment_waveform
            ),
        }

        outputs = pack_absolute_waveform_outputs(
            inputs,
            source_data=source_data,
            artery_segments=segments,
            include_per_beat=False,
            include_segments=False,
            include_quadrants=True,
        )

        north_west = (
            f"{schema.absolute_waveform_metrics_root}/artery/Quadrants/"
            "north_west/global/raw/vmax"
        )
        north_east_branch = (
            f"{schema.absolute_waveform_metrics_root}/artery/Quadrants/"
            "north_east/by_branch/branch_2/raw/vmax"
        )
        self.assertIn(north_west, outputs)
        self.assertIn(north_east_branch, outputs)
        self.assertFalse(any("/by_segment/" in key for key in outputs))
        np.testing.assert_allclose(outputs[north_west].data, [1.0, 1.0])
        np.testing.assert_allclose(
            outputs[north_east_branch].data,
            [3.0, 3.0],
        )
        quadrant_root = (
            f"{schema.absolute_waveform_metrics_root}/artery/Quadrants/"
        )
        self.assertEqual(
            {"north_east", "south_east", "north_west", "south_west"},
            {
                key.removeprefix(quadrant_root).split("/", 1)[0]
                for key in outputs
                if key.startswith(quadrant_root)
            },
        )

    def test_runner_consumes_shared_per_beat_state(self) -> None:
        schema = EyeFlowOutputPaths.active()
        waveform = np.ones((1, 8), dtype=np.float32)
        state = _State(
            {
                core_runner.VELOCITY_PER_BEAT_OUTPUTS_STATE: {
                    schema.beat_period_seconds: np.asarray(
                        [[0.8]],
                        dtype=np.float32,
                    ),
                    schema.artery_per_beat.velocity_signal: waveform,
                    schema.artery_per_beat.velocity_signal_band_limited: waveform,
                }
            }
        )
        ctx = SimpleNamespace(
            state=state,
            options_for=lambda name: (
                frozenset(("per_beat",))
                if name == "absolute_waveform_metrics"
                else frozenset()
            ),
        )

        outputs = run_absolute_waveform_metrics(ctx)

        self.assertTrue(outputs)
        self.assertEqual(
            outputs,
            state.get("absolute_waveform_metric_outputs"),
        )

    def test_constant_waveform_has_expected_absolute_metrics(self) -> None:
        metrics = AbsoluteWaveformMetricsCalculator()._compute_metrics_1d(
            np.ones(8, dtype=float),
            2.0,
        )

        self.assertAlmostEqual(metrics["vmean"], 1.0)
        self.assertAlmostEqual(metrics["vti_total"], 2.0)
        self.assertAlmostEqual(metrics["signal_energy"], 2.0)
        self.assertAlmostEqual(metrics["dc_level"], 1.0)


if __name__ == "__main__":
    unittest.main()
