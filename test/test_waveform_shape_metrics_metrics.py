"""Tests for waveform-shape metric calculation over packed EyeFlow outputs."""

import unittest
from types import SimpleNamespace

import numpy as np

import pipelines  # noqa: F401
from input_output.schema import EyeFlowOutputPaths
from pipeline_engine import PIPELINE_REGISTRY, PipelineDAG
from pipelines.waveform_shape_metrics.metrics.runner import (
    run_waveform_shape_metric_calculations,
)
from pipelines.waveform_shape_metrics.metrics.calculator import (
    WaveformShapeMetricsCalculator,
)
from pipelines.waveform_shape_metrics.outputs import pack_waveform_shape_outputs
from pipelines.waveform_shape_metrics.metrics.hemifield import pack_hemifield_metrics
from pipelines.waveform_velocity.hemifield import pack_hemifield_velocity_outputs


class WaveformShapeMetricsTests(unittest.TestCase):
    def test_waveform_pipelines_have_separate_dag_responsibilities(self):
        pipelines.load_pipeline_catalog()

        self.assertIn("waveform_velocity_core", PIPELINE_REGISTRY)
        self.assertIn("waveform_velocity", PIPELINE_REGISTRY)
        self.assertIn("waveform_shape_metrics", PIPELINE_REGISTRY)
        self.assertIn("pdf_report", PIPELINE_REGISTRY)
        self.assertEqual(
            "hidden",
            PIPELINE_REGISTRY["waveform_velocity_core"].visibility,
        )
        self.assertNotIn("waveform_shape_metrics_angioeye", PIPELINE_REGISTRY)
        self.assertNotIn("topological_metrics", PIPELINE_REGISTRY)
        self.assertEqual(
            "Processing/Metrics/waveform_shape_metrics",
            EyeFlowOutputPaths.active().waveform_shape_metrics_root,
        )

        dag = PipelineDAG(PIPELINE_REGISTRY.values())
        metrics_plan = dag.resolve_targets(["waveform_shape_metrics"])
        report_plan = dag.resolve_targets(["pdf_report"])

        self.assertEqual(
            (
                "waveform_velocity_core",
                "waveform_velocity",
                "waveform_shape_metrics",
            ),
            metrics_plan.names,
        )
        self.assertEqual("pdf_report", report_plan.names[-1])
        self.assertLess(
            report_plan.names.index("waveform_velocity_core"),
            report_plan.names.index("waveform_velocity"),
        )
        self.assertLess(
            report_plan.names.index("waveform_velocity_core"),
            report_plan.names.index("waveform_shape_metrics"),
        )
        self.assertNotIn("waveform_shape_metrics_angioeye", metrics_plan.names)

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

    def test_ri_pi_use_unrectified_mean_repeating_cycle(self):
        cycles = np.asarray(
            [
                [-2.0, 2.0],
                [4.0, 2.0],
                [0.0, 2.0],
                [2.0, 2.0],
            ],
            dtype=np.float32,
        )
        periods = np.asarray([[1.0, 1.0]], dtype=np.float32)

        result = WaveformShapeMetricsCalculator()._compute_block_global(
            cycles,
            periods,
        )

        # The unrectified mean cycle is [0, 3, 1, 2].
        np.testing.assert_allclose(result["RI"], [1.0, 1.0])
        np.testing.assert_allclose(result["PI"], [2.0, 2.0])

    def test_hemifield_metrics_are_nested_below_waveform_shape_outputs(self):
        schema = EyeFlowOutputPaths.active()
        sample_count = 2
        branch_count = 2
        radius_count = 2
        metric_values = np.ones(
            (sample_count, branch_count, radius_count),
            dtype=np.float32,
        )
        metrics = {}
        for signal_group, signal_name in (
            ("raw_segment", "raw"),
            ("bandlimited_segment", "bandlimited"),
        ):
            for metric_name in WaveformShapeMetricsCalculator._metric_names():
                metrics[
                    f"{schema.waveform_shape_metrics_root}/artery/"
                    f"by_segment/{signal_group}/{metric_name}"
                ] = (metric_values, {"signal_type": signal_name})

        labels = np.zeros((8, 8), dtype=np.int32)
        labels[1, 1] = 1
        labels[1, 6] = 2
        segments = SimpleNamespace(
            branch_ids=np.asarray([1, 2], dtype=np.int32),
            labels=labels,
            segment_center_xy=np.zeros((branch_count, radius_count, 2)),
            velocity=np.zeros((radius_count, branch_count, 3)),
        )
        source_data = SimpleNamespace(
            retinal_artery_mask=np.zeros((8, 8), dtype=bool),
            retinal_vein_mask=np.zeros((8, 8), dtype=bool),
            optic_disc_mask=np.zeros((8, 8), dtype=bool),
            optic_disc_center=np.asarray([3.0, 2.0]),
            optic_disc_width=None,
            optic_disc_height=None,
        )

        regional = pack_hemifield_metrics(
            metrics,
            source_data,
            segments,
            None,
        )

        north_west = (
            f"{schema.waveform_shape_metrics_root}/artery/hemifield/"
            "north_west/global/raw/mu_t"
        )
        north_east_branch = (
            f"{schema.waveform_shape_metrics_root}/artery/hemifield/"
            "north_east/by_branch/branch_2/raw/mu_t"
        )
        self.assertIn(north_west, regional)
        self.assertIn(north_east_branch, regional)
        np.testing.assert_allclose(regional[north_west].data, [1.0, 1.0])
        np.testing.assert_allclose(
            regional[north_east_branch].data,
            [1.0, 1.0],
        )

    def test_output_composer_contains_only_selected_metric_groups(self):
        schema = EyeFlowOutputPaths.active()
        metrics = self._global_artery_inputs(schema)
        waveform = np.asarray(
            metrics[schema.artery_per_beat.velocity_signal],
            dtype=np.float32,
        )
        segment_waveform = np.broadcast_to(
            waveform[:, :, np.newaxis, np.newaxis],
            (waveform.shape[0], waveform.shape[1], 2, 2),
        ).copy()
        metrics[schema.artery_per_beat.segment_velocity_signal] = segment_waveform
        metrics[schema.artery_per_beat.segment_velocity_signal_band_limited] = (
            segment_waveform
        )

        labels = np.zeros((8, 8), dtype=np.int32)
        labels[1, 1] = 1
        labels[1, 6] = 2
        segments = SimpleNamespace(
            branch_ids=np.asarray([1, 2], dtype=np.int32),
            labels=labels,
            segment_center_xy=np.zeros((2, 2, 2)),
            velocity=np.zeros((2, 2, waveform.shape[0])),
        )
        source_data = SimpleNamespace(
            retinal_artery_mask=np.zeros((8, 8), dtype=bool),
            retinal_vein_mask=np.zeros((8, 8), dtype=bool),
            optic_disc_mask=np.zeros((8, 8), dtype=bool),
            optic_disc_center=np.asarray([3.0, 2.0]),
            optic_disc_width=None,
            optic_disc_height=None,
            timing=SimpleNamespace(dt_seconds=np.float32(0.01)),
        )

        outputs = pack_waveform_shape_outputs(
            metrics,
            source_data,
            segments,
            None,
        )

        self.assertNotIn(schema.segmentation.artery.branch_label_map, outputs)
        self.assertIn(
            f"{schema.waveform_shape_metrics_root}/artery/global/raw/mu_t",
            outputs,
        )
        self.assertIn(
            f"{schema.waveform_shape_metrics_root}/artery/"
            "by_segment/raw_segment/mu_t",
            outputs,
        )
        self.assertIn(
            f"{schema.waveform_shape_metrics_root}/artery/hemifield/"
            "north_west/global/raw/mu_t",
            outputs,
        )

        self.assertFalse(any(key.startswith("Processing/Velocity/") for key in outputs))
        self.assertFalse(
            any(key.startswith("Processing/VelocityPerBeat/") for key in outputs)
        )

        no_hemifield = pack_waveform_shape_outputs(
            metrics,
            source_data,
            segments,
            None,
            include_hemifield=False,
        )
        self.assertFalse(any("/hemifield/" in key for key in no_hemifield))

        hemifield_without_segment_outputs = pack_waveform_shape_outputs(
            metrics,
            source_data,
            segments,
            None,
            include_segments=False,
            include_hemifield=True,
        )
        self.assertTrue(
            any("/hemifield/" in key for key in hemifield_without_segment_outputs)
        )
        self.assertFalse(
            any("/by_segment/" in key for key in hemifield_without_segment_outputs)
        )

        no_metrics = pack_waveform_shape_outputs(
            metrics,
            source_data,
            segments,
            None,
            include_per_beat=False,
            include_hemifield=False,
        )
        self.assertEqual({}, no_metrics)

    def test_hemifield_velocity_outputs_reduce_segment_waveforms(self):
        schema = EyeFlowOutputPaths.active()
        labels = np.zeros((8, 8), dtype=np.int32)
        labels[1, 1] = 1
        labels[1, 6] = 2
        segment_velocity = np.asarray(
            [
                [[1.0, 2.0, 3.0], [10.0, 20.0, 30.0]],
                [[3.0, 4.0, 5.0], [30.0, 40.0, 50.0]],
            ],
            dtype=np.float32,
        )
        segment_per_beat = np.asarray(
            [
                [[[1.0, 3.0], [10.0, 30.0]], [[2.0, 4.0], [20.0, 40.0]]],
                [[[5.0, 7.0], [50.0, 70.0]], [[6.0, 8.0], [60.0, 80.0]]],
            ],
            dtype=np.float32,
        )
        segments = SimpleNamespace(
            branch_ids=np.asarray([1, 2], dtype=np.int32),
            labels=labels,
            segment_center_xy=np.zeros((2, 2, 2)),
            velocity=segment_velocity,
        )
        source_data = SimpleNamespace(
            retinal_artery_mask=np.zeros((8, 8), dtype=bool),
            retinal_vein_mask=np.zeros((8, 8), dtype=bool),
            optic_disc_mask=np.zeros((8, 8), dtype=bool),
            optic_disc_center=np.asarray([3.0, 2.0]),
            optic_disc_width=None,
            optic_disc_height=None,
            timing=SimpleNamespace(dt_seconds=np.float32(0.01)),
        )
        metrics = {
            schema.artery_per_beat.segment_velocity_signal: segment_per_beat,
            schema.artery_per_beat.segment_velocity_signal_band_limited: (
                segment_per_beat + 100.0
            ),
        }

        outputs = pack_hemifield_velocity_outputs(
            metrics,
            source_data,
            segments,
            None,
        )

        north_west_raw = (
            "Processing/Velocity/Artery/hemifield/north_west/Raw/value"
        )
        north_east_raw = (
            "Processing/Velocity/Artery/hemifield/north_east/Raw/value"
        )
        north_west_per_beat = (
            "Processing/VelocityPerBeat/Artery/hemifield/"
            "north_west/Raw/value"
        )
        self.assertIn(north_west_raw, outputs)
        self.assertIn(north_east_raw, outputs)
        np.testing.assert_allclose(outputs[north_west_raw].data, [2.0, 3.0, 4.0])
        np.testing.assert_allclose(
            outputs[north_east_raw].data,
            [20.0, 30.0, 40.0],
        )
        np.testing.assert_allclose(
            outputs[north_west_per_beat].data,
            [[2.0, 6.0], [3.0, 7.0]],
        )

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
