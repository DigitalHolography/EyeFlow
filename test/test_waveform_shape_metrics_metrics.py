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
from pipelines.waveform_shape_metrics.metrics.quadrants import pack_quadrant_metrics
from pipelines.waveform_velocity.quadrants import pack_quadrant_velocity_outputs
from pipelines.waveform_velocity.continuous import pack_segment_velocity_outputs
from utils.logger import Logger


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
        for pipeline_name in (
            "waveform_velocity",
            "waveform_shape_metrics",
            "absolute_waveform_metrics",
            "lowrank_waveform_decomposition",
        ):
            quadrant_option = next(
                option
                for option in PIPELINE_REGISTRY[pipeline_name].options
                if option.name == "quadrants"
            )
            self.assertEqual("Quadrants", quadrant_option.label)
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

    def test_graphics_support_logs_and_skips_invalid_or_missing_beats(self):
        calculator = WaveformShapeMetricsCalculator()
        waveform = np.column_stack(
            (
                np.ones(8),
                -np.ones(8),
                np.zeros(8),
                np.full(8, np.nan),
            )
        )
        periods = np.asarray([[1.0, 0.0, 1.0, 1.0, 1.0]], dtype=np.float32)
        messages = []
        Logger.configure(on_log=messages.append)

        try:
            result = calculator._compute_graphics_support_block(waveform, periods)
        finally:
            Logger.reset_current()

        self.assertTrue(np.isfinite(result["m0"][0]))
        self.assertTrue(np.all(np.isnan(result["m0"][1:])))
        joined_messages = "\n".join(messages)
        self.assertEqual(3, joined_messages.count("invalid input"))
        self.assertIn("Tbeat=0.0", joined_messages)
        self.assertIn("m0_sum=0.0", joined_messages)
        self.assertIn("waveform column is missing", joined_messages)

    def test_quadrant_metrics_are_nested_below_waveform_shape_outputs(self):
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

        regional = pack_quadrant_metrics(
            metrics,
            source_data,
            segments,
            None,
        )

        north_west = (
            f"{schema.waveform_shape_metrics_root}/artery/Quadrants/"
            "north_west/global/raw/mu_t"
        )
        north_east_branch = (
            f"{schema.waveform_shape_metrics_root}/artery/Quadrants/"
            "north_east/by_branch/branch_2/raw/mu_t"
        )
        self.assertIn(north_west, regional)
        self.assertIn(north_east_branch, regional)
        np.testing.assert_allclose(regional[north_west].data, [1.0, 1.0])
        np.testing.assert_allclose(
            regional[north_east_branch].data,
            [1.0, 1.0],
        )
        quadrant_root = (
            f"{schema.waveform_shape_metrics_root}/artery/Quadrants/"
        )
        self.assertEqual(
            {"north_east", "south_east", "north_west", "south_west"},
            {
                key.removeprefix(quadrant_root).split("/", 1)[0]
                for key in regional
                if key.startswith(quadrant_root)
            },
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
            f"{schema.waveform_shape_metrics_root}/artery/Quadrants/"
            "north_west/global/raw/mu_t",
            outputs,
        )

        self.assertFalse(any(key.startswith("Processing/Velocity/") for key in outputs))
        self.assertFalse(
            any(key.startswith("Processing/VelocityPerBeat/") for key in outputs)
        )

        no_quadrants = pack_waveform_shape_outputs(
            metrics,
            source_data,
            segments,
            None,
            include_quadrants=False,
        )
        self.assertFalse(any("/Quadrants/" in key for key in no_quadrants))

        quadrants_without_segment_outputs = pack_waveform_shape_outputs(
            metrics,
            source_data,
            segments,
            None,
            include_segments=False,
            include_quadrants=True,
        )
        self.assertTrue(
            any("/Quadrants/" in key for key in quadrants_without_segment_outputs)
        )
        self.assertFalse(
            any("/by_segment/" in key for key in quadrants_without_segment_outputs)
        )

        no_metrics = pack_waveform_shape_outputs(
            metrics,
            source_data,
            segments,
            None,
            include_per_beat=False,
            include_quadrants=False,
        )
        self.assertEqual({}, no_metrics)

    def test_quadrant_velocity_outputs_reduce_segment_waveforms(self):
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

        outputs = pack_quadrant_velocity_outputs(
            metrics,
            source_data,
            segments,
            None,
        )
        segment_outputs = pack_segment_velocity_outputs(
            segments,
            None,
            source_data=source_data,
        )
        continuous_segment_path = schema.artery_segments.velocity_signal
        self.assertIsNotNone(continuous_segment_path)
        self.assertIn(continuous_segment_path, segment_outputs)
        np.testing.assert_allclose(
            segment_outputs[continuous_segment_path][0],
            segment_velocity.transpose(2, 1, 0),
        )
        continuous_band_limited_path = (
            schema.artery_segments.velocity_signal_band_limited
        )
        self.assertIsNotNone(continuous_band_limited_path)
        self.assertIn(continuous_band_limited_path, segment_outputs)
        np.testing.assert_allclose(
            segment_outputs[continuous_band_limited_path][0],
            segment_velocity.transpose(2, 1, 0),
        )

        north_west_raw = (
            "Processing/Velocity/Quadrants/Artery/north_west/Raw/value"
        )
        north_east_raw = (
            "Processing/Velocity/Quadrants/Artery/north_east/Raw/value"
        )
        north_west_per_beat = (
            "Processing/VelocityPerBeat/Artery/Quadrants/"
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
        continuous_quadrant_root = "Processing/Velocity/Quadrants/Artery/"
        self.assertEqual(
            {"north_east", "south_east", "north_west", "south_west"},
            {
                key.removeprefix(continuous_quadrant_root).split("/", 1)[0]
                for key in outputs
                if key.startswith(continuous_quadrant_root)
            },
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
