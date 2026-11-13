"""Tests for waveform-shape metric calculation over packed EyeFlow outputs."""

import unittest
from types import SimpleNamespace

import numpy as np

from calculations.topology import BranchIdentityResult, OpticDisc, SegmentTopology
from calculations.vessel_segments.sampling import SegmentSamplingPlan
from input_output.schema import EyeFlowOutputPaths
from pipelines.waveform_shape_metrics.metrics.calculator import (
    WaveformShapeMetricsCalculator,
)
from pipelines.waveform_shape_metrics.outputs import pack_waveform_shape_outputs
from utils.logger import Logger


def _segment_topology(labels, branch_ids, centers, optic_disc_center):
    radius_count, branch_count = centers.shape[:2]
    native = SegmentTopology(
        optic_disc_center_xy=optic_disc_center,
        branches=BranchIdentityResult(
            labels,
            branch_ids,
            np.zeros(labels.shape, dtype=bool),
        ),
        annulus_masks=np.zeros((radius_count, *labels.shape), dtype=bool),
        segment_masks=np.zeros((radius_count, branch_count, 1, 1), dtype=bool),
        segment_centers_xy=centers,
        window_bounds_xyxy=np.zeros((radius_count, branch_count, 4), dtype=int),
    )
    return SegmentSamplingPlan(
        native=native,
        rotation_degrees=np.zeros((radius_count, branch_count), dtype=np.float32),
        interpolated_masks=native.segment_masks,
        rotated_masks=native.segment_masks,
    )


class WaveformShapeMetricsTests(unittest.TestCase):
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
        periods = np.asarray([1.0, 1.0], dtype=np.float32)

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
        periods = np.asarray([1.0, 0.0, 1.0, 1.0, 1.0], dtype=np.float32)
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

    def test_output_composer_contains_only_selected_metric_groups(self):
        schema = EyeFlowOutputPaths.active()
        metrics = self._global_artery_inputs(schema)
        waveform = np.asarray(
            metrics[schema.artery_per_beat_filtered.velocity_signal],
            dtype=np.float32,
        )
        segment_waveform = np.broadcast_to(
            waveform[:, :, np.newaxis, np.newaxis],
            (waveform.shape[0], waveform.shape[1], 2, 2),
        ).copy()
        metrics[schema.artery_per_beat_filtered.segment_velocity_signal] = segment_waveform
        metrics[schema.artery_per_beat_filtered.segment_velocity_signal_band_limited] = (
            segment_waveform
        )

        labels = np.zeros((8, 8), dtype=np.int32)
        labels[1, 1] = 1
        labels[1, 6] = 2
        segments = SimpleNamespace(
            profile=SimpleNamespace(
                segment_signal=np.zeros((2, 2, waveform.shape[0])),
                topology=_segment_topology(
                    labels,
                    np.asarray([1, 2], dtype=np.int32),
                    np.zeros((2, 2, 2)),
                    (3.0, 2.0),
                ),
            ),
        )
        source_data = SimpleNamespace(
            retinal_artery_mask=np.zeros((8, 8), dtype=bool),
            retinal_vein_mask=np.zeros((8, 8), dtype=bool),
            optic_disc=OpticDisc(
                np.zeros((8, 8), dtype=bool), (3.0, 2.0), None, None
            ),
            source=SimpleNamespace(
                holodoppler=SimpleNamespace(
                    timing=SimpleNamespace(dt_seconds=np.float32(0.01)),
                ),
            ),
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

        self.assertFalse(any(key.startswith("ProcessingMoments/Velocity/") for key in outputs))
        self.assertFalse(
            any(key.startswith("ProcessingMoments/VelocityPerBeat/") for key in outputs)
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
            schema.cardiac_cycle.systolic_cycle_duration_seconds: (
                np.asarray([0.8, 0.9], dtype=np.float32),
                {"unit": "s"},
            ),
            schema.artery_per_beat_filtered.velocity_signal: waveform,
            schema.artery_per_beat_filtered.velocity_signal_band_limited: (
                waveform,
                {"unit": "mm/s"},
            ),
        }


if __name__ == "__main__":
    unittest.main()
