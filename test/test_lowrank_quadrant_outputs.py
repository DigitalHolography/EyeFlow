"""Tests for quadrant-indexed low-rank waveform decomposition outputs."""

from __future__ import annotations

import unittest
from types import SimpleNamespace

import numpy as np

from input_output.schema import EyeFlowOutputPaths
from pipeline_engine import DatasetValue
from pipelines.lowrank_waveform_decomposition.outputs import (
    pack_lowrank_waveform_decomposition_outputs,
)


class LowRankQuadrantOutputTests(unittest.TestCase):
    def test_outputs_include_four_independently_computed_quadrants(self) -> None:
        schema = EyeFlowOutputPaths.active()
        sample_count = 16
        beat_count = 2
        branch_count = 4
        radius_count = 3
        phase = np.linspace(
            0.0,
            2.0 * np.pi,
            sample_count,
            endpoint=False,
            dtype=np.float32,
        )
        amplitudes = np.arange(1, branch_count + 1, dtype=np.float32)
        waveforms = np.broadcast_to(
            np.sin(phase)[:, None, None, None]
            * amplitudes[None, None, :, None],
            (sample_count, beat_count, branch_count, radius_count),
        ).copy()
        labels = np.zeros((8, 8), dtype=np.int32)
        labels[1, 1] = 1
        labels[1, 6] = 2
        labels[6, 1] = 3
        labels[6, 6] = 4
        segments = SimpleNamespace(
            branch_ids=np.arange(1, branch_count + 1, dtype=np.int32),
            labels=labels,
            segment_center_xy=np.zeros(
                (branch_count, radius_count, 2),
                dtype=float,
            ),
        )
        velocity_outputs = {
            schema.beat_period_seconds: np.asarray([[0.8, 0.9]], dtype=np.float32),
            schema.artery_per_beat.segment_velocity_signal: waveforms,
        }

        outputs = pack_lowrank_waveform_decomposition_outputs(
            velocity_outputs,
            vein_flag=False,
            include_quadrants=True,
            source_data=SimpleNamespace(
                optic_disc_center=np.asarray([3.0, 3.0]),
            ),
            artery_segments=segments,
        )

        root = schema.lowrank_waveform_decomposition_root
        global_r0 = f"{root}/artery/raw/endpoints/joint/R0"
        self.assertIn(global_r0, outputs)
        quadrant_root = f"{root}/artery/Quadrants/"
        self.assertEqual(
            {"north_east", "south_east", "north_west", "south_west"},
            {
                key.removeprefix(quadrant_root).split("/", 1)[0]
                for key in outputs
                if key.startswith(quadrant_root)
            },
        )
        quadrant_r0 = {
            region: self._data(
                outputs[
                    f"{quadrant_root}{region}/raw/endpoints/joint/R0"
                ]
            )
            for region in ("north_west", "north_east", "south_west", "south_east")
        }
        self.assertTrue(all(np.isfinite(value) for value in quadrant_r0.values()))
        self.assertEqual(4, len({float(value) for value in quadrant_r0.values()}))

    def test_quadrant_outputs_are_omitted_when_option_is_disabled(self) -> None:
        schema = EyeFlowOutputPaths.active()
        sample_count = 8
        waveforms = np.ones((sample_count, 1, 1, 3), dtype=np.float32)
        velocity_outputs = {
            schema.beat_period_seconds: np.asarray([[0.8]], dtype=np.float32),
            schema.artery_per_beat.segment_velocity_signal: waveforms,
        }

        outputs = pack_lowrank_waveform_decomposition_outputs(
            velocity_outputs,
            vein_flag=False,
            include_quadrants=False,
        )

        self.assertTrue(outputs)
        self.assertFalse(any("/Quadrants/" in key for key in outputs))

    @staticmethod
    def _data(value) -> np.ndarray:
        return np.asarray(value.data if isinstance(value, DatasetValue) else value)


if __name__ == "__main__":
    unittest.main()
