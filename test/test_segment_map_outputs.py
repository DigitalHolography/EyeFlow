"""Tests for rotated per-segment velocity-map and mask outputs."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import h5py
import numpy as np

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from calculations.math import interpft_real  # noqa: E402
from input_output.schema import EyeFlowOutputPaths  # noqa: E402
from input_output.writers.h5 import write_value_dataset  # noqa: E402
from pipelines.waveform_velocity.segment_maps import (  # noqa: E402
    _segment_map_worker_count,
    interpolate_velocity_maps_per_beat,
    pack_segment_map_outputs,
)


class SegmentMapOutputTests(unittest.TestCase):
    def test_map_interpolation_matches_scalar_interpft(self) -> None:
        maps = np.arange(1 * 1 * 6 * 2 * 3, dtype=np.float32).reshape(
            1,
            1,
            6,
            2,
            3,
        )
        boundaries = np.asarray([0, 2, 5], dtype=np.int32)

        actual = interpolate_velocity_maps_per_beat(maps, boundaries)
        expected = interpft_real(maps[0, 0, 0:3, 1, 2], 5)[:-1]

        self.assertEqual((3, 2, 4, 2, 1, 1), actual.shape)
        np.testing.assert_allclose(actual[2, 1, :, 0, 0, 0], expected, rtol=1e-6)

    def test_threaded_segments_match_scalar_interpft(self) -> None:
        maps = np.arange(2 * 3 * 8 * 2 * 3, dtype=np.float32).reshape(
            2,
            3,
            8,
            2,
            3,
        )
        maps[1, 2, :, 0, 1] = np.nan
        boundaries = np.asarray([0, 3, 7], dtype=np.int32)

        actual = interpolate_velocity_maps_per_beat(maps, boundaries)

        self.assertEqual((3, 2, 4, 2, 3, 2), actual.shape)
        for radius_index in range(2):
            for branch_index in range(3):
                for beat_index, (start, stop) in enumerate(((0, 4), (3, 8))):
                    for y_index in range(2):
                        for x_index in range(3):
                            expected = interpft_real(
                                maps[
                                    radius_index,
                                    branch_index,
                                    start:stop,
                                    y_index,
                                    x_index,
                                ],
                                5,
                            )[:-1]
                            np.testing.assert_allclose(
                                actual[
                                    x_index,
                                    y_index,
                                    :,
                                    beat_index,
                                    branch_index,
                                    radius_index,
                                ],
                                expected,
                                rtol=1e-6,
                                equal_nan=True,
                            )

    def test_segment_worker_count_honors_parallel_job_cap(self) -> None:
        with patch(
            "pipelines.waveform_velocity.segment_maps.cap_parallel_jobs",
            return_value=4,
        ) as capped:
            self.assertEqual(1, _segment_map_worker_count(0))
            self.assertEqual(1, _segment_map_worker_count(1))
            self.assertEqual(3, _segment_map_worker_count(3))
            self.assertEqual(4, _segment_map_worker_count(20))

        capped.assert_called_with(8)

    def test_h5_outputs_have_requested_paths_shapes_and_types(self) -> None:
        artery = _segments(radius_count=2, branch_count=1)
        vein = _segments(radius_count=2, branch_count=0)
        metrics = pack_segment_map_outputs(
            artery,
            vein,
            np.asarray([0, 2, 5], dtype=np.int32),
        )
        schema = EyeFlowOutputPaths.active()

        self.assertEqual(4, len(metrics))
        self.assertEqual(
            "Processing/VelocityMapPerSegment/Artery",
            schema.artery_segments.velocity_map_per_segment,
        )
        self.assertEqual(
            "Processing/VelocityMapPerSegment/Vein",
            schema.vein_segments.velocity_map_per_segment,
        )
        self.assertEqual(
            "Segmentation/Artery/Segments",
            schema.artery_segments.segments,
        )
        self.assertEqual(
            "Segmentation/Vein/Segments",
            schema.vein_segments.segments,
        )

        with h5py.File("segment_maps.h5", "w", driver="core", backing_store=False) as h5:
            for path, value in metrics.items():
                write_value_dataset(h5, path, value)

            artery_maps = h5[schema.artery_segments.velocity_map_per_segment]
            artery_masks = h5[schema.artery_segments.segments]
            vein_maps = h5[schema.vein_segments.velocity_map_per_segment]
            vein_masks = h5[schema.vein_segments.segments]

            self.assertEqual((4, 3, 4, 2, 1, 2), artery_maps.shape)
            self.assertEqual((4, 3, 1, 2), artery_masks.shape)
            self.assertEqual((4, 3, 4, 2, 0, 2), vein_maps.shape)
            self.assertEqual((4, 3, 0, 2), vein_masks.shape)
            self.assertEqual(
                ["x", "y", "time", "beat", "branch", "radius"],
                list(artery_maps.attrs["dimDesc"]),
            )
            self.assertEqual(
                ["x", "y", "branch", "radius"],
                list(artery_masks.attrs["dimDesc"]),
            )
            self.assertEqual(np.bool_, artery_masks.dtype)
            self.assertEqual("lzf", artery_maps.compression)
            self.assertEqual((4, 3, 1, 1, 1, 1), artery_maps.chunks)


def _segments(*, radius_count: int, branch_count: int):
    maps = np.arange(
        radius_count * branch_count * 6 * 3 * 4,
        dtype=np.float32,
    ).reshape(radius_count, branch_count, 6, 3, 4)
    masks = np.zeros((radius_count, branch_count, 3, 4), dtype=bool)
    masks[..., 1:, 1:3] = True
    return SimpleNamespace(
        velocity_maps_per_segment=maps,
        segment_masks=masks,
    )


if __name__ == "__main__":
    unittest.main()
