"""Tests for the ordered topology workflow."""

from __future__ import annotations

import unittest

import numpy as np

from calculations.topology.geometry import SegmentRingSettings
from calculations.topology.workflow import prepare_segments, prepare_topology


class TopologyWorkflowTests(unittest.TestCase):
    def test_both_transform_workflows_use_the_standard_topology(self) -> None:
        from calculations.topology import (
            extract_segment, interpolate_segments, rotate_segments,
            resample_rotate_segment,
        )
        vessel = np.zeros((41, 41), dtype=bool)
        vessel[18:23, 5:36] = True
        disc = np.zeros_like(vessel)
        disc[18:23, 18:23] = True
        prepared = prepare_topology(
            vessel, disc, SegmentRingSettings(.1, .6, .25, 2)
        )
        cube = np.arange(2 * 41 * 41, dtype=np.float32).reshape(2, 41, 41)
        for mode in ("fused", "sequential"):
            items = list(prepare_segments(cube, prepared, transform_mode=mode))
            self.assertTrue(items)
            for item in items:
                index = (item.ring_index, item.branch_index)
                extracted = extract_segment(cube, prepared.topology, *index)
                angle = float(prepared.rotation_degrees[index])
                if mode == "fused":
                    expected = resample_rotate_segment(extracted, angle, 128)
                else:
                    expected = rotate_segments(
                        interpolate_segments(extracted, 128)[None, None],
                        np.array([[angle]], np.float32),
                    )[0, 0]
                np.testing.assert_allclose(item.rotated, expected, equal_nan=True)

    def test_empty_topology_allows_parallel_preparation(self) -> None:
        vessel = np.zeros((21, 21), bool)
        prepared = prepare_topology(
            vessel, vessel, SegmentRingSettings(0., .5, .5, 1)
        )
        self.assertEqual([], list(prepare_segments(
            np.zeros((2, 21, 21), np.float32), prepared, worker_count=2,
        )))

    def test_reuses_one_topology_for_multiple_maps(self) -> None:
        vessel = np.zeros((41, 41), dtype=bool)
        vessel[18:23, 5:36] = True
        optic_disc = np.zeros_like(vessel)
        optic_disc[18:23, 18:23] = True
        prepared = prepare_topology(
            vessel,
            optic_disc,
            SegmentRingSettings(0.1, 0.6, 0.25, 2),
            window_size_percentile_kept=1.0,
        )

        first = list(
            prepare_segments(np.ones((2, 41, 41), dtype=np.float32), prepared)
        )
        second = list(
            prepare_segments(
                np.full((2, 41, 41), 2.0, dtype=np.float32),
                prepared,
            )
        )

        ring_count, branch_count = prepared.topology.valid_segments.shape
        expected_count = np.count_nonzero(
            prepared.topology.valid_segments
            & np.isfinite(prepared.rotation_degrees)
        )
        self.assertEqual(expected_count, len(first))
        self.assertEqual(expected_count, len(second))
        self.assertTrue(all(item.rotated.shape == (2, 181, 181) for item in first))
        self.assertEqual(
            [(item.ring_index, item.branch_index) for item in first],
            [(item.ring_index, item.branch_index) for item in second],
        )
        self.assertEqual(
            (ring_count, branch_count, 128, 128),
            prepared.interpolated_masks.shape,
        )
        self.assertEqual(
            (ring_count, branch_count, 181, 181),
            prepared.rotated_masks.shape,
        )
        np.testing.assert_allclose(
            np.concatenate(
                [
                    item.rotated[np.isfinite(item.rotated)]
                    for item in second
                ]
            ),
            2.0,
        )


if __name__ == "__main__":
    unittest.main()
