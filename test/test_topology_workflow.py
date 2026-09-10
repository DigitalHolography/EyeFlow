"""Tests for the ordered topology workflow."""

from __future__ import annotations

import unittest

import numpy as np

from calculations.topology.geometry import SegmentRingSettings
from calculations.topology.workflow import prepare_segments, prepare_topology


class TopologyWorkflowTests(unittest.TestCase):
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
