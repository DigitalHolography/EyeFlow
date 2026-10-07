"""Tests for topology-owned anatomical quadrant assignment."""

from __future__ import annotations

import unittest

import numpy as np

from calculations.topology import (
    QUADRANT_NAMES,
    BranchIdentityResult,
    PreparedTopology,
    SegmentTopology,
    quadrant_membership,
)


class TopologyQuadrantTests(unittest.TestCase):
    def test_segment_and_prepared_topology_share_quadrant_assignment(self) -> None:
        topology = _four_quadrant_topology(radius_count=3)
        prepared = PreparedTopology(
            native=topology,
            rotation_degrees=np.zeros((3, 4), dtype=np.float32),
            interpolated_masks=np.zeros((3, 4, 1, 1), dtype=bool),
            rotated_masks=np.zeros((3, 4, 1, 1), dtype=bool),
        )

        direct = quadrant_membership(topology)
        through_prepared = quadrant_membership(prepared)

        self.assertEqual(
            ("north_west", "north_east", "south_west", "south_east"),
            QUADRANT_NAMES,
        )
        self.assertEqual((4, 4, 3), direct.shape)
        np.testing.assert_array_equal(direct, through_prepared)
        np.testing.assert_array_equal(direct[:, :, 0], np.eye(4, dtype=bool))

    def test_every_branch_must_exist_in_the_label_map(self) -> None:
        topology = _four_quadrant_topology(radius_count=1)
        topology.labels[topology.labels == 4] = 0

        with self.assertRaisesRegex(ValueError, "at least one pixel"):
            quadrant_membership(topology)


def _four_quadrant_topology(*, radius_count: int) -> SegmentTopology:
    labels = np.zeros((8, 8), dtype=np.int32)
    labels[1, 1] = 1
    labels[1, 6] = 2
    labels[6, 1] = 3
    labels[6, 6] = 4
    branch_count = 4
    return SegmentTopology(
        optic_disc_center_xy=(3.0, 3.0),
        branches=BranchIdentityResult(
            labels,
            np.arange(1, branch_count + 1, dtype=np.int32),
            labels > 0,
        ),
        annulus_masks=np.zeros((radius_count, *labels.shape), dtype=bool),
        segment_masks=np.zeros((radius_count, branch_count, 1, 1), dtype=bool),
        segment_centers_xy=np.zeros(
            (radius_count, branch_count, 2),
            dtype=np.float32,
        ),
        window_bounds_xyxy=np.zeros(
            (radius_count, branch_count, 4),
            dtype=np.int32,
        ),
    )


if __name__ == "__main__":
    unittest.main()
