from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from calculations.topology import (  # noqa: E402
    AnnulusGeometry,
    BranchIdentityResult,
    OpticDisc,
    SegmentTopology,
    build_segment_topology,
)
from calculations.vessel_segments.sampling.extraction import (  # noqa: E402
    extract_segment,
)


class SegmentTopologyTests(unittest.TestCase):
    def test_builds_map_independent_segment_windows(self) -> None:
        vessel = np.zeros((41, 41), dtype=bool)
        vessel[18:23, 5:36] = True
        optic_disc = np.zeros_like(vessel)
        optic_disc[18:23, 18:23] = True

        topology = build_segment_topology(
            vessel,
            OpticDisc(optic_disc, (20.0, 20.0), None, None),
            AnnulusGeometry(0.1, 0.6, 0.25, 2),
            window_size_percentile_kept=1.0,
        )

        self.assertEqual((20.0, 20.0), topology.optic_disc_center_xy)
        self.assertEqual((2, 2, 2), topology.segment_centers_xy.shape)
        self.assertEqual((2, 2, 7, 7), topology.segment_masks.shape)
        self.assertEqual((41, 41), topology.centerline.shape)
        self.assertEqual(np.bool_, topology.centerline.dtype)
        self.assertEqual(7, topology.window_side_pixels)
        self.assertEqual(np.float32, topology.delta_radius.dtype)
        np.testing.assert_allclose(
            topology.delta_radius,
            np.full(2, np.hypot(20.0, 20.0) * 0.25, dtype=np.float32),
        )
        self.assertTrue(np.all(topology.valid_segments))

    def test_explicit_center_drives_the_centered_circle_exclusion(self) -> None:
        vessel = np.ones((31, 35), dtype=bool)
        disc = np.zeros_like(vessel)
        disc[4:9, 5:12] = True
        settings = AnnulusGeometry(0.0, 0.8, 0.4, 2)
        topology = build_segment_topology(
            vessel,
            OpticDisc(disc, (24.0, 20.0), None, None),
            settings,
        )
        self.assertEqual((24.0, 20.0), topology.optic_disc_center_xy)
        self.assertIs(topology.ring_settings, settings)
        circle = topology.optic_disc_mask
        self.assertFalse(np.any(topology.labels[circle]))
        self.assertFalse(np.any(topology.annulus_masks[:, circle]))
        self.assertTrue(np.any(topology.branch_identity.stages.vessel[circle]))
        self.assertFalse(
            np.any(topology.branch_identity.stages.cleaned_skeleton[circle])
        )
        self.assertTrue(np.any(topology.annulus_masks[:, disc]))
        self.assertTrue(np.any(topology.branch_identity.stages.vessel[disc]))

    def test_explicit_center_is_used_for_nonempty_and_empty_masks(self) -> None:
        vessel = np.zeros((21, 31), dtype=bool)
        disc = np.zeros_like(vessel)
        disc[4:7, 8:11] = True
        settings = AnnulusGeometry(0.0, 0.8, 0.4, 1)
        center = (13.0, 8.0)
        from_mask = build_segment_topology(
            vessel, OpticDisc(disc, center, None, None), settings
        )
        from_empty = build_segment_topology(
            vessel, OpticDisc(np.zeros_like(disc), center, None, None), settings
        )
        self.assertEqual(center, from_mask.optic_disc_center_xy)
        self.assertEqual(center, from_empty.optic_disc_center_xy)

    def test_extracts_vector_maps_with_explicit_spatial_axes(self) -> None:
        topology = _edge_topology()
        vector_map = np.arange(2 * 4 * 5 * 2, dtype=np.float32).reshape(2, 4, 5, 2)

        segment = extract_segment(vector_map, topology, 0, 0, spatial_axes=(1, 2))

        self.assertEqual((2, 2, 3, 3), segment.shape)
        self.assertTrue(np.all(np.isnan(segment[:, :, 0, :])))
        self.assertTrue(np.all(np.isnan(segment[:, :, :, 0])))
        np.testing.assert_array_equal(
            segment[:, :, 1:, 1:],
            np.moveaxis(vector_map[:, :2, :2, :], (1, 2), (-2, -1)),
        )

    def test_extracts_scalar_stacks_without_topology_recalculation(self) -> None:
        topology = _edge_topology()
        scalar_stack = np.arange(2 * 4 * 5, dtype=np.float32).reshape(2, 4, 5)

        first = extract_segment(scalar_stack, topology, 0, 0)
        second = extract_segment(scalar_stack + 100.0, topology, 0, 0)

        self.assertEqual((2, 3, 3), first.shape)
        np.testing.assert_array_equal(
            second[..., 1:, 1:],
            first[..., 1:, 1:] + 100.0,
        )


def _edge_topology() -> SegmentTopology:
    labels = np.ones((4, 5), dtype=np.int32)
    return SegmentTopology(
        optic_disc_center_xy=(2.0, 2.0),
        branches=BranchIdentityResult(
            labels,
            np.asarray([1], dtype=np.int32),
            np.ones((4, 5), dtype=bool),
        ),
        annulus_masks=np.ones((1, 4, 5), dtype=bool),
        segment_masks=np.ones((1, 1, 3, 3), dtype=bool),
        segment_centers_xy=np.asarray([[[0.0, 0.0]]], dtype=np.float32),
        window_bounds_xyxy=np.asarray([[[0, 2, 0, 2]]], dtype=np.int32),
    )


if __name__ == "__main__":
    unittest.main()
