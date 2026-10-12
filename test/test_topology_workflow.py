"""Tests for the ordered topology workflow."""

from __future__ import annotations

import unittest
from unittest.mock import patch

import numpy as np

from calculations.topology.geometry import AnnulusGeometry
from calculations.topology.optic_disc import OpticDisc
from calculations.vessel_segments.sampling import prepare_segments, prepare_sampling_plan


class TopologyWorkflowTests(unittest.TestCase):
    def test_both_transform_workflows_use_the_standard_topology(self) -> None:
        from calculations.vessel_segments.sampling.extraction import extract_segment
        from calculations.vessel_segments.sampling.transforms import (
            interpolate_segments,
            resample_rotate_segment,
            rotate_segments,
        )
        vessel = np.zeros((41, 41), dtype=bool)
        vessel[18:23, 5:36] = True
        disc = np.zeros_like(vessel)
        disc[18:23, 18:23] = True
        prepared = prepare_sampling_plan(
            vessel, OpticDisc(disc, (20.0, 20.0), None, None),
            AnnulusGeometry(.1, .6, .25, 2)
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
        prepared = prepare_sampling_plan(
            vessel, OpticDisc(vessel, (10.0, 10.0), None, None),
            AnnulusGeometry(0., .5, .5, 1)
        )
        self.assertEqual([], list(prepare_segments(
            np.zeros((2, 21, 21), np.float32), prepared, worker_count=2,
        )))

    def test_reuses_one_topology_for_multiple_maps(self) -> None:
        vessel = np.zeros((41, 41), dtype=bool)
        vessel[18:23, 5:36] = True
        optic_disc = np.zeros_like(vessel)
        optic_disc[18:23, 18:23] = True
        prepared = prepare_sampling_plan(
            vessel,
            OpticDisc(optic_disc, (20.0, 20.0), None, None),
            AnnulusGeometry(0.1, 0.6, 0.25, 2),
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

    def test_velocity_and_generic_profiles_share_authoritative_topology(self) -> None:
        from calculations.vessel_segments.measurement import (
            SegmentMeasurementSettings,
            analyze_segment_profiles,
        )
        from pipelines.velocity_analysis.analysis.segments import (
            analyze_velocity_segment_profiles,
        )

        vessel = np.zeros((31, 31), dtype=bool)
        vessel[13:18, 3:28] = True
        optic_disc = OpticDisc(None, (15.0, 15.0), 6.0, 6.0)
        rings = AnnulusGeometry(0.1, 0.6, 0.25, 2)
        prepared = prepare_sampling_plan(
            vessel,
            optic_disc,
            rings,
            window_size_percentile_kept=1.0,
        )
        topologies = {"artery": prepared}
        signal = np.ones((3, 31, 31), dtype=np.float32)

        with patch(
            "calculations.vessel_segments.measurement.runner.resolve_segment_rotations"
        ) as resolve:
            generic = analyze_segment_profiles(
                signal,
                {"artery": vessel},
                optic_disc,
                rings,
                SegmentMeasurementSettings(0.01),
                prepared_topologies=topologies,
            )["artery"]
            velocity = analyze_velocity_segment_profiles(
                signal,
                {"artery": vessel},
                optic_disc,
                rings,
                SegmentMeasurementSettings(0.01),
                prepared_topologies=topologies,
            )["artery"]

        resolve.assert_not_called()
        self.assertIs(generic.topology, prepared)
        self.assertIs(velocity.profile.topology, prepared)
        np.testing.assert_array_equal(
            generic.topology.native.branch_ids,
            velocity.profile.topology.native.branch_ids,
        )
        np.testing.assert_allclose(
            generic.topology.native.segment_centers_xy,
            velocity.profile.topology.native.segment_centers_xy,
            equal_nan=True,
        )
        np.testing.assert_allclose(
            generic.topology.rotation_degrees,
            velocity.profile.topology.rotation_degrees,
            equal_nan=True,
        )


if __name__ == "__main__":
    unittest.main()
