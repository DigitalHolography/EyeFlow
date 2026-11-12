"""Tests for the ordered topology workflow."""

from __future__ import annotations

import unittest
from unittest.mock import patch

import numpy as np

from calculations.topology.geometry import AnnulusGeometry
from calculations.topology.optic_disc import OpticDisc
from calculations.vessel_segments.sampling import prepare_segment_chunks, prepare_sampling_plan


class TopologyWorkflowTests(unittest.TestCase):
    def test_both_transform_workflows_use_the_standard_topology(self) -> None:
        from calculations.vessel_segments.sampling.extraction import extract_segment
        from calculations.vessel_segments.sampling.transforms import (
            interpolate_segments,
            resample_rotate_segment,
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
            items = list(prepare_segment_chunks(cube, prepared, transform_mode=mode))
            self.assertTrue(items)
            for item in items:
                index = (item.ring_index, item.branch_index)
                extracted = extract_segment(cube, prepared.native, *index)
                angle = float(prepared.rotation_degrees[index])
                if mode == "fused":
                    expected = resample_rotate_segment(extracted, angle, 128)
                else:
                    expected = resample_rotate_segment(
                        interpolate_segments(extracted, 128), angle, 128,
                    )
                np.testing.assert_allclose(item.rotated, expected, equal_nan=True)

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
