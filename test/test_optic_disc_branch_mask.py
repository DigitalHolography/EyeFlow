"""Regression tests for branch extraction at the actual optic-disc boundary."""

import sys
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from calculations.topology import (
    AnnulusGeometry,
    OpticDisc,
    build_segment_topology,
    label_vessel_branches,
)
from input_output.schema import PixelPitch
from pipelines.topology_core.outputs import pack_topology_outputs


class OpticDiscBranchMaskTests(unittest.TestCase):
    def test_near_disc_branch_points_use_the_complete_vessel_context(self):
        vessel = np.zeros((61, 61), dtype=bool)
        vessel[30, 5:56] = True
        vessel[5:31, 39] = True
        optic_disc = OpticDisc(None, (30.0, 30.0), 16.0, 16.0)

        branches = label_vessel_branches(
            vessel,
            optic_disc,
            AnnulusGeometry(0.0, 1.0, 0.2, 4),
        )
        circle = optic_disc.centered_circle_mask_for(vessel.shape)

        np.testing.assert_array_equal(branches.stages.vessel, vessel)
        self.assertTrue(np.any(branches.stages.skeleton[circle]))
        self.assertTrue(branches.stages.branch_points[30, 39])
        self.assertFalse(np.any(branches.stages.cleaned_skeleton[circle]))
        upper_marker = branches.stages.marker_labels[15, 39]
        right_marker = branches.stages.marker_labels[30, 50]
        self.assertGreater(upper_marker, 0)
        self.assertGreater(right_marker, 0)
        self.assertNotEqual(upper_marker, right_marker)

    def test_actual_mask_does_not_replace_the_centered_topology_circle(self):
        vessel, disc, source, settings = self._inputs()
        # A nonelliptical extension reaches beyond the nominal disc dimensions.
        disc[40:49, 90:97] = True
        topology = build_segment_topology(
            vessel,
            OpticDisc(disc, source.source.segmentation.optic_disc.center, 20.0, 50.0),
            settings,
        )

        topology_disc = topology.optic_disc_mask
        self.assertTrue(np.any(topology.branch_identity.stages.vessel[topology_disc]))
        self.assertFalse(
            np.any(topology.branch_identity.stages.cleaned_skeleton[topology_disc])
        )
        self.assertFalse(np.any(topology.labels[topology_disc]))
        self.assertFalse(np.any(topology.annulus_masks[:, topology_disc]))
        self.assertTrue(np.any(topology.branch_identity.stages.vessel[disc]))
        self.assertTrue(np.any(topology.annulus_masks[:, disc]))
        self.assertGreater(topology.labels[43, 73], 0)
        self.assertTrue(np.any(topology.annulus_masks[:, 43, 105]))

    def test_wrong_disc_mask_shape_is_ignored_by_topology_but_rejected_on_output(self):
        vessel, _, source, settings = self._inputs()
        bad_disc = OpticDisc(
            np.zeros(vessel.shape[::-1], dtype=bool),
            source.source.segmentation.optic_disc.center,
            20.0,
            50.0,
        )
        branches = label_vessel_branches(vessel, bad_disc, settings)
        self.assertGreater(branches.branch_ids.size, 0)
        source.source.segmentation.optic_disc = bad_disc
        with self.assertRaisesRegex(ValueError, "same orientation and shape"):
            self._pack_source(source, branches, branches)

    @staticmethod
    def _pack_source(source_data, artery_segments, vein_segments):
        source = source_data.source
        return pack_topology_outputs(
            source.segmentation.vessels.artery,
            source.segmentation.vessels.vein,
            source.segmentation.optic_disc,
            {"artery": artery_segments, "vein": vein_segments},
            source.holodoppler.pixel_pitch,
        )

    @staticmethod
    def _inputs():
        shape = (101, 151)
        vessel = np.zeros(shape, dtype=bool)
        vessel[41:48, 5:146] = True
        vessel[5:96, 57:64] = True
        y, x = np.indices(shape)
        disc = ((x - 60) / 10.0) ** 2 + ((y - 45) / 25.0) ** 2 <= 1.0
        source = SimpleNamespace(
            source=SimpleNamespace(
                segmentation=SimpleNamespace(
                    vessels=SimpleNamespace(
                        artery=vessel,
                        vein=vessel.copy(),
                    ),
                    optic_disc=OpticDisc(disc, (60.0, 45.0), 20.0, 50.0),
                ),
                holodoppler=SimpleNamespace(
                    pixel_pitch=PixelPitch(20e-6, 20e-6),
                ),
            ),
            profile_settings="settings",
            provenance={"beat_index_base": 0},
        )
        settings = source.source.segmentation.optic_disc.annulus_geometry(shape)
        return vessel, disc, source, settings


if __name__ == "__main__":
    unittest.main()
