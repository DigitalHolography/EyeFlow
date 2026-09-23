"""Regression tests for branch extraction at the actual optic-disc boundary."""

import sys
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

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
from input_output.schema import EyeFlowOutputPaths
from pipelines.waveform_velocity_core import runner
from pipelines.waveform_velocity_core.segmentation import (
    OPTIC_DISC_LABEL,
    REGION_AXIS_LABEL,
    pack_segmentation_outputs,
)


class OpticDiscBranchMaskTests(unittest.TestCase):
    def test_empty_mask_retains_the_configured_inner_boundary(self):
        shape = (21, 21)
        vessel = np.ones(shape, dtype=bool)
        center = (10.0, 10.0)
        settings = AnnulusGeometry(0.4, 0.9, 0.5, 1)
        empty = label_vessel_branches(
            vessel,
            OpticDisc(np.zeros(shape, dtype=bool), center, None, None),
            settings,
        )
        exact = np.zeros(shape, dtype=bool)
        exact[10, 10] = True
        nonempty = label_vessel_branches(
            vessel,
            OpticDisc(exact, center, None, None),
            settings,
        )

        self.assertFalse(empty.stages.section[10, 13])
        self.assertTrue(empty.stages.section[10, 18])
        self.assertTrue(nonempty.stages.section[10, 13])

    def test_elliptical_disc_keeps_vessels_inside_the_old_circular_cutoff(self):
        vessel, disc, source, settings = self._inputs()
        original_vessel = vessel.copy()
        original_disc = disc.copy()
        circular_disc = OpticDisc(
            np.zeros_like(disc), source.optic_disc.center, 20.0, 50.0
        )
        circular = label_vessel_branches(vessel, circular_disc, settings)
        branches = label_vessel_branches(
            vessel, source.optic_disc, settings
        )

        # This vessel lies outside the ellipse but inside its bounding circle.
        self.assertFalse(disc[43, 73])
        self.assertEqual(0, circular.labels[43, 73])
        self.assertGreater(branches.labels[43, 73], 0)
        self.assertEqual(4, branches.branch_ids.size)
        self.assertFalse(np.any(branches.labels[disc]))
        self.assertFalse(np.any(branches.stages.skeleton[disc]))
        np.testing.assert_array_equal(vessel, original_vessel)
        np.testing.assert_array_equal(disc, original_disc)

        outputs = pack_segmentation_outputs(source, branches, branches)
        schema = EyeFlowOutputPaths.active()
        mask, _ = outputs[schema.segmentation.optic_disc.mask]
        for paths in (schema.segmentation.artery, schema.segmentation.vein):
            image, _ = outputs[paths.branch_label_map]
            self.assertGreater(image[73, 57], 0)
            self.assertFalse(np.any(image[mask] > 0))
            np.testing.assert_array_equal(
                mask[image != REGION_AXIS_LABEL],
                (image == OPTIC_DISC_LABEL)[image != REGION_AXIS_LABEL],
            )

    def test_actual_mask_also_excludes_disc_pixels_from_measurement_rings(self):
        vessel, disc, source, settings = self._inputs()
        # A nonelliptical extension reaches beyond the nominal disc dimensions.
        disc[40:49, 90:97] = True
        topology = build_segment_topology(
            vessel,
            OpticDisc(disc, source.optic_disc.center, 20.0, 50.0),
            settings,
        )

        self.assertFalse(np.any(topology.branch_identity.stages.vessel[disc]))
        self.assertFalse(np.any(topology.labels[disc]))
        self.assertFalse(np.any(topology.annulus_masks[:, disc]))
        self.assertGreater(topology.labels[43, 73], 0)
        self.assertTrue(np.any(topology.annulus_masks[:, 43, 105]))

    def test_pipeline_uses_the_published_disc_mask_with_geometry_fallback(self):
        vessel, disc, source, settings = self._inputs()
        ctx = SimpleNamespace(
            output=SimpleNamespace(available=False),
            pipeline_scheduled=lambda name: False,
            option_enabled=lambda *args, **kwargs: False,
            inputs=SimpleNamespace(
                hd=SimpleNamespace(filename="hd.h5"),
                dv=SimpleNamespace(filename="dv.h5"),
            ),
            state=SimpleNamespace(raw={}),
        )
        schema = EyeFlowOutputPaths.active()
        for source_mask in (disc, None):
            with self.subTest(has_source_mask=source_mask is not None):
                source.optic_disc = OpticDisc(
                    source_mask,
                    (60.0, 45.0),
                    20.0,
                    50.0,
                )
                published, _ = pack_segmentation_outputs(source, None, None)[
                    schema.segmentation.optic_disc.mask
                ]
                expected_disc = np.flip(published.T, axis=0)
                with patch.object(
                    runner,
                    "analyze_velocity_segment_profiles",
                    return_value={"artery": "artery", "vein": "vein"},
                ) as extract:
                    result = runner._segment_velocity_inputs(
                        np.zeros((1, *vessel.shape), dtype=np.float32),
                        source,
                        settings,
                        ctx,
                        cycle_boundary_indexes=np.asarray([0, 0]),
                    )
                self.assertEqual(("artery", "vein"), result)
                np.testing.assert_array_equal(
                    extract.call_args.args[2].mask_for(vessel.shape), expected_disc
                )
                branches = label_vessel_branches(
                    vessel, source.optic_disc, settings,
                )
                self.assertGreater(branches.labels[43, 73], 0)

    def test_wrong_disc_mask_shape_is_rejected(self):
        vessel, _, source, settings = self._inputs()
        with self.assertRaisesRegex(ValueError, "same orientation and shape"):
            label_vessel_branches(
                vessel,
                OpticDisc(
                    np.zeros(vessel.shape[::-1], dtype=bool),
                    source.optic_disc.center,
                    20.0,
                    50.0,
                ),
                settings,
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
            retinal_artery_mask=vessel,
            retinal_vein_mask=vessel.copy(),
            optic_disc=OpticDisc(disc, (60.0, 45.0), 20.0, 50.0),
            cross_section_settings="settings",
            provenance={"beat_index_base": 0},
        )
        settings = source.optic_disc.annulus_geometry(shape)
        return vessel, disc, source, settings


if __name__ == "__main__":
    unittest.main()
