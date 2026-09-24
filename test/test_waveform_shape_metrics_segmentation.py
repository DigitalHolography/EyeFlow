"""Tests for canonical segmentation outputs published by waveform metrics."""

import unittest
from types import SimpleNamespace

import numpy as np

from calculations.topology import AnnulusGeometry, OpticDisc

from input_output.schema import EyeFlowOutputPaths
from pipelines.waveform_velocity_core.segmentation import (
    ANNULUS_OUTLINE_LABEL,
    BACKGROUND_LABEL,
    INNER_R0_VESSEL_LABEL,
    _annulus_outlines,
    pack_segmentation_outputs,
)


class SegmentationOutputTests(unittest.TestCase):
    def test_outputs_use_lower_left_xy_frame_and_distinct_label_classes(self):
        image_shape = (16, 16)
        optic_disc_mask = np.zeros(image_shape, dtype=bool)
        optic_disc_mask[4, 4] = True
        artery_mask = np.zeros(image_shape, dtype=bool)
        artery_mask[2, 2] = True
        artery_mask[7, 8] = True
        artery_mask[12, 12] = True
        vein_mask = np.zeros(image_shape, dtype=bool)
        vein_mask[:] = artery_mask

        labels = np.zeros(image_shape, dtype=np.int32)
        labels[2, 2] = 7
        labels[12, 12] = 9
        segments = SimpleNamespace(labels=labels)
        source_data = SimpleNamespace(
            retinal_artery_mask=artery_mask,
            retinal_vein_mask=vein_mask,
            optic_disc=OpticDisc(optic_disc_mask, (8.0, 7.0), 4.0, 4.0),
        )

        outputs = pack_segmentation_outputs(source_data, segments, segments)
        schema = EyeFlowOutputPaths.active()

        self.assertEqual(
            {
                schema.segmentation.optic_disc.mask,
                schema.segmentation.artery.mask,
                schema.segmentation.artery.branch_label_map,
                schema.segmentation.artery.segment_map,
                schema.segmentation.vein.mask,
                schema.segmentation.vein.branch_label_map,
                schema.segmentation.vein.segment_map,
            },
            set(outputs),
        )

        optic_disc, optic_disc_attrs = outputs[schema.segmentation.optic_disc.mask]
        np.testing.assert_array_equal(optic_disc[4, 11], True)
        self._assert_xy_lower_left_attrs(optic_disc_attrs)

        artery_mask_output, _ = outputs[schema.segmentation.artery.mask]
        np.testing.assert_array_equal(artery_mask_output[2, 13], True)

        branch_map, branch_attrs = outputs[
            schema.segmentation.artery.branch_label_map
        ]
        self.assertEqual((16, 16), branch_map.shape)
        self.assertEqual(0, branch_map[2, 13])
        self.assertEqual(1, branch_map[12, 3])
        self.assertEqual(INNER_R0_VESSEL_LABEL, branch_map[8, 8])
        self.assertEqual(ANNULUS_OUTLINE_LABEL, branch_map[10, 8])
        self.assertEqual(BACKGROUND_LABEL, branch_map[0, 0])
        np.testing.assert_array_equal(
            np.flatnonzero(branch_map[:, 8] == ANNULUS_OUTLINE_LABEL),
            [6, 10],
        )
        self._assert_xy_lower_left_attrs(branch_attrs)
        self.assertEqual(BACKGROUND_LABEL, branch_attrs["background_label"])
        self.assertEqual(
            ANNULUS_OUTLINE_LABEL,
            branch_attrs["annulus_outline_label"],
        )
        self.assertEqual(
            INNER_R0_VESSEL_LABEL,
            branch_attrs["inner_r0_vessel_label"],
        )

        segment_map, segment_attrs = outputs[schema.segmentation.artery.segment_map]
        self.assertGreater(
            np.count_nonzero(segment_map == ANNULUS_OUTLINE_LABEL),
            np.count_nonzero(branch_map == ANNULUS_OUTLINE_LABEL),
        )
        self.assertEqual(
            "all calculated annulus outlines",
            segment_attrs["annulus_outlines"],
        )

    def test_final_annulus_outline_keeps_the_regular_spacing(self):
        shape = (101, 101)
        settings = AnnulusGeometry(
            inner_radius_frac=0.1,
            outer_radius_frac=0.5,
            ring_width_frac=0.2,
            ring_count=3,
        )

        outlines = _annulus_outlines(
            shape,
            np.asarray([50.0, 50.0]),
            7,
            settings,
        )

        # The last outline is at 0.7 of the half-diagonal, not clipped to 0.5.
        self.assertTrue(outlines[50, 99])

    @staticmethod
    def _assert_xy_lower_left_attrs(attrs):
        assert attrs["dimDesc"] == ["x", "y"]
        assert attrs["image_origin"] == "lower_left"
        assert attrs["y_axis_direction"] == "increasing_toward_north"


if __name__ == "__main__":
    unittest.main()
