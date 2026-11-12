"""Tests for canonical segmentation outputs published by waveform metrics."""

import unittest
from types import SimpleNamespace

import numpy as np

from calculations.topology import AnnulusGeometry, OpticDisc
from input_output.schema import EyeFlowOutputPaths, PixelPitch
from pipelines.topology_core.outputs import (
    ANNULUS_OUTLINE_LABEL,
    BACKGROUND_LABEL,
    INNER_R0_VESSEL_LABEL,
    pack_topology_outputs,
)


class SegmentationOutputTests(unittest.TestCase):
    def test_missing_dopplerview_optic_disc_publishes_nan_vein_maps(self):
        image_shape = (8, 10)
        artery = np.zeros(image_shape, dtype=bool)
        artery[2, 3] = True
        vein = np.zeros(image_shape, dtype=bool)
        optic_disc = OpticDisc(
            np.zeros(image_shape, dtype=bool),
            (5.0, 4.0),
            2.0,
            2.0,
            is_fallback=True,
        )
        source_data = SimpleNamespace(
            source=SimpleNamespace(
                segmentation=SimpleNamespace(
                    vessels=SimpleNamespace(artery=artery, vein=vein),
                    optic_disc=optic_disc,
                ),
                holodoppler=SimpleNamespace(
                    pixel_pitch=PixelPitch(20e-6, 20e-6),
                ),
            ),
        )

        outputs = pack_topology_outputs(
            artery, vein, optic_disc, {}, source_data.source.holodoppler.pixel_pitch,
        )
        schema = EyeFlowOutputPaths.active()

        artery_output, _ = outputs[schema.segmentation.artery.mask]
        vein_output, vein_attrs = outputs[schema.segmentation.vein.mask]
        vein_branches, _ = outputs[schema.segmentation.vein.branch_label_map]
        vein_segments, _ = outputs[schema.segmentation.vein.segment_map]
        self.assertTrue(artery_output[3, 5])
        self.assertTrue(np.all(np.isnan(vein_output)))
        self.assertTrue(np.all(np.isnan(vein_branches)))
        self.assertTrue(np.all(np.isnan(vein_segments)))
        self.assertEqual(
            "unavailable_without_dopplerview_optic_disc",
            vein_attrs["source"],
        )

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
            source=SimpleNamespace(
                segmentation=SimpleNamespace(
                    vessels=SimpleNamespace(
                        artery=artery_mask,
                        vein=vein_mask,
                    ),
                    optic_disc=OpticDisc(
                        optic_disc_mask,
                        (8.0, 7.0),
                        4.0,
                        4.0,
                    ),
                ),
                holodoppler=SimpleNamespace(
                    pixel_pitch=PixelPitch(20e-6, 20e-6),
                ),
            ),
        )

        outputs = pack_topology_outputs(
            artery_mask, vein_mask, source_data.source.segmentation.optic_disc,
            {"artery": segments, "vein": segments},
            source_data.source.holodoppler.pixel_pitch,
        )
        schema = EyeFlowOutputPaths.active()

        self.assertEqual(
            {
                schema.segmentation.optic_disc.mask,
                schema.segmentation.optic_disc.height,
                schema.segmentation.optic_disc.width,
                schema.segmentation.optic_disc.center,
                schema.segmentation.pixel_pitch_m,
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

        height, height_attrs = outputs[schema.segmentation.optic_disc.height]
        width, width_attrs = outputs[schema.segmentation.optic_disc.width]
        center, center_attrs = outputs[schema.segmentation.optic_disc.center]
        self.assertEqual("Meta/PixelPitch_m", schema.segmentation.pixel_pitch_m)
        self.assertNotIn("Segmentation/PixelPitch_m/value", outputs)
        pixel_pitch_m, pixel_pitch_attrs = outputs[schema.segmentation.pixel_pitch_m]
        self.assertEqual(np.float32, height.dtype)
        self.assertEqual(np.float32, width.dtype)
        self.assertEqual(np.float32, center.dtype)
        self.assertEqual(np.float32, pixel_pitch_m.dtype)
        self.assertEqual(np.float32(4.0), height)
        self.assertEqual(np.float32(4.0), width)
        np.testing.assert_array_equal(center, [8.0, 8.0])
        self.assertAlmostEqual(float(pixel_pitch_m), 20e-6)
        self.assertEqual("pixels", height_attrs["unit"])
        self.assertEqual("pixels", width_attrs["unit"])
        self.assertEqual(["coordinate"], center_attrs["dimDesc"])
        self.assertEqual(["x", "y"], center_attrs["coordinates"])
        self.assertEqual("lower_left", center_attrs["image_origin"])
        self.assertEqual("m", pixel_pitch_attrs["unit"])

        artery_mask_output, _ = outputs[schema.segmentation.artery.mask]
        np.testing.assert_array_equal(artery_mask_output[2, 13], True)

        branch_map, branch_attrs = outputs[schema.segmentation.artery.branch_label_map]
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

    @staticmethod
    def _assert_xy_lower_left_attrs(attrs):
        assert attrs["dimDesc"] == ["x", "y"]
        assert attrs["image_origin"] == "lower_left"
        assert attrs["y_axis_direction"] == "increasing_toward_north"


if __name__ == "__main__":
    unittest.main()
