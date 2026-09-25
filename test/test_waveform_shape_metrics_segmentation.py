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
    _topology_delta_radius,
    pack_segmentation_outputs,
    pack_topology_outputs,
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
        pixel_pitch_m, pixel_pitch_attrs = outputs[
            schema.segmentation.pixel_pitch_m
        ]
        self.assertEqual(np.float32, height.dtype)
        self.assertEqual(np.float32, width.dtype)
        self.assertEqual(np.float32, center.dtype)
        self.assertEqual(np.float32, pixel_pitch_m.dtype)
        self.assertEqual(np.float32(4.0), height)
        self.assertEqual(np.float32(4.0), width)
        np.testing.assert_array_equal(center, [8.0, 8.0])
        self.assertAlmostEqual(float(pixel_pitch_m), 1.91e-3 / 4.0)
        self.assertEqual("pixels", height_attrs["unit"])
        self.assertEqual("pixels", width_attrs["unit"])
        self.assertEqual(["coordinate"], center_attrs["dimDesc"])
        self.assertEqual(["x", "y"], center_attrs["coordinates"])
        self.assertEqual("lower_left", center_attrs["image_origin"])
        self.assertEqual("m", pixel_pitch_attrs["unit"])

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

    def test_lumen_diameter_uses_each_vessels_area_and_float32_delta_radius(self):
        image_shape = (4, 4)
        annuli = np.zeros((4, *image_shape), dtype=bool)
        annuli[0, :2] = True
        annuli[1, 2:] = True
        annuli[2, 2:] = True
        artery_labels = np.zeros(image_shape, dtype=np.int32)
        artery_labels[0, :2] = 1
        artery_labels[2, :3] = 1
        vein_labels = np.zeros(image_shape, dtype=np.int32)
        vein_labels[0, :4] = 1
        vein_labels[2, :1] = 1
        settings = AnnulusGeometry(0.0, 1.0, 0.25, 4)

        def topology(labels):
            return SimpleNamespace(
                labels=labels,
                branch_ids=np.asarray([1, 2], dtype=np.int32),
                annulus_masks=annuli,
                ring_settings=settings,
                delta_radius=np.asarray(
                    [2.0, 0.0, np.nan, 2.0],
                    dtype=np.float32,
                ),
            )

        outputs = pack_topology_outputs(
            artery_labels > 0,
            vein_labels > 0,
            OpticDisc(None, (1.5, 1.5), 1.0, 1.0),
            {
                "artery": topology(artery_labels),
                "vein": topology(vein_labels),
            },
        )
        schema = EyeFlowOutputPaths.active()

        artery_diameter, artery_attrs = outputs[
            schema.segmentation.artery.lumen_diameter
        ]
        vein_diameter, _ = outputs[schema.segmentation.vein.lumen_diameter]
        delta_radius, delta_attrs = outputs[
            schema.segmentation.artery.delta_radius
        ]

        self.assertEqual(np.float32, artery_diameter.dtype)
        self.assertEqual(np.float32, vein_diameter.dtype)
        self.assertEqual(np.float32, delta_radius.dtype)
        np.testing.assert_allclose(
            delta_radius,
            [2.0, np.nan, np.nan, np.nan],
            equal_nan=True,
        )
        np.testing.assert_allclose(
            artery_diameter[:, 0],
            [1.0, np.nan],
            equal_nan=True,
        )
        np.testing.assert_allclose(
            vein_diameter[:, 0],
            [2.0, np.nan],
            equal_nan=True,
        )
        self.assertTrue(np.all(np.isnan(artery_diameter[:, 1:])))
        self.assertTrue(np.all(np.isnan(vein_diameter[:, 1:])))
        self.assertEqual("pixels", artery_attrs["unit"])
        self.assertEqual(["branch", "radius"], artery_attrs["dimDesc"])
        self.assertEqual("pixels", delta_attrs["unit"])
        self.assertEqual(["radius"], delta_attrs["dimDesc"])

    def test_missing_delta_radius_is_float32_nan(self):
        delta_radius = _topology_delta_radius(SimpleNamespace(), 2)

        self.assertEqual(np.float32, delta_radius.dtype)
        self.assertTrue(np.all(np.isnan(delta_radius)))

    @staticmethod
    def _assert_xy_lower_left_attrs(attrs):
        assert attrs["dimDesc"] == ["x", "y"]
        assert attrs["image_origin"] == "lower_left"
        assert attrs["y_axis_direction"] == "increasing_toward_north"


if __name__ == "__main__":
    unittest.main()
