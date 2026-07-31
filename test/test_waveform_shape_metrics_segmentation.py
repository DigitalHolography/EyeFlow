"""Tests for canonical segmentation outputs published by waveform metrics."""

import unittest
from types import SimpleNamespace

import numpy as np

from input_output.schema import EyeFlowOutputPaths
from pipelines.waveform_shape_metrics.velocity.segmentation import (
    OPTIC_DISC_LABEL,
    REGION_AXIS_LABEL,
    pack_segmentation_outputs,
)


class SegmentationOutputTests(unittest.TestCase):
    def test_outputs_use_lower_left_xy_frame_and_enriched_branch_labels(self):
        image_shape = (16, 16)
        optic_disc_mask = np.zeros(image_shape, dtype=bool)
        optic_disc_mask[4, 4] = True
        artery_mask = np.zeros(image_shape, dtype=bool)
        artery_mask[1, 1] = True
        vein_mask = np.zeros(image_shape, dtype=bool)
        vein_mask[14, 14] = True

        labels = np.zeros(image_shape, dtype=np.int32)
        labels[2, 2] = 7
        labels[12, 12] = 9
        segments = SimpleNamespace(labels=labels)
        source_data = SimpleNamespace(
            retinal_artery_mask=artery_mask,
            retinal_vein_mask=vein_mask,
            optic_disc_mask=optic_disc_mask,
            optic_disc_center=np.asarray([8.0, 7.0], dtype=np.float32),
            optic_disc_width=np.float32(4.0),
            optic_disc_height=np.float32(4.0),
        )

        outputs = pack_segmentation_outputs(source_data, segments, segments)
        schema = EyeFlowOutputPaths.active()

        self.assertEqual(
            {
                schema.segmentation.optic_disc.mask,
                schema.segmentation.artery.mask,
                schema.segmentation.artery.branch_label_map,
                schema.segmentation.vein.mask,
                schema.segmentation.vein.branch_label_map,
            },
            set(outputs),
        )

        optic_disc, optic_disc_attrs = outputs[schema.segmentation.optic_disc.mask]
        np.testing.assert_array_equal(optic_disc[4, 11], True)
        self._assert_xy_lower_left_attrs(optic_disc_attrs)

        artery_mask_output, _ = outputs[schema.segmentation.artery.mask]
        np.testing.assert_array_equal(artery_mask_output[1, 14], True)

        branch_map, branch_attrs = outputs[
            schema.segmentation.artery.branch_label_map
        ]
        self.assertEqual((16, 16), branch_map.shape)
        self.assertEqual(7, branch_map[2, 13])
        self.assertEqual(9, branch_map[12, 3])
        self.assertEqual(OPTIC_DISC_LABEL, branch_map[4, 11])
        self.assertEqual(REGION_AXIS_LABEL, branch_map[8, 0])
        self.assertEqual(REGION_AXIS_LABEL, branch_map[0, 8])
        self._assert_xy_lower_left_attrs(branch_attrs)
        self.assertEqual(OPTIC_DISC_LABEL, branch_attrs["optic_disc_label"])
        self.assertEqual(REGION_AXIS_LABEL, branch_attrs["axis_label"])

    @staticmethod
    def _assert_xy_lower_left_attrs(attrs):
        assert attrs["dimDesc"] == ["x", "y"]
        assert attrs["image_origin"] == "lower_left"
        assert attrs["y_axis_direction"] == "increasing_toward_north"


if __name__ == "__main__":
    unittest.main()
