"""Tests for spatial segment geometry used to build segment waveforms."""

from __future__ import annotations

import sys
import unittest
from types import SimpleNamespace
from unittest.mock import patch
from pathlib import Path

import numpy as np

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from calculations.blood_flow_velocity import CrossSectionSignalSettings  # noqa: E402
from calculations.blood_flow_velocity.cross_section.generate_cross_section_signals import (  # noqa: E402
    _CrossSectionBuffers,
    _cross_section_velocity,
    _fill_cross_section_buffers,
    _subimage_stack,
    _subimage_stack2,
)
from calculations.blood_flow_velocity.cross_section.branch_identity import (  # noqa: E402
    _branch_identity_stages,
)
from calculations.blood_flow_velocity.cross_section.segment_geometry import (  # noqa: E402
    annulus_mask,
    ring_masks,
    section_masks,
)
from pipelines.waveform_velocity_core.runner import _segment_ring_settings  # noqa: E402
from pipelines.waveform_velocity_core.branch_identity_debug import (  # noqa: E402
    _labels_with_substack_boxes,
)
from utils.logger import Logger  # noqa: E402


class SegmentCenterTests(unittest.TestCase):
    def tearDown(self) -> None:
        Logger.reset_current()

    def test_watershed_boundaries_remain_split_during_annulus_relabeling(self) -> None:
        vessel = np.ones((9, 9), dtype=bool)
        section = np.ones_like(vessel)
        skeleton = np.zeros_like(vessel)
        skeleton[4, 2] = True
        skeleton[4, 6] = True
        watershed_labels = np.zeros(vessel.shape, dtype=np.int32)
        watershed_labels[:, 1:4] = 1
        watershed_labels[:, 4:8] = 2
        settings = SimpleNamespace(
            inner_radius_frac=0.0,
            outer_radius_frac=0.5,
            ring_width_frac=0.5,
            ring_count=1,
        )

        with patch(
            "calculations.blood_flow_velocity.cross_section.branch_identity."
            "segment_annulus_mask",
            return_value=section,
        ), patch(
            "calculations.blood_flow_velocity.cross_section.branch_identity."
            "skeletonize",
            return_value=skeleton,
        ), patch(
            "calculations.blood_flow_velocity.cross_section.branch_identity."
            "_branch_points",
            return_value=np.zeros_like(vessel),
        ), patch(
            "calculations.blood_flow_velocity.cross_section.branch_identity."
            "_remove_small",
            side_effect=lambda mask, _min_area: mask,
        ), patch(
            "calculations.blood_flow_velocity.cross_section.branch_identity."
            "watershed",
            return_value=watershed_labels,
        ):
            stages = _branch_identity_stages(vessel, None, settings)

        self.assertEqual(2, int(stages.annulus_refined_labels.max()))

    def test_subimage_uses_matlab_fixed_scale_factor_width(self) -> None:
        velocity = np.ones((2, 101, 101), dtype=np.float32)
        segment = np.zeros((101, 101), dtype=bool)
        segment[49:52, 60:71] = True
        messages: list[str] = []
        settings = CrossSectionSignalSettings(3.0, False, 0.5, False, 0.01)

        Logger.configure(on_log=messages.append)
        stack, mask = _subimage_stack(velocity, segment, (65, 50), settings)

        self.assertEqual((2, 5, 5), stack.shape)
        self.assertEqual((5, 5), mask.shape)
        self.assertEqual(1, len(messages))
        self.assertIn("5x5 px window for a 11x3 px segment", messages[0])

    def test_subimage_stack2_contains_segment_with_rotation_margin(self) -> None:
        velocity = np.ones((2, 101, 101), dtype=np.float32)
        segment = np.zeros((101, 101), dtype=bool)
        segment[49:52, 60:71] = True
        messages: list[str] = []
        settings = CrossSectionSignalSettings(3.0, False, 0.5, False, 0.01)

        Logger.configure(on_log=messages.append)
        stack, mask = _subimage_stack2(velocity, segment, (65, 50), settings)

        self.assertEqual((2, 15, 15), stack.shape)
        self.assertEqual((15, 15), mask.shape)
        self.assertEqual(33, int(np.count_nonzero(mask)))
        self.assertEqual([], messages)

    def test_substack_debug_overlay_marks_shared_box_edges(self) -> None:
        labels = np.ones((20, 20), dtype=np.int32)
        ring_settings = SimpleNamespace(
            inner_radius_frac=0.0,
            outer_radius_frac=1.0,
            ring_width_frac=0.5,
            ring_count=2,
        )
        cross_section_settings = CrossSectionSignalSettings(
            3.0,
            False,
            0.5,
            False,
            0.01,
        )
        centers = np.asarray(
            [[[5.0, 5.0], [np.nan, np.nan]], [[7.0, 5.0], [np.nan, np.nan]]],
            dtype=np.float32,
        )

        image = _labels_with_substack_boxes(
            labels,
            None,
            ring_settings,
            centers,
            cross_section_settings,
        )

        self.assertEqual((20, 20, 3), image.shape)
        np.testing.assert_array_equal(image[4, 4], [255, 255, 0])
        np.testing.assert_array_equal(image[4, 6], [255, 0, 255])

    def test_subimage_logs_error_only_when_segment_reaches_crop_edge(self) -> None:
        velocity = np.zeros((1, 100, 100), dtype=np.float32)
        settings = CrossSectionSignalSettings(3.0, False, 0.5, False, 0.01)
        messages: list[str] = []
        Logger.configure(on_log=messages.append)

        contained = np.zeros((100, 100), dtype=bool)
        contained[50, 50] = True
        _subimage_stack(velocity, contained, (50, 50), settings)
        self.assertEqual([], messages)

        clipped = contained.copy()
        clipped[50, 48] = True
        _subimage_stack(velocity, clipped, (50, 50), settings)
        self.assertEqual(1, len(messages))
        self.assertTrue(messages[0].startswith("[WARNING] Cross-section window"))
        self.assertIn("3x1 px segment", messages[0])

    def test_missing_annulus_intersections_do_not_select_another_crop_path(self) -> None:
        velocity = np.ones((2, 101, 101), dtype=np.float32)
        segment = np.zeros((101, 101), dtype=bool)
        segment[49:52, 60:71] = True
        settings = CrossSectionSignalSettings(3.0, False, 0.5, False, 0.01)

        with patch(
            "calculations.blood_flow_velocity.cross_section."
            "generate_cross_section_signals._circle_tilt_geometry",
            return_value=None,
        ), patch(
            "calculations.blood_flow_velocity.cross_section."
            "generate_cross_section_signals._estimate_orientation",
            return_value=0.0,
        ):
            measurement = _cross_section_velocity(
                velocity,
                segment,
                segment,
                (65, 50),
                (50, 50),
                settings,
            )

        self.assertEqual((2, 5), measurement[2].shape)

    def test_segment_analysis_uses_disc_extent_and_corner_spacing(self) -> None:
        settings = _segment_ring_settings(
            55,
            69,
            image_shape=(512, 512),
            optic_disc_center=np.asarray([267.0, 230.0]),
        )

        corner_radius = np.hypot(267.0, 281.0)
        self.assertEqual(37, settings.ring_count)
        self.assertAlmostEqual(
            (69.0 / 2.0) / corner_radius,
            settings.inner_radius_frac,
        )
        self.assertEqual(1.0, settings.outer_radius_frac)
        self.assertEqual(0.025, settings.ring_width_frac)
        self.assertEqual(0.025, settings.segment_length_frac)

    def test_segment_settings_carry_detected_optic_disc_dimensions(self) -> None:
        settings = _segment_ring_settings(55, 69)

        self.assertEqual(55.0, settings.optic_disc_width)
        self.assertEqual(69.0, settings.optic_disc_height)

    def test_annulus_starts_at_optic_disc_and_uses_corner_extent(self) -> None:
        center = np.asarray([267.0, 230.0])
        settings = _segment_ring_settings(
            55,
            69,
            image_shape=(512, 512),
            optic_disc_center=center,
        )
        mask = section_masks(
            (512, 512),
            center,
            settings,
        )[0]

        # The first annulus starts at max(width, height) / 2 = 34.5 px.
        self.assertFalse(mask[230, 301])
        self.assertTrue(mask[230, 302])
        self.assertFalse(mask[264, 267])
        self.assertTrue(mask[265, 267])

        # The first annulus has the configured 0.025 radial spacing.
        self.assertTrue(mask[230, 311])
        self.assertFalse(mask[230, 312])
        self.assertTrue(mask[274, 267])
        self.assertFalse(mask[275, 267])

        rings = ring_masks((512, 512), center, settings)
        self.assertTrue(rings[-1, 511, 0])
        self.assertFalse(rings[-1, 0, 0])

    def test_annulus_is_independent_of_invalid_disc_dimensions(self) -> None:
        legacy = annulus_mask((100, 200), None, 0.10, 0.35)
        fallback = annulus_mask(
            (100, 200),
            None,
            0.10,
            0.35,
            optic_disc_width=0,
            optic_disc_height=np.nan,
            optic_disc_boundary_radius_frac=0.10,
        )

        np.testing.assert_array_equal(legacy, fallback)

    def test_cross_section_centroids_are_recorded_in_waveform_axis_order(self) -> None:
        labels = np.zeros((7, 7), dtype=np.int32)
        labels[2, 1] = 1
        labels[5, 4] = 2
        branches = SimpleNamespace(
            labels=labels,
            branch_ids=np.asarray([1, 2], dtype=np.int32),
        )
        sections = np.zeros((2, 7, 7), dtype=bool)
        sections[0, 2, 1] = True
        sections[1, 5, 4] = True
        velocity = np.zeros((3, 7, 7), dtype=np.float32)
        buffers = _CrossSectionBuffers.allocate(
            frame_count=velocity.shape[0],
            ring_count=sections.shape[0],
            branch_count=branches.branch_ids.size,
        )
        settings = CrossSectionSignalSettings(3.0, False, 0.5, False, 0.01)

        signal = np.arange(3, dtype=np.float32)
        with patch(
            "calculations.blood_flow_velocity.cross_section."
            "generate_cross_section_signals._tilt_angle",
            return_value=0.0,
        ), patch(
            "calculations.blood_flow_velocity.cross_section."
            "generate_cross_section_signals._cross_section_velocity",
            return_value=(signal, signal),
        ):
            _fill_cross_section_buffers(
                buffers,
                velocity,
                sections,
                branches,
                np.asarray([3.0, 3.0], dtype=np.float32),
                settings,
            )

        np.testing.assert_array_equal(
            buffers.segment_center_xy[0, 0],
            [1.0, 2.0],
        )
        np.testing.assert_array_equal(
            buffers.segment_center_xy[1, 1],
            [4.0, 5.0],
        )
        self.assertTrue(np.all(np.isnan(buffers.segment_center_xy[0, 1])))
        self.assertTrue(np.all(np.isnan(buffers.segment_center_xy[1, 0])))
        np.testing.assert_array_equal(buffers.velocity[0, 0], signal)
        np.testing.assert_array_equal(buffers.velocity[1, 1], signal)


if __name__ == "__main__":
    unittest.main()
