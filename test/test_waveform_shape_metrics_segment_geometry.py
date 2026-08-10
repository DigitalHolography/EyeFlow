"""Tests for spatial segment geometry used to build segment waveforms."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from calculations.blood_flow_velocity import CrossSectionSignalSettings  # noqa: E402
from calculations.blood_flow_velocity.cross_section.branch_identity import (  # noqa: E402
    _branch_identity_stages,
)
from calculations.blood_flow_velocity.cross_section.generate_cross_section_signals import (  # noqa: E402
    _cross_section_velocity,
    _CrossSectionBuffers,
    _fill_cross_section_buffers,
    _frame_velocities,
    _subimage_stack,
    _subimage_stack2,
)
from calculations.blood_flow_velocity.cross_section.reusable_cross_section_signals import (  # noqa: E402
    generate_cross_section_signals_for_cubes,
    project_cross_section_cube,
)
from calculations.blood_flow_velocity.cross_section.segment_geometry import (  # noqa: E402
    annulus_mask,
    ring_masks,
    section_masks,
)
from pipelines.waveform_velocity_core.branch_identity_debug import (  # noqa: E402
    _labels_with_substack_boxes,
)
from pipelines.waveform_velocity_core.runner import _segment_ring_settings  # noqa: E402
from utils.logger import Logger  # noqa: E402


class SegmentCenterTests(unittest.TestCase):
    def tearDown(self) -> None:
        Logger.reset_current()

    def test_raw_segment_velocity_is_the_unrotated_masked_mean(self) -> None:
        sub_stack = np.asarray(
            [
                [[1.0, np.nan, 9.0], [3.0, np.nan, np.nan]],
                [[2.0, 4.0, np.nan], [6.0, np.nan, np.nan]],
            ],
            dtype=np.float32,
        )

        with patch(
            "calculations.blood_flow_velocity.cross_section."
            "generate_cross_section_signals.rotate_image_with_nan",
            side_effect=AssertionError("raw segment velocity must not rotate"),
        ) as rotate:
            raw, safe, profiles = _frame_velocities(
                sub_stack,
                angle=90.0,
                c1=0,
                c2=0,
            )

        rotate.assert_not_called()
        np.testing.assert_allclose(raw, [13.0 / 3.0, 4.0])
        np.testing.assert_allclose(
            profiles,
            [[2.0, np.nan, 9.0], [4.0, 4.0, np.nan]],
            equal_nan=True,
        )
        np.testing.assert_allclose(safe, [5.5, 4.0])

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

        self.assertEqual((2, 15), measurement[2].shape)

    def test_segment_analysis_uses_disc_extent_and_corner_spacing(self) -> None:
        settings = _segment_ring_settings(
            55,
            69,
            image_shape=(512, 512),
            optic_disc_center=np.asarray([267.0, 230.0]),
        )

        corner_radius = np.hypot(267.0, 281.0)
        self.assertEqual(23, settings.ring_count)
        self.assertAlmostEqual(
            (69.0 / 2.0) / corner_radius,
            settings.inner_radius_frac,
        )
        self.assertEqual(1.0, settings.outer_radius_frac)
        self.assertEqual(0.04, settings.ring_width_frac)
        self.assertEqual(0.04, settings.segment_length_frac)

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

        # The first annulus has the configured 0.04 radial spacing.
        self.assertTrue(mask[230, 311])
        self.assertTrue(mask[230, 317])
        self.assertFalse(mask[230, 318])
        self.assertTrue(mask[274, 267])

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


class ReusableCrossSectionProjectionTests(unittest.TestCase):
    def test_named_cubes_reuse_adaptive_window_angle_and_limits(self) -> None:
        reference, vessel, branches, sections, ring_settings, settings = (
            self._projection_inputs()
        )
        additional = {
            "second": reference.copy(),
            "third": reference * np.float32(3.0),
            "fourth": reference + np.float32(7.0),
        }

        with (
            patch(
                "calculations.blood_flow_velocity.cross_section."
                "generate_cross_section_signals.label_vessel_branches",
                return_value=branches,
            ),
            patch(
                "calculations.blood_flow_velocity.cross_section."
                "generate_cross_section_signals.section_masks",
                return_value=sections,
            ),
            patch(
                "calculations.blood_flow_velocity.cross_section."
                "reusable_cross_section_signals.section_masks",
                return_value=sections,
            ),
            patch(
                "calculations.blood_flow_velocity.cross_section."
                "generate_cross_section_signals._estimate_orientation",
                return_value=37.0,
            ) as estimate_orientation,
            patch(
                "calculations.blood_flow_velocity.cross_section."
                "generate_cross_section_signals._cross_section_limits",
                return_value=(1, 5),
            ) as cross_section_limits,
        ):
            result = generate_cross_section_signals_for_cubes(
                reference,
                vessel,
                (10, 10),
                ring_settings,
                settings,
                additional_cubes=additional,
                reference_name="velocity",
            )

        self.assertEqual(
            ["velocity", "second", "third", "fourth"],
            list(result.passes),
        )
        self.assertEqual(1, estimate_orientation.call_count)
        self.assertEqual(1, cross_section_limits.call_count)
        np.testing.assert_array_equal(
            result.plan.rotation_angles_degrees,
            [[37.0]],
        )
        np.testing.assert_array_equal(
            result.plan.profile_window_bounds_xyxy,
            [[[6, 15, 6, 15]]],
        )
        np.testing.assert_array_equal(
            result.plan.profile_integration_limits_pixels,
            [[[1, 5]]],
        )
        for name in ("second", "third", "fourth"):
            projected = result.passes[name]
            np.testing.assert_array_equal(
                projected.profile_rotation_degrees,
                result.plan.profile_rotation_degrees,
            )
            np.testing.assert_array_equal(
                projected.profile_window_bounds_xyxy,
                result.plan.profile_window_bounds_xyxy,
            )
            np.testing.assert_array_equal(
                projected.profile_integration_limits_pixels,
                result.plan.profile_integration_limits_pixels,
            )
        np.testing.assert_array_equal(
            result.passes["second"].velocity,
            result.reference.velocity,
        )
        self.assertFalse(
            np.array_equal(
                result.passes["third"].velocity,
                result.reference.velocity,
            )
        )
        self.assertIs(
            result.passes["second"].projected_signal,
            result.passes["second"].velocity,
        )

    def test_per_cube_mode_refits_only_profile_limits(self) -> None:
        reference, vessel, branches, sections, ring_settings, settings = (
            self._projection_inputs()
        )
        with (
            patch(
                "calculations.blood_flow_velocity.cross_section."
                "generate_cross_section_signals.label_vessel_branches",
                return_value=branches,
            ),
            patch(
                "calculations.blood_flow_velocity.cross_section."
                "generate_cross_section_signals.section_masks",
                return_value=sections,
            ),
            patch(
                "calculations.blood_flow_velocity.cross_section."
                "generate_cross_section_signals._estimate_orientation",
                return_value=21.0,
            ),
            patch(
                "calculations.blood_flow_velocity.cross_section."
                "generate_cross_section_signals._cross_section_limits",
                return_value=(1, 5),
            ),
        ):
            multi = generate_cross_section_signals_for_cubes(
                reference,
                vessel,
                (10, 10),
                ring_settings,
                settings,
            )

        with (
            patch(
                "calculations.blood_flow_velocity.cross_section."
                "reusable_cross_section_signals.section_masks",
                return_value=sections,
            ),
            patch(
                "calculations.blood_flow_velocity.cross_section."
                "reusable_cross_section_signals._cross_section_limits",
                return_value=(2, 4),
            ) as cross_section_limits,
        ):
            projected = project_cross_section_cube(
                reference,
                multi.plan,
                limits_mode="per_cube",
            )

        self.assertEqual(1, cross_section_limits.call_count)
        np.testing.assert_array_equal(
            projected.profile_rotation_degrees,
            [[21.0]],
        )
        np.testing.assert_array_equal(
            projected.profile_window_bounds_xyxy,
            multi.plan.profile_window_bounds_xyxy,
        )
        np.testing.assert_array_equal(
            projected.profile_integration_limits_pixels,
            [[[2, 4]]],
        )
        np.testing.assert_array_equal(
            multi.plan.profile_integration_limits_pixels,
            [[[1, 5]]],
        )

    def test_registered_cube_spatial_shape_is_enforced(self) -> None:
        reference, vessel, branches, sections, ring_settings, settings = (
            self._projection_inputs()
        )
        with (
            patch(
                "calculations.blood_flow_velocity.cross_section."
                "generate_cross_section_signals.label_vessel_branches",
                return_value=branches,
            ),
            patch(
                "calculations.blood_flow_velocity.cross_section."
                "generate_cross_section_signals.section_masks",
                return_value=sections,
            ),
            patch(
                "calculations.blood_flow_velocity.cross_section."
                "generate_cross_section_signals._estimate_orientation",
                return_value=0.0,
            ),
        ):
            multi = generate_cross_section_signals_for_cubes(
                reference,
                vessel,
                (10, 10),
                ring_settings,
                settings,
            )

        with self.assertRaisesRegex(ValueError, "spatial shape"):
            project_cross_section_cube(
                np.ones((4, 20, 21), dtype=np.float32),
                multi.plan,
            )

    @staticmethod
    def _projection_inputs():
        frame_count = 4
        vessel = np.zeros((21, 21), dtype=bool)
        vessel[8:13, 9:12] = True
        labels = vessel.astype(np.int32)
        branches = SimpleNamespace(
            labels=labels,
            branch_ids=np.asarray([1], dtype=np.int32),
            stages=SimpleNamespace(),
        )
        sections = vessel[None, ...]
        reference = np.zeros((frame_count, 21, 21), dtype=np.float32)
        transverse = np.asarray([1.0, 3.0, 1.0], dtype=np.float32)
        for frame_index in range(frame_count):
            reference[frame_index, 8:13, 9:12] = (
                np.float32(frame_index + 1) * transverse[None, :]
            )
        ring_settings = SimpleNamespace(
            inner_radius_frac=0.0,
            outer_radius_frac=0.5,
            ring_width_frac=0.5,
            ring_count=1,
            segment_length_frac=0.5,
        )
        settings = CrossSectionSignalSettings(30.0, True, 0.5, False, 0.01)
        return reference, vessel, branches, sections, ring_settings, settings


if __name__ == "__main__":
    unittest.main()
