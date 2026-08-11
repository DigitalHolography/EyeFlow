"""Tests for spatial segment geometry used to build segment waveforms."""

from __future__ import annotations

import sys
import unittest
import warnings
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from calculations.blood_flow_velocity import (  # noqa: E402
    CrossSectionSignalSettings,
    segment_velocity_results,
)
from calculations.blood_flow_velocity.cross_section.branch_identity import (  # noqa: E402
    _branch_identity_stages,
)
from calculations.blood_flow_velocity.cross_section.generate_cross_section_signals import (  # noqa: E402
    _ROTATED_SUBSTACK_SIDE,
    _CrossSectionBuffers,
    _CrossSectionMeasurement,
    _center_pad_for_rotation,
    _fill_cross_section_buffers,
    _fixed_subimage_stack,
    _fixed_substack_side_pixels,
    _frame_velocities,
    _PreparedCrossSectionGeometry,
    _resize_subimage_stack,
    _resize_submask,
    _rotated_profile_sample_count,
    _rotate_stack_with_nan,
    _sample_nanstd_axis0,
)
from calculations.blood_flow_velocity.cross_section.reusable_cross_section_signals import (  # noqa: E402
    generate_cross_section_signals_for_cubes,
    project_cross_section_cube,
)
from calculations.blood_flow_velocity.cross_section.segment_geometry import (  # noqa: E402
    ring_masks,
    section_masks,
)
from calculations.math import rotate_image_with_nan  # noqa: E402
from pipelines.waveform_velocity_core.branch_identity_debug import (  # noqa: E402
    _labels_with_substack_boxes,
)
from pipelines.waveform_velocity_core.cross_section_images import (  # noqa: E402
    export_rotated_mean_pngs,
)
from pipelines.waveform_velocity_core.runner import _segment_ring_settings  # noqa: E402
from utils.logger import Logger  # noqa: E402


class SegmentCenterTests(unittest.TestCase):
    def tearDown(self) -> None:
        Logger.reset_current()

    def test_sample_nanstd_does_not_warn_for_undersampled_columns(self) -> None:
        values = np.asarray(
            [
                [np.nan, 3.0, 1.0, 2.0],
                [np.nan, np.nan, 3.0, 4.0],
                [np.nan, np.nan, np.nan, 6.0],
            ],
            dtype=np.float32,
        )

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            actual = _sample_nanstd_axis0(values)

        self.assertEqual([], caught)
        np.testing.assert_allclose(
            actual,
            [np.nan, np.nan, np.sqrt(2.0), 2.0],
            rtol=1e-6,
            equal_nan=True,
        )

    def test_frame_velocities_rotate_the_complete_batch_before_projection(self) -> None:
        sub_stack = np.asarray(
            [
                [[1.0, np.nan, 9.0], [3.0, np.nan, np.nan]],
                [[2.0, 4.0, np.nan], [6.0, np.nan, np.nan]],
            ],
            dtype=np.float32,
        )
        rotated_frames = [
            np.asarray([[1.0, 2.0, np.nan], [3.0, 4.0, 5.0]], dtype=np.float32),
            np.asarray([[2.0, 4.0, 6.0], [4.0, 8.0, np.nan]], dtype=np.float32),
        ]

        with patch(
            "calculations.blood_flow_velocity.cross_section."
            "generate_cross_section_signals._rotate_stack_with_nan",
            return_value=np.stack(rotated_frames),
        ) as rotate:
            raw, safe, profiles = _frame_velocities(
                sub_stack,
                angle=90.0,
                c1=1,
                c2=2,
            )

        rotate.assert_called_once()
        np.testing.assert_array_equal(rotate.call_args.args[0], sub_stack)
        self.assertEqual(90.0, rotate.call_args.args[1])
        np.testing.assert_allclose(
            profiles,
            [[2.0, 3.0, 5.0], [3.0, 6.0, 6.0]],
            equal_nan=True,
        )
        np.testing.assert_allclose(raw, [4.0, 6.0])
        np.testing.assert_allclose(safe, [10.0 / 3.0, 5.0])

    def test_batched_nan_rotation_matches_independent_frames(self) -> None:
        rng = np.random.default_rng(5)
        stack = rng.random((4, 17, 17), dtype=np.float32)
        stack[:, :3, :] = np.nan

        expected = np.stack(
            [rotate_image_with_nan(frame, 31.0) for frame in stack],
            axis=0,
        )
        actual = _rotate_stack_with_nan(stack, 31.0)

        np.testing.assert_allclose(actual, expected, rtol=1e-6, equal_nan=True)

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
            "annulus_mask",
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

    def test_fixed_substack_side_uses_higher_joint_percentiles_and_next_odd(self) -> None:
        sections = np.ones((1, 20, 20), dtype=bool)
        artery_labels = np.zeros((20, 20), dtype=np.int32)
        vein_labels = np.zeros((20, 20), dtype=np.int32)
        artery_labels[1:4, 1:5] = 1
        artery_labels[6:11, 1:7] = 2
        vein_labels[1:9, 10:12] = 1
        artery = _PreparedCrossSectionGeometry(
            artery_labels > 0,
            SimpleNamespace(
                labels=artery_labels,
                branch_ids=np.asarray([1, 2], dtype=np.int32),
            ),
            sections,
        )
        vein = _PreparedCrossSectionGeometry(
            vein_labels > 0,
            SimpleNamespace(
                labels=vein_labels,
                branch_ids=np.asarray([1], dtype=np.int32),
            ),
            sections,
        )

        self.assertEqual(9, _fixed_substack_side_pixels((artery, vein), 0.95))

    def test_artery_and_vein_measurements_receive_one_joint_fixed_side(self) -> None:
        sections = np.ones((1, 20, 20), dtype=bool)
        artery_labels = np.zeros((20, 20), dtype=np.int32)
        vein_labels = np.zeros((20, 20), dtype=np.int32)
        artery_labels[1:4, 1:5] = 1
        vein_labels[1:9, 10:12] = 1
        geometries = (
            _PreparedCrossSectionGeometry(
                artery_labels > 0,
                SimpleNamespace(
                    labels=artery_labels,
                    branch_ids=np.asarray([1], dtype=np.int32),
                ),
                sections,
            ),
            _PreparedCrossSectionGeometry(
                vein_labels > 0,
                SimpleNamespace(
                    labels=vein_labels,
                    branch_ids=np.asarray([1], dtype=np.int32),
                ),
                sections,
            ),
        )
        settings = CrossSectionSignalSettings(False, 0.5, False, 0.01)

        with patch(
            "calculations.blood_flow_velocity.cross_section."
            "segment_velocity_signals._prepare_cross_section_geometry",
            side_effect=geometries,
        ), patch(
            "calculations.blood_flow_velocity.cross_section."
            "segment_velocity_signals._generate_cross_section_signals_from_geometry",
            side_effect=("artery", "vein"),
        ) as generate:
            results = segment_velocity_results(
                np.zeros((1, 20, 20), dtype=np.float32),
                artery_labels > 0,
                vein_labels > 0,
                (10, 10),
                SimpleNamespace(ring_count=1),
                settings,
            )

        self.assertEqual(("artery", "vein"), results)
        self.assertEqual([9, 9], [call.args[-1] for call in generate.call_args_list])

    def test_fixed_subimage_is_centroid_centered_and_padded_at_periphery(self) -> None:
        velocity = np.arange(2 * 5 * 6, dtype=np.float32).reshape(2, 5, 6)
        segment = np.ones((5, 6), dtype=bool)

        stack, mask, bounds = _fixed_subimage_stack(
            velocity,
            segment,
            (0, 1),
            5,
        )

        self.assertEqual((2, 5, 5), stack.shape)
        self.assertEqual((5, 5), mask.shape)
        self.assertEqual((0, 3, 0, 4), bounds)
        np.testing.assert_array_equal(mask[:, :2], False)
        np.testing.assert_array_equal(mask[0], False)
        self.assertTrue(np.all(np.isnan(stack[:, :, :2])))
        self.assertTrue(np.all(np.isnan(stack[:, 0])))
        np.testing.assert_array_equal(stack[:, 1:5, 2:5], velocity[:, 0:4, 0:3])

    def test_fixed_subimage_retains_in_bounds_pixels_outside_segment_mask(self) -> None:
        velocity = np.arange(25, dtype=np.float32).reshape(1, 5, 5)
        segment = np.zeros((5, 5), dtype=bool)
        segment[2, 2] = True

        stack, mask, _ = _fixed_subimage_stack(
            velocity,
            segment,
            (2, 2),
            3,
        )

        np.testing.assert_array_equal(stack[0], velocity[0, 1:4, 1:4])
        self.assertEqual(1, int(np.count_nonzero(mask)))
        self.assertTrue(mask[1, 1])

    def test_image_and_mask_resize_use_separate_128_square_paths(self) -> None:
        image = np.full((3, 3), np.nan, dtype=np.float32)
        image[1, 1] = 7.0
        mask = np.zeros((3, 3), dtype=bool)
        mask[1, 1] = True

        resized_image = _resize_subimage_stack(image[None, ...])[0]
        resized_mask = _resize_submask(mask)

        self.assertEqual((128, 128), resized_image.shape)
        self.assertEqual((128, 128), resized_mask.shape)
        self.assertEqual(np.bool_, resized_mask.dtype)
        self.assertTrue(np.allclose(resized_image[resized_mask], 7.0))
        resized_image[~resized_mask] = np.nan
        self.assertTrue(np.all(np.isnan(resized_image[~resized_mask])))

    def test_rotation_canvas_and_unpadded_sample_count_cover_all_angles(self) -> None:
        values = np.ones((2, 128, 128), dtype=np.float32)
        mask = np.ones((128, 128), dtype=bool)

        padded_values = _center_pad_for_rotation(values, np.nan)
        padded_mask = _center_pad_for_rotation(mask, False)

        self.assertEqual(181, _ROTATED_SUBSTACK_SIDE)
        self.assertEqual((2, 181, 181), padded_values.shape)
        self.assertEqual((181, 181), padded_mask.shape)
        np.testing.assert_array_equal(padded_values[:, 26:154, 26:154], values)
        np.testing.assert_array_equal(padded_mask[26:154, 26:154], mask)
        self.assertTrue(np.all(np.isnan(padded_values[:, :26])))
        self.assertFalse(np.any(padded_mask[:26]))
        self.assertEqual(128, _rotated_profile_sample_count(0.0))
        self.assertEqual(181, _rotated_profile_sample_count(45.0))
        self.assertEqual(128, _rotated_profile_sample_count(90.0))
        self.assertEqual(179, _rotated_profile_sample_count(37.0))

    def test_substack_debug_overlay_marks_shared_box_edges(self) -> None:
        labels = np.ones((20, 20), dtype=np.int32)
        ring_settings = SimpleNamespace(
            inner_radius_frac=0.0,
            outer_radius_frac=1.0,
            ring_width_frac=0.5,
            ring_count=2,
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
            np.asarray(
                [
                    [[4, 7, 4, 7], [6, 9, 4, 7]],
                    [[-1, -1, -1, -1], [-1, -1, -1, -1]],
                ],
                dtype=np.int32,
            ),
        )

        self.assertEqual((20, 20, 3), image.shape)
        np.testing.assert_array_equal(image[4, 4], [255, 255, 0])
        np.testing.assert_array_equal(image[4, 6], [255, 0, 255])

    def test_rotated_means_export_to_separate_vessel_subfolders(self) -> None:
        class RecordingOutput:
            def __init__(self) -> None:
                self.writes = []

            def write_png(self, image, filename):
                self.writes.append((np.asarray(image).copy(), filename))
                return Path(filename)

        rotated_means = np.arange(
            2 * 2 * 181 * 181,
            dtype=np.float32,
        ).reshape(2, 2, 181, 181)
        result = SimpleNamespace(
            rotated_mean_images=rotated_means,
            rotated_mean_images_masked=-rotated_means,
            profile_sample_count=np.asarray([[128, 0], [128, 128]]),
            branch_ids=np.asarray([4, 9]),
        )
        output = RecordingOutput()

        artery_paths = export_rotated_mean_pngs(output, result, "arteries")
        vein_paths = export_rotated_mean_pngs(output, result, "veins")

        self.assertEqual(
            [
                Path("rotated_mean/arteries/ring_001_branch_004.png"),
                Path("rotated_mean/arteries/ring_002_branch_004.png"),
                Path("rotated_mean/arteries/ring_002_branch_009.png"),
                Path("rotated_mean_masked/arteries/ring_001_branch_004.png"),
                Path("rotated_mean_masked/arteries/ring_002_branch_004.png"),
                Path("rotated_mean_masked/arteries/ring_002_branch_009.png"),
            ],
            artery_paths,
        )
        self.assertEqual(
            [
                Path("rotated_mean/veins/ring_001_branch_004.png"),
                Path("rotated_mean/veins/ring_002_branch_004.png"),
                Path("rotated_mean/veins/ring_002_branch_009.png"),
                Path("rotated_mean_masked/veins/ring_001_branch_004.png"),
                Path("rotated_mean_masked/veins/ring_002_branch_004.png"),
                Path("rotated_mean_masked/veins/ring_002_branch_009.png"),
            ],
            vein_paths,
        )
        np.testing.assert_array_equal(output.writes[0][0], rotated_means[0, 0])

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
        settings = CrossSectionSignalSettings(False, 0.5, False, 0.01)

        signal = np.arange(3, dtype=np.float32)
        rotated_mean = np.full((181, 181), 7.0, dtype=np.float32)
        rotated_mean_masked = np.full((181, 181), 5.0, dtype=np.float32)
        profiles = np.ones((3, 181), dtype=np.float32)
        profiles_masked = np.full((3, 181), 2.0, dtype=np.float32)
        with patch(
            "calculations.blood_flow_velocity.cross_section."
            "generate_cross_section_signals._tilt_angle",
            return_value=0.0,
        ), patch(
            "calculations.blood_flow_velocity.cross_section."
            "generate_cross_section_signals._cross_section_velocity_from_substack",
            return_value=_CrossSectionMeasurement(
                unmasked=(
                    signal,
                    signal,
                    profiles,
                    0.0,
                    np.zeros((181,), dtype=np.float32),
                ),
                masked=(
                    signal,
                    signal,
                    profiles_masked,
                    0.0,
                    np.zeros((181,), dtype=np.float32),
                ),
                rotated_mean=rotated_mean,
                rotated_mean_masked=rotated_mean_masked,
                limits=(0, 127),
                sample_count=128,
            ),
        ):
            _fill_cross_section_buffers(
                buffers,
                velocity,
                sections,
                branches,
                np.asarray([3.0, 3.0], dtype=np.float32),
                settings,
                1,
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
        self.assertEqual(128, buffers.profile_sample_count[0, 0])
        self.assertEqual(128, buffers.profile_sample_count[1, 1])
        self.assertEqual(
            (2, 2, 3, 181),
            buffers.velocity_profiles.shape,
        )
        np.testing.assert_array_equal(buffers.rotated_mean_images[0, 0], rotated_mean)
        np.testing.assert_array_equal(buffers.rotated_mean_images[1, 1], rotated_mean)
        np.testing.assert_array_equal(buffers.velocity_profiles[0, 0], profiles)
        np.testing.assert_array_equal(
            buffers.velocity_profiles_masked[0, 0],
            profiles_masked,
        )
        np.testing.assert_array_equal(
            buffers.rotated_mean_images_masked[0, 0],
            rotated_mean_masked,
        )


class ReusableCrossSectionProjectionTests(unittest.TestCase):
    def test_named_cubes_reuse_fixed_window_angle_and_limits(self) -> None:
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
                return_value=(80, 100),
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
            [[[8, 13, 8, 13]]],
        )
        self.assertEqual(5, result.plan.profile_window_side_pixels)
        self.assertAlmostEqual(0.01 * 5.0 / 128.0, result.plan.profile_pixel_size_mm)
        self.assertEqual((1, 1, 4, 181), result.reference.velocity_profiles.shape)
        self.assertEqual(
            (1, 1, 4, 181),
            result.reference.velocity_profiles_masked.shape,
        )
        self.assertEqual((1, 1, 181, 181), result.reference.rotated_mean_images.shape)
        self.assertEqual(
            (1, 1, 181, 181),
            result.reference.rotated_mean_images_masked.shape,
        )
        self.assertEqual((181,), result.reference.profile_x_micrometers.shape)
        np.testing.assert_array_equal(result.reference.profile_sample_count, [[179]])
        np.testing.assert_array_equal(
            result.plan.profile_integration_limits_pixels,
            [[[80, 100]]],
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
        settings = CrossSectionSignalSettings(True, 0.5, False, 0.01)
        return reference, vessel, branches, sections, ring_settings, settings


if __name__ == "__main__":
    unittest.main()
