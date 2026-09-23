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

from calculations.blood_flow_velocity import (  # noqa: E402
    CrossSectionSignalSettings,
)
from calculations.topology.branch_identity import (  # noqa: E402
    _branch_identity_stages,
)
from calculations.blood_flow_velocity.cross_section.generate_cross_section_signals import (
    _dilate_profile_mask,
)
from calculations.topology import (  # noqa: E402
    OpticDisc,
    ring_masks,
    section_masks,
)
from pipelines.waveform_velocity_core.branch_identity_debug import (  # noqa: E402
    _labels_with_substack_boxes,
)
from pipelines.waveform_velocity_core.segments import (  # noqa: E402
    analyze_velocity_segment_profiles,
)
from pipelines.waveform_velocity_core.cross_section_images import (  # noqa: E402
    export_rotated_mean_pngs,
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
            "calculations.topology.branch_identity."
            "annulus_mask",
            return_value=section,
        ), patch(
            "calculations.topology.branch_identity."
            "skeletonize",
            return_value=skeleton,
        ), patch(
            "calculations.topology.branch_identity."
            "_branch_points",
            return_value=np.zeros_like(vessel),
        ), patch(
            "calculations.topology.branch_identity."
            "_remove_small",
            side_effect=lambda mask, _min_area: mask,
        ), patch(
            "calculations.topology.branch_identity."
            "watershed",
            return_value=watershed_labels,
        ):
            stages = _branch_identity_stages(
                vessel,
                (4.0, 4.0),
                settings,
                optic_disc_mask=np.zeros_like(vessel),
            )

        self.assertEqual(2, int(stages.annulus_refined_labels.max()))

    def test_segment_velocity_results_runs_with_prepared_topology(self) -> None:
        velocity = np.ones((3, 61, 61), dtype=np.float32)
        artery_mask = np.zeros((61, 61), dtype=bool)
        artery_mask[27:34, 5:56] = True
        vein_mask = np.zeros_like(artery_mask)
        vein_mask[5:56, 37:44] = True
        optic_disc_mask = np.zeros_like(artery_mask)
        optic_disc_mask[27:34, 27:34] = True
        ring_settings = SimpleNamespace(
            inner_radius_frac=0.1,
            outer_radius_frac=0.7,
            ring_width_frac=0.25,
            ring_count=2,
            segment_length_frac=None,
        )
        settings = CrossSectionSignalSettings(0.01, 1.0)

        results = analyze_velocity_segment_profiles(
            velocity,
            {
                "artery": artery_mask,
                "vein": vein_mask,
            },
            OpticDisc(optic_disc_mask, (30.0, 30.0), None, None),
            ring_settings,
            settings,
            retain_velocity_maps=False,
            cycle_boundary_indexes=np.asarray([0, 2], dtype=np.int32),
            velocity_profile_fft=True,
        )

        for result in results.values():
            self.assertGreater(result.branch_ids.size, 0)
            self.assertIsNone(result.velocity_maps_per_segment)
            self.assertEqual((0, 2), result.velocity_map_segment_indexes.shape)
            self.assertEqual(
                (181, 2, 1, result.branch_ids.size, 2),
                result.transverse_velocity_fft_profiles_unmasked.shape,
            )
            self.assertEqual(
                result.transverse_velocity_fft_profiles_unmasked.shape,
                result.transverse_velocity_fft_profiles_masked.shape,
            )
            valid = result.topology.valid_segments
            self.assertTrue(np.any(valid))
            np.testing.assert_allclose(result.velocity[valid], 1.0)

    def test_profile_mask_dilation_expands_ten_pixels_horizontally(self) -> None:
        mask = np.zeros((51, 51), dtype=bool)
        mask[25, 25] = True

        dilated = _dilate_profile_mask(mask, 10)

        self.assertEqual(np.bool_, dilated.dtype)
        self.assertEqual(21, int(np.count_nonzero(dilated)))
        self.assertTrue(np.all(dilated[25, 15:36]))
        self.assertFalse(np.any(dilated[:25]))
        self.assertFalse(np.any(dilated[26:]))
        self.assertFalse(np.any(dilated[:, :15]))
        self.assertEqual(1, int(np.count_nonzero(mask)))

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
            (10.0, 10.0),
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

    def test_segment_analysis_uses_disc_extent_and_fixed_fov_spacing(self) -> None:
        settings = OpticDisc(None, (267.0, 230.0), 55, 69).annulus_geometry(
            (512, 512)
        )

        corner_radius = np.hypot(255.5, 255.5)
        self.assertEqual(16, settings.ring_count)
        self.assertAlmostEqual(
            (69.0 / 2.0) / corner_radius,
            settings.inner_radius_frac,
        )
        self.assertEqual(1.0, settings.outer_radius_frac)
        expected_width = (512.0 / 25.0) / corner_radius
        self.assertAlmostEqual(expected_width, settings.ring_width_frac)
        self.assertAlmostEqual(expected_width, settings.segment_length_frac)

    def test_number_of_radii_in_fov_controls_pixel_width_and_count(self) -> None:
        settings = OpticDisc(None, (50.0, 50.0), 40, 40).annulus_geometry(
            (101, 101),
            number_of_radii_in_fov=10,
        )

        radius_scale = np.hypot(50.0, 50.0)
        self.assertAlmostEqual(101.0 / 10.0, settings.ring_width_frac * radius_scale)
        self.assertEqual(settings.ring_width_frac, settings.segment_length_frac)
        self.assertEqual(6, settings.ring_count)

    def test_number_of_radii_in_fov_must_be_positive(self) -> None:
        with self.assertRaisesRegex(ValueError, "must be positive"):
            OpticDisc(None, (1.0, 1.0), 1, 1).annulus_geometry(
                (10, 10), number_of_radii_in_fov=0
            )

    def test_annulus_starts_at_optic_disc_and_uses_fixed_fov_extent(self) -> None:
        center = np.asarray([267.0, 230.0])
        settings = OpticDisc(None, center, 55, 69).annulus_geometry((512, 512))
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
        self.assertTrue(mask[230, 321])
        self.assertFalse(mask[230, 322])
        self.assertTrue(mask[284, 267])

        rings = ring_masks((512, 512), center, settings)
        self.assertTrue(rings[-1, 0, 0])
        self.assertFalse(rings[-1, 511, 511])

    def test_annulus_pixel_radius_is_independent_of_optic_disc_position(self) -> None:
        settings = SimpleNamespace(
            inner_radius_frac=0.0,
            outer_radius_frac=0.5,
            ring_width_frac=0.25,
            ring_count=1,
            segment_length_frac=0.25,
        )
        centered = section_masks((101, 101), (50, 50), settings)[0]
        shifted = section_masks((101, 101), (30, 40), settings)[0]

        centered_y, centered_x = np.nonzero(centered)
        shifted_y, shifted_x = np.nonzero(shifted)
        self.assertEqual(centered_x.max() - 50, shifted_x.max() - 30)
        self.assertEqual(centered_y.max() - 50, shifted_y.max() - 40)

if __name__ == "__main__":
    unittest.main()
