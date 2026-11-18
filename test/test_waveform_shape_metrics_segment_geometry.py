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

from calculations.vessel_segments.measurement import (  # noqa: E402
    MaskedArrays,
    SegmentMeasurementSettings,
)
from calculations.topology import AnnulusGeometry, OpticDisc, section_masks  # noqa: E402
from calculations.topology.branch_identity import _branch_identity_stages  # noqa: E402
from pipelines.velocity_analysis.artifacts.branch_identity import (  # noqa: E402
    _label_image,
    _labels_with_ring_overlay,
    _labels_with_substack_boxes,
)
from pipelines.topology_core.annulus_visualization import (  # noqa: E402
    annulus_gray_levels, annulus_outline_image,
)
from pipelines.velocity_analysis.artifacts.cross_sections import (  # noqa: E402
    export_rotated_mean_pngs,
)
from pipelines.velocity_analysis.analysis.segments import (  # noqa: E402
    analyze_velocity_segment_profiles,
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
        settings = AnnulusGeometry(
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
        ring_settings = AnnulusGeometry(
            inner_radius_frac=0.1,
            outer_radius_frac=0.7,
            ring_width_frac=0.25,
            ring_count=2,
            segment_length_frac=None,
        )
        settings = SegmentMeasurementSettings(0.01, 1.0)

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
            profile = result.profile
            self.assertGreater(profile.topology.native.branch_ids.size, 0)
            self.assertIsNone(profile.maps)
            fft = result.require_fft()
            self.assertEqual(
                (181, 2, 1, profile.topology.native.branch_ids.size, 2),
                fft.unmasked.shape,
            )
            self.assertEqual(
                fft.unmasked.shape,
                fft.masked.shape,
            )
            valid = profile.topology.valid_segments
            self.assertTrue(np.any(valid))
            np.testing.assert_allclose(profile.segment_signal[valid], 1.0)

    def test_substack_debug_overlay_uses_gray_box_edges(self) -> None:
        labels = np.ones((20, 20), dtype=np.int32)
        ring_settings = AnnulusGeometry(
            inner_radius_frac=0.0,
            outer_radius_frac=1.0,
            ring_width_frac=0.5,
            ring_count=2,
        )
        centers = np.asarray(
            [[[5.0, 5.0], [7.0, 5.0]], [[12.0, 5.0], [np.nan, np.nan]]],
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
                    [[11, 14, 4, 7], [-1, -1, -1, -1]],
                ],
                dtype=np.int32,
            ),
        )

        self.assertEqual((20, 20, 3), image.shape)
        np.testing.assert_array_equal(image[4, 4], [128, 128, 128])
        np.testing.assert_array_equal(image[4, 6], [128, 128, 128])
        np.testing.assert_array_equal(image[4, 11], [128, 128, 128])
        base = _labels_with_ring_overlay(labels, (10.0, 10.0), ring_settings)
        np.testing.assert_array_equal(image[5, 5], base[5, 5])

    def test_annulus_outlines_alternate_white_and_gray_at_half_width_stride(self):
        shape = (101, 121)
        disc = OpticDisc(None, (55, 40), 12, 12)
        settings = disc.annulus_geometry(shape, annulus_size=20)
        labels = np.ones(shape, dtype=np.int32)
        outlines = annulus_outline_image(shape, disc.center, settings)
        image = _labels_with_ring_overlay(labels, disc.center, settings)

        np.testing.assert_array_equal(annulus_gray_levels(settings), [255, 128])
        for radius, gray in ((6, 255), (16, 128), (26, 255), (36, 128), (46, 255)):
            self.assertEqual(gray, outlines[40, 55 + radius])
            np.testing.assert_array_equal(image[40, 55 + radius], [gray] * 3)
        np.testing.assert_array_equal(image[outlines == 0], _label_image(labels)[outlines == 0])

    def test_annulus_shades_follow_overlap_phase_count_and_shared_clipped_boundaries(self):
        disc = OpticDisc(None, (50, 50), 12, 12)
        for stride, expected in ((20, [255]), (30, [255]), (10, [255, 128]),
                                 (5, [255, 213, 170, 128])):
            with self.subTest(stride=stride):
                settings = disc.annulus_geometry((101, 101), annulus_size=20, annulus_stride=stride)
                np.testing.assert_array_equal(annulus_gray_levels(settings), expected)

        settings = AnnulusGeometry(.1, .5, .35, 4, ring_stride_frac=.1)
        outlines = annulus_outline_image((101, 101), (50, 50), settings)
        self.assertEqual(213, outlines[50, 85])  # First section clipped to R sets its shade.
        self.assertEqual(0, outlines[50, 99])

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
            profile=SimpleNamespace(
                mean_images=MaskedArrays(rotated_means, -rotated_means),
                topology=SimpleNamespace(
                    valid_segments=np.asarray([[True, False], [True, True]]),
                    native=SimpleNamespace(branch_ids=np.asarray([4, 9])),
                ),
            ),
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

    def test_segment_analysis_uses_disc_extent_and_half_diagonal_spacing(self) -> None:
        settings = OpticDisc(None, (267.0, 230.0), 55, 69).annulus_geometry(
            (512, 512)
        )

        corner_radius = np.hypot(255.5, 255.5)
        farthest_radius = np.hypot(267, 281)
        self.assertEqual(int(np.ceil((farthest_radius - 28) / (corner_radius / 32))),
                         settings.ring_count)
        self.assertAlmostEqual(
            28.0 / corner_radius,
            settings.inner_radius_frac,
        )
        self.assertAlmostEqual(farthest_radius / corner_radius, settings.outer_radius_frac, places=6)
        expected_width = 2 / 32
        self.assertAlmostEqual(expected_width, settings.ring_width_frac)
        self.assertAlmostEqual(expected_width, settings.segment_length_frac)
        self.assertAlmostEqual(expected_width / 2, settings.stride_frac)

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

    def test_annulus_starts_at_optic_disc_and_uses_half_diagonal_extent(self) -> None:
        center = np.asarray([267.0, 230.0])
        settings = OpticDisc(None, center, 55, 69).annulus_geometry((512, 512))
        mask = section_masks(
            (512, 512),
            center,
            settings,
        )[0]

        # The first annulus starts at ceil(min(width, height) / 2) = 28 px.
        self.assertFalse(mask[230, 294])
        self.assertTrue(mask[230, 295])
        self.assertFalse(mask[257, 267])
        self.assertTrue(mask[258, 267])

        # Width is twice the half-diagonal / 32 stride: about 22.58 pixels.
        self.assertTrue(mask[230, 305])
        self.assertTrue(mask[230, 315])
        self.assertTrue(mask[230, 317])
        self.assertFalse(mask[230, 318])
        self.assertFalse(mask[281, 267])

    def test_annulus_pixel_radius_is_independent_of_optic_disc_position(self) -> None:
        settings = AnnulusGeometry(
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
