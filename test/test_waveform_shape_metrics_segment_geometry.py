"""Tests for spatial segment geometry used to build segment waveforms."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from calculations.vessel_segments.measurement import SegmentMeasurementSettings  # noqa: E402
from calculations.topology import AnnulusGeometry, OpticDisc, section_masks
from pipelines.topology_core.annulus_visualization import annulus_gray_levels, annulus_outline_image
from pipelines.velocity_analysis.artifacts.branch_identity import (
    _label_image,
    _labels_with_ring_overlay,
    _labels_with_substack_boxes,
)  # noqa: E402
from pipelines.velocity_analysis.analysis.segments import (  # noqa: E402
    analyze_velocity_segment_profiles,
)
class SegmentCenterTests(unittest.TestCase):
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

if __name__ == "__main__":
    unittest.main()
