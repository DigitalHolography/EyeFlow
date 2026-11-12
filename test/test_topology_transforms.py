"""Tests for map-independent topology interpolation and rotation."""

from __future__ import annotations

import unittest

import numpy as np

from calculations.topology.branch_identity import BranchIdentityResult
from calculations.topology.geometry import annulus_mask
from calculations.topology.segments import SegmentTopology
from calculations.topology.orientation import determine_segment_rotations
from calculations.vessel_segments.sampling.transforms import (
    dilate_segment_masks,
    interpolate_segment_masks,
    interpolate_segments,
    resample_rotate_segment,
)


class TestTopologyTransforms(unittest.TestCase):
    def test_rotation_uses_a_truncated_multiline_centerline(self) -> None:
        shape = (21, 21)
        center_xy = (3.0, 10.0)
        labels = np.zeros(shape, dtype=np.int32)
        labels[8:13, 8:16] = 1
        centerline = np.zeros(shape, dtype=bool)
        centerline[10, 8:11] = True
        centerline[10, 13:16] = True
        centerline[9, 10:14] = True
        centerline[11, 10:14] = True
        centerline[9:12, 10] = True
        centerline[9:12, 13] = True
        annuli = annulus_mask(shape, center_xy, 0.1, 0.95)[None, ...]
        topology = SegmentTopology(
            optic_disc_center_xy=center_xy,
            branches=BranchIdentityResult(
                labels,
                np.asarray([1], dtype=np.int32),
                centerline,
            ),
            annulus_masks=annuli,
            segment_masks=np.ones((1, 1, 9, 9), dtype=bool),
            segment_centers_xy=np.asarray([[[12.0, 10.0]]], dtype=np.float32),
            window_bounds_xyxy=np.asarray([[[8, 17, 6, 15]]], dtype=np.int32),
        )

        rotations = determine_segment_rotations(topology)

        np.testing.assert_allclose(rotations, [[90.0]])

    def test_values_and_masks_use_separate_interpolation(self) -> None:
        values = np.full((1, 1, 2, 2), np.float32(7.0))
        values[..., 0, 0] = np.nan
        masks = np.asarray([[[[False, True], [True, True]]]], dtype=bool)

        interpolated_values = interpolate_segments(values, output_side_pixels=5)
        interpolated_masks = interpolate_segment_masks(
            masks,
            output_side_pixels=5,
        )

        self.assertEqual((1, 1, 5, 5), interpolated_values.shape)
        self.assertEqual(np.float32, interpolated_values.dtype)
        self.assertEqual((1, 1, 5, 5), interpolated_masks.shape)
        self.assertEqual(np.bool_, interpolated_masks.dtype)
        np.testing.assert_allclose(
            interpolated_values[np.isfinite(interpolated_values)],
            7.0,
        )

    def test_mask_dilation_removes_competing_vessel_pixels(self) -> None:
        masks = np.zeros((1, 1, 7, 7), dtype=bool)
        masks[0, 0, 3, 3] = True
        competing = np.zeros_like(masks)
        competing[0, 0, 3, 3] = True
        competing[0, 0, 3, 4:6] = True

        analysis_mask = dilate_segment_masks(
            masks,
            iterations=2,
            exclusion_masks=competing,
        )

        self.assertTrue(analysis_mask[0, 0, 3, 2])
        self.assertTrue(analysis_mask[0, 0, 3, 3])
        self.assertFalse(analysis_mask[0, 0, 3, 4])
        self.assertFalse(analysis_mask[0, 0, 3, 5])

    def test_fused_transform_preserves_alignment_with_bounded_scientific_delta(
        self,
    ) -> None:
        y, x = np.mgrid[-1:1:29j, -1:1:29j]
        base = np.exp(-4.0 * (x**2 + y**2)).astype(np.float32)
        values = np.stack(
            [(1.0 + 0.15 * frame) * base for frame in range(4)]
        ).astype(np.float32)
        angle = np.float32(27.0)

        staged = resample_rotate_segment(interpolate_segments(values, 128), float(angle), 128)
        fused = resample_rotate_segment(values, float(angle), 128)

        self.assertEqual(staged.shape, fused.shape)
        staged_y, staged_x = np.nonzero(np.isfinite(staged[0]))
        fused_y, fused_x = np.nonzero(np.isfinite(fused[0]))
        np.testing.assert_allclose(
            (np.mean(fused_y), np.mean(fused_x)),
            (np.mean(staged_y), np.mean(staged_x)),
            atol=0.1,
        )
        self.assertLessEqual(
            abs(int(np.isfinite(fused[0]).sum()) - int(np.isfinite(staged[0]).sum())),
            16,
        )

        overlap = np.isfinite(staged) & np.isfinite(fused)
        pixel_rmse = np.sqrt(np.mean((staged[overlap] - fused[overlap]) ** 2))
        pixel_scale = np.sqrt(np.mean(staged[overlap] ** 2))
        pixel_nrmse = float(pixel_rmse / pixel_scale)
        self.assertGreater(pixel_nrmse, 1e-5)
        self.assertLess(pixel_nrmse, 0.002)

        staged_signal = np.nanmean(staged, axis=(-2, -1))
        fused_signal = np.nanmean(fused, axis=(-2, -1))
        signal_relative_error = np.max(
            np.abs(staged_signal - fused_signal)
            / np.maximum(np.abs(staged_signal), 1e-6)
        )
        self.assertLess(float(signal_relative_error), 0.002)

        def transverse_profile(stack: np.ndarray) -> np.ndarray:
            finite = np.isfinite(stack)
            count = np.sum(finite, axis=-2)
            total = np.sum(stack, axis=-2, where=finite)
            return np.divide(
                total,
                count,
                out=np.full(total.shape, np.nan, dtype=np.float32),
                where=count > 0,
            )

        staged_profile = transverse_profile(staged)
        fused_profile = transverse_profile(fused)
        profile_overlap = np.isfinite(staged_profile) & np.isfinite(fused_profile)
        profile_rmse = np.sqrt(
            np.mean(
                (staged_profile[profile_overlap] - fused_profile[profile_overlap]) ** 2
            )
        )
        profile_scale = np.sqrt(np.mean(staged_profile[profile_overlap] ** 2))
        self.assertLess(float(profile_rmse / profile_scale), 0.003)

if __name__ == "__main__":
    unittest.main()
