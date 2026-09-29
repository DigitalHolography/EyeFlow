"""Tests for map-independent topology interpolation and rotation."""

from __future__ import annotations

import unittest
from unittest.mock import patch

import numpy as np

from calculations.compute_backend import optional_cupy_backend
from calculations.topology.geometry import annulus_mask
from calculations.topology.segments import SegmentTopology
from calculations.topology.transforms import (
    determine_segment_rotations,
    dilate_segment_masks,
    interpolate_segment_masks,
    interpolate_segments,
    resample_rotate_segment,
    rotate_segment_masks,
    rotate_segments,
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
            spatial_shape=shape,
            optic_disc_center_xy=center_xy,
            labels=labels,
            centerline=centerline,
            branch_ids=np.asarray([1], dtype=np.int32),
            annulus_masks=annuli,
            segment_masks=np.ones((1, 1, 9, 9), dtype=bool),
            segment_centers_xy=np.asarray([[[12.0, 10.0]]], dtype=np.float32),
            window_bounds_xyxy=np.asarray([[[8, 17, 6, 15]]], dtype=np.int32),
            window_side_pixels=9,
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

    def test_mask_dilation_treats_each_segment_independently(self) -> None:
        masks = np.zeros((1, 2, 7, 7), dtype=bool)
        masks[0, 0, 3, 3] = True
        masks[0, 1, 0, 0] = True

        dilated = dilate_segment_masks(masks, iterations=2)

        self.assertEqual(masks.shape, dilated.shape)
        self.assertEqual(np.bool_, dilated.dtype)
        self.assertEqual(25, np.count_nonzero(dilated[0, 0]))
        self.assertEqual(9, np.count_nonzero(dilated[0, 1]))
        self.assertFalse(dilated[0, 0, 0, 0])
        self.assertFalse(dilated[0, 1, 4, 4])

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

    def test_rotation_occurs_after_interpolation_on_a_larger_canvas(self) -> None:
        values = np.arange(18, dtype=np.float32).reshape(1, 1, 2, 3, 3)
        masks = np.ones((1, 1, 3, 3), dtype=bool)
        rotations = np.asarray([[0.0]], dtype=np.float32)

        interpolated_values = interpolate_segments(values)
        interpolated_masks = interpolate_segment_masks(masks)
        rotated_values = rotate_segments(interpolated_values, rotations)
        rotated_masks = rotate_segment_masks(interpolated_masks, rotations)

        self.assertEqual((1, 1, 2, 128, 128), interpolated_values.shape)
        self.assertEqual((1, 1, 128, 128), interpolated_masks.shape)
        self.assertEqual((1, 1, 2, 181, 181), rotated_values.shape)
        self.assertEqual((1, 1, 181, 181), rotated_masks.shape)
        np.testing.assert_allclose(
            rotated_values[..., 26:154, 26:154],
            interpolated_values,
            equal_nan=True,
        )
        np.testing.assert_array_equal(
            rotated_masks[..., 26:154, 26:154],
            interpolated_masks,
        )

    def test_fused_transform_preserves_alignment_with_bounded_scientific_delta(
        self,
    ) -> None:
        y, x = np.mgrid[-1:1:29j, -1:1:29j]
        base = np.exp(-4.0 * (x**2 + y**2)).astype(np.float32)
        values = np.stack(
            [(1.0 + 0.15 * frame) * base for frame in range(4)]
        ).astype(np.float32)
        angle = np.float32(27.0)

        legacy = rotate_segments(
            interpolate_segments(values[None, None], 128),
            np.asarray([[angle]], dtype=np.float32),
        )[0, 0]
        fused = resample_rotate_segment(values, float(angle), 128)

        self.assertEqual(legacy.shape, fused.shape)
        legacy_y, legacy_x = np.nonzero(np.isfinite(legacy[0]))
        fused_y, fused_x = np.nonzero(np.isfinite(fused[0]))
        np.testing.assert_allclose(
            (np.mean(fused_y), np.mean(fused_x)),
            (np.mean(legacy_y), np.mean(legacy_x)),
            atol=0.1,
        )
        self.assertLessEqual(
            abs(int(np.isfinite(fused[0]).sum()) - int(np.isfinite(legacy[0]).sum())),
            16,
        )

        overlap = np.isfinite(legacy) & np.isfinite(fused)
        pixel_rmse = np.sqrt(np.mean((legacy[overlap] - fused[overlap]) ** 2))
        pixel_scale = np.sqrt(np.mean(legacy[overlap] ** 2))
        pixel_nrmse = float(pixel_rmse / pixel_scale)
        self.assertGreater(pixel_nrmse, 1e-5)
        self.assertLess(pixel_nrmse, 0.002)

        legacy_signal = np.nanmean(legacy, axis=(-2, -1))
        fused_signal = np.nanmean(fused, axis=(-2, -1))
        signal_relative_error = np.max(
            np.abs(legacy_signal - fused_signal)
            / np.maximum(np.abs(legacy_signal), 1e-6)
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

        legacy_profile = transverse_profile(legacy)
        fused_profile = transverse_profile(fused)
        profile_overlap = np.isfinite(legacy_profile) & np.isfinite(fused_profile)
        profile_rmse = np.sqrt(
            np.mean(
                (legacy_profile[profile_overlap] - fused_profile[profile_overlap]) ** 2
            )
        )
        profile_scale = np.sqrt(np.mean(legacy_profile[profile_overlap] ** 2))
        self.assertLess(float(profile_rmse / profile_scale), 0.003)

    def test_invalid_rotation_keeps_segment_outputs_empty(self) -> None:
        values = np.ones((1, 1, 4, 4), dtype=np.float32)
        masks = np.ones((1, 1, 4, 4), dtype=bool)
        rotations = np.asarray([[np.nan]], dtype=np.float32)

        rotated_values = rotate_segments(values, rotations)
        rotated_masks = rotate_segment_masks(masks, rotations)

        self.assertTrue(np.all(np.isnan(rotated_values)))
        self.assertFalse(np.any(rotated_masks))

    def test_cupy_fused_transform_does_not_fall_back_to_scipy(self) -> None:
        if optional_cupy_backend() is None:
            self.skipTest("CuPy/CUDA is unavailable.")
        rng = np.random.default_rng(91)
        values = rng.normal(size=(8, 29, 29)).astype(np.float32)
        values[:, :4, :] = np.nan
        values[:, 20:, 25:] = np.nan

        with patch(
            "calculations.topology.transforms.ndi.affine_transform",
            side_effect=AssertionError("unexpected CPU fallback"),
        ):
            gpu = resample_rotate_segment(values, 31.7)
        with patch(
            "calculations.topology.transforms.optional_cupy_backend",
            return_value=None,
        ):
            cpu = resample_rotate_segment(values, 31.7)

        self.assertEqual(cpu.shape, gpu.shape)
        self.assertEqual(cpu.dtype, gpu.dtype)
        np.testing.assert_array_equal(np.isnan(gpu), np.isnan(cpu))
        np.testing.assert_allclose(gpu, cpu, rtol=5e-5, atol=2e-5, equal_nan=True)
        np.testing.assert_allclose(
            np.nanmean(gpu, axis=(-2, -1)),
            np.nanmean(cpu, axis=(-2, -1)),
            rtol=1e-5,
            atol=1e-6,
        )


if __name__ == "__main__":
    unittest.main()
