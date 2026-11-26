"""Tests for map-independent topology interpolation and rotation."""

from __future__ import annotations

import unittest
from unittest.mock import patch

import numpy as np
import pytest
from scipy import ndimage as ndi

from calculations.topology.branch_identity import BranchIdentityResult
from calculations.topology.geometry import annulus_mask
from calculations.topology.segments import SegmentTopology
from calculations.topology.orientation import determine_segment_rotations
from calculations.vessel_segments.sampling import transforms
from calculations.vessel_segments.sampling.transforms import (
    dilate_segment_masks,
    interpolate_segment_masks,
    interpolate_segments,
    resample_rotate_segment,
)


def _original_cpu_affine(values, matrix, offset, output_shape):
    """Original N-D SciPy call, including interpolation and border settings."""
    return ndi.affine_transform(
        values, matrix, offset, output_shape=output_shape,
        order=1, mode="grid-constant", cval=0.0, prefilter=False,
    )


@pytest.mark.parametrize("shape,side,angle", [
    ((9, 9), 128, 31.7),
    ((1, 9, 9), 9, 0.0),
    ((4, 9, 9), 13, -23.5),
    ((4, 9, 9), 5, 90.0),
    ((3, 2, 9, 9), 17, 127.3),
    ((2, 1, 1), 4, 180.0),
    ((3, 9, 9), 1, -90.0),
    ((2, 9, 9), 16, 89.999),
])
@pytest.mark.parametrize("missing", ["none", "shared", "changing", "all", "masked", "zeros"])
def test_cpu_weighted_affine_matches_original_nd_call(shape, side, angle, missing):
    rng = np.random.default_rng(117)
    values = rng.normal(size=shape).astype(np.float64 if side == 5 else np.float32)
    if missing == "shared":
        values[..., 0, :] = np.nan
        values[..., -1, -1] = np.inf
    elif missing == "changing":
        for frame, index in enumerate(np.ndindex(shape[:-2])):
            values[index][frame % shape[-2], :] = np.nan
            values[index][-1, frame % shape[-1]] = -np.inf
    elif missing == "all":
        values.fill(np.nan)
    elif missing == "masked":
        y, x = np.indices(shape[-2:])
        values[..., (x - shape[-1] / 2)**2 + (y - shape[-2] / 2)**2 > (shape[-1] / 3)**2] = np.nan
    elif missing == "zeros":
        values.fill(0.0)
        values[..., ::2, ::2] = -0.0
    # Exercise read-only, noncontiguous input without allowing either path to edit it.
    values = values[..., ::-1, ::-1]
    values.flags.writeable = False
    original = values.copy()
    with patch.object(transforms, "optional_cupy_backend", return_value=None):
        with patch.object(transforms, "_affine_transform_cpu", side_effect=_original_cpu_affine):
            expected = resample_rotate_segment(values, angle, side)
        actual = resample_rotate_segment(values, angle, side)
    assert actual.dtype == expected.dtype == np.dtype(np.float32)
    assert actual.shape == expected.shape
    # Check bit patterns, including signed zeros and the NaN mask.
    np.testing.assert_array_equal(actual.view(np.uint32), expected.view(np.uint32))
    np.testing.assert_array_equal(values, original)


def test_cpu_spatial_affine_uses_views_into_one_output_allocation():
    values = np.arange(3 * 2 * 7 * 7, dtype=np.float32).reshape(3, 2, 7, 7)
    values = values[..., ::-1, ::-1]
    values.flags.writeable = False
    matrix, offset = transforms._fused_affine_mapping(4, 7, 11, 16, 31.7)
    output_shape = (3, 2, 16, 16)
    expected = _original_cpu_affine(values, matrix, offset, output_shape)
    original_affine = ndi.affine_transform
    output_planes = []

    def sample(plane, *args, **kwargs):
        assert plane.ndim == 2
        assert np.shares_memory(plane, values)
        output_planes.append(kwargs["output"])
        return original_affine(plane, *args, **kwargs)

    with patch.object(transforms.ndi, "affine_transform", side_effect=sample):
        actual = transforms._affine_transform_cpu(values, matrix, offset, output_shape)
    assert len(output_planes) == 6
    assert actual.flags.owndata
    assert all(np.shares_memory(plane, actual) for plane in output_planes)
    np.testing.assert_array_equal(actual.view(np.uint32), expected.view(np.uint32))


@pytest.mark.parametrize("change", [
    "time_shift", "time_scale", "time_from_y", "y_from_time", "tiny_coupling", "frame_count",
])
def test_cpu_affine_retains_nd_path_for_temporal_transform(change):
    values = np.arange(4 * 7 * 7, dtype=np.float32).reshape(4, 7, 7)
    matrix, offset = transforms._fused_affine_mapping(3, 7, 11, 16, 23.0)
    output_shape = (4, 16, 16)
    if change == "time_shift":
        offset[0] = 0.25
    elif change == "time_scale":
        matrix[0, 0] = 0.75
    elif change == "time_from_y":
        matrix[0, 1] = 0.1
    elif change == "y_from_time":
        matrix[1, 0] = 0.1
    elif change == "tiny_coupling":
        matrix[1, 0] = 1e-15
    else:
        output_shape = (6, 16, 16)
    expected = _original_cpu_affine(values, matrix, offset, output_shape)
    with patch.object(transforms.ndi, "affine_transform", wraps=ndi.affine_transform) as sample:
        actual = transforms._affine_transform_cpu(values, matrix, offset, output_shape)
    sample.assert_called_once()
    assert sample.call_args.args[0].ndim == 3
    np.testing.assert_array_equal(actual.view(np.uint32), expected.view(np.uint32))


def test_cpu_spatial_affine_preserves_empty_batch_shape():
    values = np.empty((0, 7, 7), dtype=np.float32)
    matrix, offset = transforms._fused_affine_mapping(3, 7, 11, 16, 23.0)
    expected = _original_cpu_affine(values, matrix, offset, (0, 16, 16))
    actual = transforms._affine_transform_cpu(values, matrix, offset, (0, 16, 16))
    assert actual.shape == expected.shape == (0, 16, 16)
    assert actual.dtype == expected.dtype


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
