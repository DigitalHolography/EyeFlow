"""Regression coverage for bounded fused and staged segment transforms."""

from __future__ import annotations

import numpy as np
import pytest
from unittest.mock import patch

from calculations.math.spatial_gradient import (
    GAUSSIAN_BLUR_RADIUS,
    UNSHARP_MASK_RADIUS,
)
from calculations.topology import (
    PreparedTopology,
    AnnulusGeometry,
    SegmentTopology,
    interpolate_segment_masks,
    interpolate_segments,
    prepare_segment_chunks,
    resample_rotate_segment,
    resolve_segment_rotations,
    rotate_segment_masks,
)
from pipelines.spatial_gradient_moment0.profiles import (
    _spatial_gradient_chain,
)


def _prepared(*, angle: float = 23.0, centerline_points: int = 9) -> PreparedTopology:
    side = 9
    mask = np.zeros((side, side), dtype=bool)
    mask[:, 3:6] = True
    centerline = np.zeros_like(mask)
    centerline[:centerline_points, 4] = True
    topology = SegmentTopology(
        spatial_shape=(side, side),
        optic_disc_center_xy=(4.0, 8.0),
        labels=mask.astype(np.int32),
        centerline=centerline,
        branch_ids=np.asarray([1], dtype=np.int32),
        annulus_masks=np.ones((1, side, side), dtype=bool),
        segment_masks=mask[None, None],
        segment_centers_xy=np.asarray([[[4.0, 4.0]]], dtype=np.float32),
        window_bounds_xyxy=np.asarray([[[0, side, 0, side]]], dtype=np.int32),
        window_side_pixels=side,
        optic_disc_mask=np.zeros((side, side), dtype=bool),
        ring_settings=AnnulusGeometry(0.0, 1.0, 1.0, 1),
    )
    rotations = np.asarray([[angle]], dtype=np.float32)
    interpolated = interpolate_segment_masks(topology.segment_masks)
    return PreparedTopology(
        topology=topology,
        rotation_degrees=rotations,
        interpolated_masks=interpolated,
        rotated_masks=rotate_segment_masks(interpolated, rotations),
    )


def _assembled(chunks) -> np.ndarray:
    values = list(chunks)
    assert values
    assert {(item.ring_index, item.branch_index) for item in values} == {(0, 0)}
    return np.concatenate([np.asarray(item.rotated) for item in values], axis=0)


def test_fused_chunks_match_one_chunk_across_budget_boundary() -> None:
    prepared = _prepared()
    rng = np.random.default_rng(90210)
    cube = rng.normal(size=(11, 9, 9)).astype(np.float32)
    cube[2, 1, 1] = np.nan
    whole = _assembled(
        prepare_segment_chunks(cube, prepared, working_memory_mb=64.0)
    )
    chunked = _assembled(
        prepare_segment_chunks(cube, prepared, working_memory_mb=0.30)
    )
    np.testing.assert_allclose(chunked, whole, rtol=1e-6, atol=1e-6, equal_nan=True)


def test_masked_companion_applies_mask_after_interpolation_before_rotation() -> None:
    prepared = _prepared()
    cube = np.arange(3 * 9 * 9, dtype=np.float32).reshape(3, 9, 9)
    chunks = list(
        prepare_segment_chunks(
            cube,
            prepared,
            working_memory_mb=64.0,
            include_masked_before_rotation=True,
        )
    )
    actual = np.concatenate(
        [np.asarray(chunk.rotated_masked) for chunk in chunks],
        axis=0,
    )
    interpolated = interpolate_segments(cube)
    interpolated[..., ~prepared.interpolated_masks[0, 0]] = np.nan
    expected = resample_rotate_segment(
        interpolated,
        float(prepared.rotation_degrees[0, 0]),
    )
    np.testing.assert_allclose(
        actual,
        expected,
        rtol=1e-6,
        atol=1e-6,
        equal_nan=True,
    )


@pytest.mark.parametrize("frame_count", [5, 17])
def test_staged_gradient_chunks_match_whole_recording(frame_count: int) -> None:
    prepared = _prepared()
    rng = np.random.default_rng(frame_count)
    cube = rng.uniform(0.5, 4.0, size=(frame_count, 9, 9)).astype(np.float32)
    cube[frame_count // 2, 4, 4] = np.nan
    whole = _assembled(
        prepare_segment_chunks(
            cube,
            prepared,
            working_memory_mb=64.0,
            transform_mode="staged",
            post_interpolation=_spatial_gradient_chain,
            temporal_halo=6,
            scratch_array_count=7,
        )
    )
    chunked = _assembled(
        prepare_segment_chunks(
            cube,
            prepared,
            working_memory_mb=7.0 if frame_count > 6 else 3.0,
            transform_mode="staged",
            post_interpolation=_spatial_gradient_chain,
            temporal_halo=6,
            scratch_array_count=7,
        )
    )
    np.testing.assert_allclose(chunked, whole, rtol=1e-5, atol=1e-5, equal_nan=True)


def test_staged_gradient_rejects_budget_below_one_frame_context() -> None:
    with pytest.raises(MemoryError, match="one output frame"):
        list(
            prepare_segment_chunks(
                np.ones((17, 9, 9), dtype=np.float32),
                _prepared(),
                working_memory_mb=1.0,
                transform_mode="staged",
                post_interpolation=_spatial_gradient_chain,
                temporal_halo=6,
                scratch_array_count=7,
            )
        )


def test_image_projection_fallback_resolves_once_and_zero_reference_does_not() -> None:
    prepared = _prepared(angle=np.nan, centerline_points=1)
    reference = np.zeros((4, 9, 9), dtype=np.float32)
    reference[:, :, 3:6] = 3.0
    resolved = resolve_segment_rotations(prepared, reference)
    assert np.isfinite(resolved.rotation_degrees[0, 0])
    assert resolve_segment_rotations(resolved, reference) is resolved
    unresolved = resolve_segment_rotations(prepared, np.zeros_like(reference))
    assert np.isnan(unresolved.rotation_degrees[0, 0])


def test_documented_imagej_radii_are_authoritative() -> None:
    assert GAUSSIAN_BLUR_RADIUS == 6.0
    assert UNSHARP_MASK_RADIUS == 8.0


def test_spatial_gradient_chain_has_the_exact_scientific_order() -> None:
    calls: list[str] = []

    def record(name):
        def operation(values):
            calls.append(name)
            return values

        return operation

    module = "pipelines.spatial_gradient_moment0.profiles"
    with (
        patch(f"{module}.moving_avg_window", side_effect=record("moving_average")),
        patch(f"{module}.sobel_spatial_gradient", side_effect=record("sobel")),
        patch(f"{module}.gaussian2d_blur", side_effect=record("gaussian_6")),
        patch(f"{module}.unsharpen", side_effect=record("unsharp_8_weight_0.6")),
    ):
        _spatial_gradient_chain(np.ones((3, 5, 5), dtype=np.float32))

    assert calls == [
        "moving_average",
        "sobel",
        "gaussian_6",
        "unsharp_8_weight_0.6",
        "moving_average",
    ]
