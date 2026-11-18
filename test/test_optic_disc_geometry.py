"""Tests for the authoritative optic-disc geometry object."""

from __future__ import annotations

import numpy as np
import pytest

from calculations.topology import AnnulusGeometry, OpticDisc
from calculations.topology.geometry import image_half_diagonal


def test_validates_center_dimensions_and_mask() -> None:
    with pytest.raises(ValueError, match="two finite"):
        OpticDisc(None, (np.nan, 2.0), 4.0, 6.0)
    with pytest.raises(ValueError, match="both be present or absent"):
        OpticDisc(None, (1.0, 2.0), 4.0, None)
    with pytest.raises(ValueError, match="positive scalar"):
        OpticDisc(None, (1.0, 2.0), 0.0, 6.0)
    with pytest.raises(ValueError, match="2-D"):
        OpticDisc(np.zeros((2, 3, 1)), (1.0, 2.0), 4.0, 6.0)


def test_supplied_mask_takes_precedence_and_is_owned_by_the_disc() -> None:
    supplied = np.zeros((5, 7), dtype=bool)
    supplied[1, 2] = True
    disc = OpticDisc(supplied, (3.0, 2.0), 100.0, 100.0)
    supplied[1, 2] = False

    resolved = disc.mask_for((5, 7))

    assert resolved.sum() == 1
    assert resolved[1, 2]
    resolved[1, 2] = False
    assert disc.mask_for((5, 7))[1, 2]


def test_reconstructs_ellipse_and_requires_dimensions_only_when_needed() -> None:
    disc = OpticDisc(None, (3.0, 2.0), 4.0, 2.0)
    mask = disc.mask_for((5, 7))

    assert mask[2, 3]
    assert mask[2, 1]
    assert mask[1, 3]
    assert not mask[1, 1]

    mask_only = OpticDisc(np.zeros((5, 7), dtype=bool), (3.0, 2.0), None, None)
    assert not mask_only.mask_for((5, 7)).any()
    with pytest.raises(ValueError, match="required to reconstruct"):
        OpticDisc(None, (3.0, 2.0), None, None).mask_for((5, 7))
    with pytest.raises(ValueError, match="required to derive"):
        mask_only.annulus_geometry((5, 7))


def test_transpose_changes_mask_and_all_geometry_atomically() -> None:
    mask = np.zeros((2, 3), dtype=bool)
    mask[0, 2] = True
    transposed = OpticDisc(mask, (2.0, 1.0), 6.0, 4.0).transposed()

    assert transposed.mask.shape == (3, 2)
    assert transposed.mask[2, 0]
    assert transposed.center == (1.0, 2.0)
    assert transposed.width == 4.0
    assert transposed.height == 6.0

    square = np.eye(3, dtype=bool)
    square_disc = OpticDisc(square, (2.0, 1.0), 6.0, 4.0).transposed()
    np.testing.assert_array_equal(square_disc.mask, square.T)
    assert square_disc.center == (1.0, 2.0)
    assert (square_disc.width, square_disc.height) == (4.0, 6.0)


def test_shape_mismatch_is_rejected_without_implicit_transpose() -> None:
    disc = OpticDisc(np.zeros((2, 3), dtype=bool), (1.0, 1.0), None, None)
    with pytest.raises(ValueError, match="same orientation and shape"):
        disc.mask_for((3, 2))


def test_subtraction_is_boolean_exact_and_does_not_mutate_input() -> None:
    disc_mask = np.zeros((4, 5), dtype=bool)
    disc_mask[1, 2] = True
    disc = OpticDisc(disc_mask, (2.0, 1.0), None, None)
    vessel = np.ones((4, 5), dtype=bool)
    original = vessel.copy()

    result = disc.subtract_from(vessel)

    np.testing.assert_array_equal(vessel, original)
    assert result is not vessel
    assert not result[1, 2]
    assert result.sum() == vessel.sum() - 1
    with pytest.raises(TypeError, match="2-D boolean"):
        disc.subtract_from(vessel.astype(np.uint8))
    with pytest.raises(TypeError, match="2-D boolean"):
        disc.subtract_from(vessel[None, ...])


def test_annulus_geometry_preserves_established_calculation() -> None:
    disc = OpticDisc(None, (7.0, 4.0), 6.0, 4.0)
    geometry = disc.annulus_geometry((10, 20), number_of_radii_in_fov=10)
    scale = image_half_diagonal(10, 20)
    expected_step = 20.0 / 10.0 / scale
    expected_inner = 2.0 / scale

    assert isinstance(geometry, AnnulusGeometry)
    assert geometry.inner_radius_frac == pytest.approx(expected_inner)
    assert geometry.outer_radius_frac == 1.0
    assert geometry.ring_width_frac == pytest.approx(expected_step)
    assert geometry.segment_length_frac == pytest.approx(expected_step)
    assert geometry.ring_count == int(np.ceil((1.0 - expected_inner) / expected_step))
    with pytest.raises(ValueError, match="positive"):
        disc.annulus_geometry((10, 20), number_of_radii_in_fov=0)
    assert OpticDisc(None, (0.0, 0.0), 1.0, 1.0).annulus_geometry((1, 1)).inner_radius_frac == 1.0


def test_default_annuli_insert_exact_intermediate_support_at_alternate_indices():
    from dataclasses import replace
    from calculations.topology.geometry import section_bounds, section_masks

    shape = (512, 512)
    disc = OpticDisc(None, (267.0, 230.0), 55, 69)
    overlap = disc.annulus_geometry(shape)
    baseline = replace(overlap, ring_count=(overlap.ring_count + 1) // 2,
                       ring_stride_frac=overlap.ring_width_frac)
    scale = image_half_diagonal(*shape)
    assert overlap.ring_width_frac * scale == pytest.approx(2 * scale / 32)
    assert overlap.stride_frac * scale == pytest.approx(scale / 32)
    for i in range(baseline.ring_count):
        np.testing.assert_allclose(section_bounds(overlap, 2 * i), section_bounds(baseline, i))
    np.testing.assert_array_equal(
        section_masks(shape, disc.center, overlap)[::2], section_masks(shape, disc.center, baseline)
    )
    np.testing.assert_allclose(
        np.asarray(section_bounds(overlap, 1)),
        (np.asarray(section_bounds(overlap, 0)) + np.asarray(section_bounds(overlap, 2))) / 2,
    )


@pytest.mark.parametrize("count", [1, 2, 16, 32, 65])
@pytest.mark.parametrize("center", [(50.0, 50.0), (15.0, 32.5), (84.5, 90.0)])
def test_half_diagonal_count_keeps_full_outer_annuli_to_every_image_corner(count, center):
    from calculations.topology.geometry import annulus_mask, section_masks
    from calculations.topology.mask_area import annulus_widths_pixels

    shape = (101, 101)
    disc = OpticDisc(None, center, 12, 16)
    geometry = disc.annulus_geometry(shape, number_of_radii_in_half_diagonal=count)
    scale = image_half_diagonal(*shape)
    stride = scale / count
    bounds = geometry.metadata(shape)["annulus_bounds_pixels"]
    farthest = np.hypot(max(center[0], 100 - center[0]), max(center[1], 100 - center[1]))
    assert not geometry.clip_outer_radius
    assert geometry.ring_count == int(np.ceil((farthest - 6) / stride))
    np.testing.assert_allclose(bounds[:, 0], 6 + np.arange(geometry.ring_count) * stride)
    np.testing.assert_allclose(bounds[:, 1] - bounds[:, 0], 2 * stride)
    assert bounds[-1, 0] < farthest < bounds[-1, 1]
    np.testing.assert_allclose(annulus_widths_pixels(shape, geometry, geometry.ring_count),
                               2 * stride, rtol=1e-6)
    support = section_masks(shape, center, geometry).any(axis=0)
    expected = ~disc.centered_circle_mask_for(shape)
    # Topology subtracts the exact disc mask after float32 annulus rasterization.
    np.testing.assert_array_equal(support & expected, expected)
    # Branch identification must also retain the furthest corner's support.
    branch_support = annulus_mask(shape, center, geometry.inner_radius_frac,
                                 geometry.outer_radius_frac)
    np.testing.assert_array_equal(branch_support & expected, expected)


@pytest.mark.parametrize("count", [None, 0, -1, 1.5, True, np.inf, "32"])
def test_half_diagonal_count_must_be_a_positive_integer(count):
    with pytest.raises(ValueError, match="number_of_radii_in_half_diagonal"):
        OpticDisc(None, (50, 50), 12, 16).annulus_geometry(
            (101, 101), number_of_radii_in_half_diagonal=count,
        )


@pytest.mark.parametrize("size,stride", [(12.25, None), (9.5, 2.25), (10.0, 10.0), (4.5, 7.25)])
def test_native_pixel_size_and_stride_define_every_bound(size, stride):
    from calculations.topology.mask_area import annulus_widths_pixels

    shape = (101, 151)
    disc = OpticDisc(None, (40.5, 60.0), 12, 16)
    geometry = disc.annulus_geometry(shape, annulus_size=size, annulus_stride=stride)
    actual_stride = size / 2 if stride is None else stride
    metadata = geometry.metadata(shape)
    bounds = metadata["annulus_bounds_pixels"]
    scale = image_half_diagonal(*shape)
    np.testing.assert_allclose(bounds[:, 0], 6 + np.arange(geometry.ring_count) * actual_stride)
    np.testing.assert_allclose(bounds[:, 1], np.minimum(scale, bounds[:, 0] + size))
    np.testing.assert_allclose(metadata["annulus_bounds_fraction"] * scale, bounds)
    np.testing.assert_allclose(annulus_widths_pixels(shape, geometry, geometry.ring_count),
                               np.diff(bounds, axis=1).ravel(), rtol=1e-6)
    assert metadata["annulus_size_pixels"] == pytest.approx(size)
    assert metadata["annulus_stride_pixels"] == pytest.approx(actual_stride)
    assert np.all(bounds[:, 1] > bounds[:, 0])


@pytest.mark.parametrize("extent,size,stride,count", [
    (60.0, 20.0, 10.0, 6),  # Five complete overlapping sections and a useful half-width tail.
    (65.0, 20.0, 10.0, 7),  # Clipped legacy section retained.
    (50.1, 20.0, 10.0, 5),  # Negligible EXTRA half-stride tail omitted.
    (60.1, 20.0, 10.0, 7),  # Equally thin LEGACY tail must survive.
    (60.00000000000001, 20.0, 10.0, 6),  # Floating-point boundary noise.
    (20.3, 4.0, 1.25, 17),  # Fractional-pixel stride and a meaningful partial tail.
    (20.01, 4.0, 1.25, 16),
    (60.1, 20.0, 20.0, 4),  # Historical non-overlapping clipping.
    (65.0, 10.0, 30.0, 3),  # Gapped sections.
    (0.0, 20.0, 10.0, 0),
])
def test_terminal_annulus_policy(extent, size, stride, count):
    from calculations.topology.geometry import annulus_count

    assert annulus_count(extent, size, stride) == count


def test_full_width_sampling_does_not_drop_tiny_terminal_annuli():
    from calculations.topology.geometry import annulus_count

    assert annulus_count(50.1, 20, 10, retain_all=True) == 6
    assert annulus_count(50, 20, 10, retain_all=True) == 5


def test_exact_division_has_two_n_minus_one_complete_overlapping_annuli():
    disc = OpticDisc(None, (50.0, 50.0), 4.0, 4.0)
    size = (image_half_diagonal(101, 101) - 2.0) / 3
    bounds = disc.annulus_geometry((101, 101), annulus_size=size).metadata((101, 101))["annulus_bounds_pixels"]
    widths = np.diff(bounds, axis=1).ravel()
    assert np.count_nonzero(np.isclose(widths, size)) == 2 * 3 - 1
    assert widths[-1] == pytest.approx(size / 2)


@pytest.mark.parametrize("name", ["annulus_size", "annulus_stride"])
@pytest.mark.parametrize("value", [0, -1, np.nan, np.inf, -np.inf, True, "2"])
def test_annulus_dimensions_must_be_finite_positive_native_pixels(name, value):
    with pytest.raises(ValueError, match=name):
        OpticDisc(None, (2.0, 2.0), 2, 2).annulus_geometry((10, 10), **{name: value})


def test_explicit_new_geometry_overrides_even_invalid_legacy_count():
    geometry = OpticDisc(None, (2.0, 2.0), 2, 2).annulus_geometry(
        (10, 10), number_of_radii_in_fov=0, annulus_size=2.5
    )
    assert geometry.stride_frac == pytest.approx(geometry.ring_width_frac / 2)


@pytest.mark.parametrize("pixel_pitch_mm", [0.02, 0.1, 0.2])
def test_missing_measurements_fall_back_to_centered_physical_circle(pixel_pitch_mm) -> None:
    shape = (101, 201)
    disc = OpticDisc.from_measurements(
        None, None, None, None, shape, pixel_pitch_mm=pixel_pitch_mm
    )
    radius_scale = image_half_diagonal(*shape)
    diameter = 1.91 / pixel_pitch_mm

    assert disc.is_fallback
    assert disc.center == (shape[1] / 2.0, shape[0] / 2.0)
    assert disc.width == pytest.approx(diameter)
    assert disc.height == pytest.approx(diameter)
    assert disc.annulus_geometry(shape).inner_radius_frac == pytest.approx(
        np.ceil(diameter / 2.0) / radius_scale
    )
    y, x = np.ogrid[: shape[0], : shape[1]]
    expected = (x - shape[1] / 2.0) ** 2 + (y - shape[0] / 2.0) ** 2 <= (diameter / 2.0) ** 2
    np.testing.assert_array_equal(disc.mask_for(shape), expected)


@pytest.mark.parametrize("pixel_pitch_mm", [None, 0, -0.02, np.nan, np.inf])
def test_physical_fallback_requires_valid_pixel_pitch(pixel_pitch_mm) -> None:
    with pytest.raises(ValueError, match="pixel.pitch"):
        OpticDisc.from_measurements(
            None, None, None, None, (20, 30), pixel_pitch_mm=pixel_pitch_mm
        )


def test_physical_fallback_is_clipped_without_rescaling_the_diameter() -> None:
    disc = OpticDisc.centered_fallback((3, 4), pixel_pitch_mm=0.02)
    assert disc.width == pytest.approx(95.5)
    assert disc.mask_for((3, 4)).all()


def test_nonfinite_mask_uses_full_fallback() -> None:
    shape = (21, 31)
    mask = np.full(shape, np.nan, dtype=np.float32)
    disc = OpticDisc.from_measurements(
        mask, (2.0, 3.0), 8.0, 10.0, shape, pixel_pitch_mm=0.02
    )

    assert disc.is_fallback
    assert disc.center == (15.5, 10.5)


def test_mask_supplies_missing_center_and_bounding_box_dimensions() -> None:
    mask = np.zeros((20, 30), dtype=np.float32)
    mask[4:10, 8:17] = 1.0
    disc = OpticDisc.from_measurements(
        mask,
        np.asarray([np.nan, np.nan]),
        None,
        np.asarray(np.nan),
        mask.shape,
    )

    assert not disc.is_fallback
    assert disc.center == pytest.approx((12.0, 6.5))
    assert disc.width == 9.0
    assert disc.height == 6.0
    np.testing.assert_array_equal(disc.mask, mask.astype(bool))
