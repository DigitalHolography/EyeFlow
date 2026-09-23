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
    expected_inner = 3.0 / scale

    assert isinstance(geometry, AnnulusGeometry)
    assert geometry.inner_radius_frac == pytest.approx(expected_inner)
    assert geometry.outer_radius_frac == 1.0
    assert geometry.ring_width_frac == pytest.approx(expected_step)
    assert geometry.segment_length_frac == pytest.approx(expected_step)
    assert geometry.ring_count == int(np.ceil((1.0 - expected_inner) / expected_step))
    with pytest.raises(ValueError, match="positive"):
        disc.annulus_geometry((10, 20), number_of_radii_in_fov=0)
    assert OpticDisc(None, (0.0, 0.0), 1.0, 1.0).annulus_geometry(
        (1, 1)
    ).inner_radius_frac == 0.5
