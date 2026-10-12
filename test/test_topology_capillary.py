"""Contracts for the native capillary exclusion-mask calculation."""

import numpy as np
import pytest

from calculations.topology import capillary_analysis_zones


def test_complements_every_combination_without_mutating_inputs():
    combinations = np.arange(8, dtype=np.uint8).reshape(2, 4)
    masks = [(combinations & bit) != 0 for bit in (1, 2, 4)]
    originals = [mask.copy() for mask in masks]

    zones = capillary_analysis_zones(*masks)

    assert zones.dtype == np.bool_
    np.testing.assert_array_equal(zones, combinations == 0)
    for mask, original in zip(masks, originals):
        np.testing.assert_array_equal(mask, original)
        assert not np.shares_memory(zones, mask)


@pytest.mark.parametrize("invalid_index", [0, 1, 2])
def test_requires_aligned_boolean_masks(invalid_index):
    masks = [np.zeros((3, 5), dtype=bool) for _ in range(3)]
    masks[invalid_index] = np.zeros((3, 5), dtype=np.uint8)
    with pytest.raises(TypeError, match="2-D boolean"):
        capillary_analysis_zones(*masks)


def test_rejects_rank_and_shape_mismatch_without_transposing():
    mask = np.zeros((3, 5), dtype=bool)
    with pytest.raises(TypeError, match="2-D boolean"):
        capillary_analysis_zones(mask[None], mask, mask)
    with pytest.raises(ValueError, match="aligned mask shape"):
        capillary_analysis_zones(mask, mask.T, mask)


@pytest.mark.parametrize("excluded", [False, True])
def test_empty_and_full_exclusions(excluded):
    mask = np.full((4, 6), excluded, dtype=bool)
    zones = capillary_analysis_zones(mask, mask, mask)
    assert np.all(zones == (not excluded))
