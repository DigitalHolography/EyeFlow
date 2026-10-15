"""Contracts for the native capillary exclusion-mask calculation."""

import numpy as np
import pytest

from calculations.topology import capillary_analysis_zones, downsample_gap_zones


def test_complements_every_combination_without_mutating_inputs():
    centers = [(y, x) for y in (24, 32, 40) for x in (24, 32, 40)][:8]
    masks = [np.zeros((64, 64), dtype=bool) for _ in range(3)]
    for combination, (y, x) in enumerate(centers):
        for mask, bit in zip(masks, (1, 2, 4)):
            mask[y, x] = bool(combination & bit)
    originals = [mask.copy() for mask in masks]

    zones = capillary_analysis_zones(*masks)

    assert zones.dtype == np.bool_
    for combination, (y, x) in enumerate(centers):
        assert zones[y, x] == (combination == 0)
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
    mask = np.full((20, 20), excluded, dtype=bool)
    zones = capillary_analysis_zones(mask, mask, mask)
    if excluded:
        assert not zones.any()
    else:
        assert zones[10, 10]
        assert not zones[0, 0]


@pytest.mark.parametrize("size", [20, 21])
def test_cutoff_is_a_centered_circle_with_eighty_percent_width_diameter(size):
    empty = np.zeros((size, size), dtype=bool)
    zones = capillary_analysis_zones(empty, empty, empty)
    center = size / 2.0
    radius = size * 0.8 / 2.0
    expected = np.array([
        [(x - center) ** 2 + (y - center) ** 2 <= radius**2 for x in range(size)]
        for y in range(size)
    ])
    np.testing.assert_array_equal(zones, expected)
    if size == 20:
        assert zones[10, 2] and zones[10, 18]  # The circle includes its boundary.
        assert not zones[10, 1] and not zones[10, 19]
        assert not zones[4, 4] and zones[5, 5]


@pytest.mark.parametrize("vessel_index", [0, 2])
def test_dilates_both_vessel_classes_with_three_pixel_square_including_diagonals(vessel_index):
    masks = [np.zeros((64, 64), dtype=bool) for _ in range(3)]
    masks[vessel_index][32, 32] = True
    empty_zones = capillary_analysis_zones(*[np.zeros((64, 64), bool) for _ in range(3)])
    expected = empty_zones.copy()
    expected[29:36, 29:36] = False

    zones = capillary_analysis_zones(*masks)

    np.testing.assert_array_equal(zones, expected)
    assert not zones[29, 32] and not zones[35, 32]
    assert not zones[30, 31] and not zones[34, 33]
    assert not zones[29, 29] and not zones[35, 35]
    assert not zones[30, 30] and not zones[34, 34]
    assert zones[28, 32] and zones[36, 32]
    assert zones[32, 28] and zones[32, 36]


def test_does_not_dilate_the_optic_disc_or_erode_the_cutoff():
    empty = np.zeros((20, 20), dtype=bool)
    optic_disc = empty.copy()
    optic_disc[10, 10] = True
    zones = capillary_analysis_zones(empty, optic_disc, empty)
    assert not zones[10, 10]
    assert zones[10, 9] and zones[9, 10]
    assert zones[10, 2] and zones[10, 18]


def test_dilates_original_vessels_before_clipping_to_the_cutoff():
    empty = np.zeros((64, 64), dtype=bool)
    vessel = empty.copy()
    vessel[32, 5] = True  # Outside the circle, close enough to exclude interior pixels.
    zones = capillary_analysis_zones(vessel, empty, empty)
    assert not zones[32, 7] and not zones[32, 8]
    assert zones[32, 9]


def test_dilation_clips_at_image_borders_without_wrapping():
    empty = np.zeros((8, 8), dtype=bool)
    vessel = empty.copy()
    vessel[0, 4] = True
    zones = capillary_analysis_zones(vessel, empty, empty)
    assert not zones[3, 4]
    assert zones[4, 4]
    assert zones[7, 4]


def test_rejects_empty_spatial_dimensions():
    empty = np.zeros((0, 5), dtype=bool)
    with pytest.raises(ValueError, match="positive image dimensions"):
        capillary_analysis_zones(empty, empty, empty)


def test_interpolation_samples_pixel_centers_without_all_true_block_exclusion():
    source = np.zeros((6, 8), dtype=bool)
    expected = np.array([
        [True, False, False, True],
        [False, True, False, False],
        [True, False, True, False],
    ])
    source[1::2, 1::2] = expected
    original = source.copy()
    # The other three pixels of each 2x2 block are false; selected pixels survive.
    reduced = downsample_gap_zones(source, target_shape=(3, 4))
    np.testing.assert_array_equal(reduced, expected)
    np.testing.assert_array_equal(source, original)
    assert not np.shares_memory(reduced, source)
    assert reduced.dtype == np.bool_


@pytest.mark.parametrize("dtype", [np.float32, np.uint8])
def test_downsampling_requires_boolean_support(dtype):
    with pytest.raises(TypeError, match="2-D boolean"):
        downsample_gap_zones(np.ones((4, 4), dtype=dtype), target_shape=(2, 2))


@pytest.mark.parametrize("shape", [(0, 4), (4, 0)])
def test_interpolation_rejects_empty_source_dimensions(shape):
    with pytest.raises(ValueError, match="positive source spatial dimensions"):
        downsample_gap_zones(np.zeros(shape, dtype=bool), target_shape=(2, 3))


@pytest.mark.parametrize("target_shape", [(128, 128), (256, 256), (170, 129), (64, 64), (700, 300)])
def test_velocity_application_interpolates_to_video_dimensions_without_changing_export(target_shape):
    from pipelines.low_freq_spectra.masks import gap_velocity_mask

    original = np.ones((512, 512), dtype=bool)
    original[100, 100] = False
    original[101, 103] = False
    original[102, 102] = False
    preserved = original.copy()
    reduced = gap_velocity_mask(original, target_shape)
    # Center coordinates select floor((output_index + 0.5) * input_size / output_size).
    rows = np.floor((np.arange(target_shape[0]) + .5) * 512 / target_shape[0]).astype(int)
    columns = np.floor((np.arange(target_shape[1]) + .5) * 512 / target_shape[1]).astype(int)
    expected = original[np.ix_(rows, columns)]
    assert reduced.shape == target_shape
    np.testing.assert_array_equal(reduced, expected)
    np.testing.assert_array_equal(original, preserved)
    assert original.shape == (512, 512)
    assert not np.shares_memory(reduced, original)


def test_velocity_application_preserves_matching_grid():
    from pipelines.low_freq_spectra.masks import gap_velocity_mask

    source = np.ones((256, 256), dtype=bool)
    same_grid = gap_velocity_mask(source, (256, 256))
    np.testing.assert_array_equal(same_grid, source)
    assert not np.shares_memory(same_grid, source)


def test_noninteger_resize_preserves_native_yx_axes_and_boolean_values():
    source = np.zeros((5, 7), dtype=bool)
    source[0, 6] = source[2, 2] = source[4, 0] = True
    original = source.copy()
    expected = np.array([
        [False, False, False, True],
        [False, True, False, False],
        [True, False, False, False],
    ])
    reduced = downsample_gap_zones(source, target_shape=(3, 4))
    np.testing.assert_array_equal(reduced, expected)
    np.testing.assert_array_equal(source, original)
    assert not np.shares_memory(reduced, source)


@pytest.mark.parametrize("target_shape", [(0, 2), (-1, 2), (2, 0), (1.5, 2), (True, 2), (np.bool_(True), 2), (2,), (2, 3, 4)])
def test_interpolation_rejects_invalid_target_dimensions(target_shape):
    with pytest.raises(ValueError, match="target_shape"):
        downsample_gap_zones(np.ones((8, 8), dtype=bool), target_shape=target_shape)
