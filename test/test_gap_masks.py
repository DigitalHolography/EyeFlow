"""Moment0ff gap contrast, temporal robustness, polarity, and morphology."""

import numpy as np
import pytest

from pipelines.topology_core import gap_masks


def test_temporal_median_is_tiled_nan_aware_and_does_not_mutate_input(monkeypatch):
    volume = np.arange(3 * 7 * 9, dtype=float).reshape(3, 7, 9)
    volume[:, 0, 0] = [np.nan, 4, np.inf]
    volume[:, 1, 1] = np.nan
    original = volume.copy()

    class LazyVolume:
        shape = volume.shape
        dtype = volume.dtype
        reads = []

        def __getitem__(self, key):
            self.reads.append(key)
            return volume[key]

    monkeypatch.setattr(gap_masks, "_MEDIAN_TILE_BYTES", 3 * 8 * 4)
    source = LazyVolume()
    actual = gap_masks.temporal_median_image(source)
    expected = original[1].copy()
    expected[0, 0] = 4
    np.testing.assert_array_equal(actual, expected)
    np.testing.assert_array_equal(volume, original)
    assert len(source.reads) > 1
    assert all(key[1].stop - key[1].start <= 2 for key in source.reads)


def test_dark_gaps_have_no_central_boundary_and_ignore_single_frame_outlier():
    image = np.full((64, 80), 100, dtype=np.float32)
    image[:20, :24] = 5
    volume = np.stack([image, image, image])
    volume[0, 40:50, 55:65] = 0
    mask, attrs = gap_masks.gaps_from_m0ff(volume)
    reference, _ = gap_masks.gaps_from_m0ff(np.stack([image] * 3))
    np.testing.assert_array_equal(mask, reference)
    assert mask.dtype == bool
    assert mask[0, 0]  # No field-of-view or disc exclusion.
    assert mask[10, 10]
    assert not mask[45, 60]
    assert not mask[-1, -1]
    assert attrs["gap_polarity"] == "dark"
    assert 0 < attrs["otsu_threshold"] < 1


@pytest.mark.parametrize("vessel_yx", [(10, 10), (0, 0)])
def test_dilates_bright_vessels_with_7x7_square_before_inverting(monkeypatch, vessel_yx):
    monkeypatch.setattr(gap_masks, "equalize_adapthist", lambda image, **kwargs: image)
    monkeypatch.setattr(gap_masks.ndi, "gaussian_filter", lambda image, sigma: image)
    monkeypatch.setattr(gap_masks, "threshold_otsu", lambda image: 0.5)
    image = np.zeros((21, 21))
    image[vessel_yx] = 1
    mask, attrs = gap_masks.gaps_from_m0ff(image[None])
    y, x = np.indices(image.shape)
    expected = (np.abs(y - vessel_yx[0]) > 3) | (np.abs(x - vessel_yx[1]) > 3)
    np.testing.assert_array_equal(mask, expected)
    assert attrs["dilation_footprint"] == "7x7_square"
    assert attrs["dilation_applies_to"] == "bright_vessels"


@pytest.mark.parametrize("value", [0, 8, np.nan, np.inf])
def test_uniform_or_invalid_image_has_no_distinguishable_gaps(value):
    mask, _ = gap_masks.gaps_from_m0ff(np.full((3, 24, 32), value))
    assert not mask.any()


def test_invalid_pixels_remain_background_after_dilation():
    image = np.ones((32, 40))
    image[:16, :20] = 0
    image[5, 5] = np.nan
    mask, _ = gap_masks.gaps_from_m0ff(image[None])
    assert mask[5, 6]
    assert not mask[5, 5]


@pytest.mark.parametrize("dtype", [complex, str, bool])
def test_non_real_numeric_source_is_rejected(dtype):
    with pytest.raises(TypeError, match="real numeric"):
        gap_masks.gaps_from_m0ff(np.ones((1, 4, 4)).astype(dtype))
