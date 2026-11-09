"""Small direct spatial references for the frequency-resolved Doppler CVI proxy."""

import warnings

import numpy as np
import pytest
from skimage.filters import threshold_otsu

from pipelines.low_freq_spectra.cvi import (
    calculate_cvi, disc_mask, calculate_binary_occupancy_maps, mean_frequency_occupancy,
)
from pipelines.low_freq_spectra.per_beat import gap_spectrograms_per_beat, align_power
from pipelines.low_freq_spectra.spectrograms import preprocess_spectrograms


def direct(source, boundaries, phases, bins=256):
    """Materialize logical S only for tiny synthetic test inputs."""
    nt, nf, h, w = source.shape
    support = disc_mask(h, w)
    s = np.full((phases, len(boundaries) - 1, nf, h, w), np.nan)
    for b, (start, stop) in enumerate(zip(boundaries[:-1], boundaries[1:])):
        for f in range(nf):
            for y in range(h):
                for x in range(w):
                    values = source[start:stop + 1, f, y, x].copy()
                    values[~np.isfinite(values)] = np.nan
                    s[:, b, f, y, x] = np.interp(
                        np.arange(phases) / phases, np.linspace(0, 1, stop - start + 1), values)
    valid = np.isfinite(s[..., support]).all(axis=-1)
    masked = s.copy()
    masked[~valid] = np.nan
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        a = np.nanmedian(masked, axis=1)
    thresholds = np.full(nf, np.nan)
    for f in range(nf):
        reference = a[:, f, support]
        reference = reference[np.isfinite(reference)]
        if reference.size:
            thresholds[f] = threshold_otsu(reference, nbins=bins)
    fractions = (s[..., support] > thresholds[None, None, :, None]).sum(axis=-1) / support.sum()
    fractions[~valid] = np.nan
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        median = np.nanmedian(fractions, axis=1)
    return thresholds, fractions, median, valid.sum(axis=1)


@pytest.mark.parametrize("height,width", [(7, 7), (3, 10), (10, 3), (1, 5)])
def test_disc_geometry(height, width):
    expected = np.array([[(x - (width - 1)/2)**2 + (y - (height - 1)/2)**2 <= (0.4*width)**2
                          for x in range(width)] for y in range(height)])
    np.testing.assert_array_equal(disc_mask(height, width), expected)
    if (height, width) == (3, 10):
        assert expected.sum() > np.pi * (height / 2)**2  # radius is not fitted to height


def test_fractions_are_beat_medians_not_fraction_of_median_image():
    source = np.zeros((7, 2, 1, 5))
    for b, pixels in enumerate(([0, 1], [2, 3], [0, 4])):
        source[2*b:2*b + 2, 0, 0, pixels] = 10
    source[:, 1] = -3  # signed bins are independent, finite negative power retained
    result = calculate_cvi(source, [0, 2, 4, 6], [1.5, -2.5], sample_count=2)
    np.testing.assert_allclose(result["per_beat"][0, :, 0], 0.4)
    assert result["cvi"][0, 0] == pytest.approx(0.4)
    assert np.count_nonzero(np.median(source[[0, 2, 4], 0], axis=0) > result["thresholds"][0]) / 5 == 0.2
    assert result["thresholds"][1] == -3
    assert result["flags"][1] == 1
    assert np.all(result["cvi"][:, 1] == 0)


def test_interpolation_precedes_binarization_and_matches_existing_phase_conventions():
    source = np.array([0., 10., 0.])[:, None, None, None]
    result = calculate_cvi(source, [0, 2], [1.], sample_count=8)
    power = np.interp(np.arange(8)/8, [0, .5, 1], [0, 10, 0])
    threshold = threshold_otsu(power, nbins=256)
    assert result["thresholds"][0] == threshold
    np.testing.assert_array_equal(result["cvi"][:, 0], power > threshold)
    interpolated_binary = np.interp(np.arange(8)/8, [0, .5, 1], [0, 1, 0])
    assert not np.array_equal(result["cvi"][:, 0], interpolated_binary)
    values = np.arange(42, dtype=float).reshape(7, 6)
    values[2, 1] = np.inf
    actual = align_power(values, np.array([0, 3, 6]), 128)
    expected = gap_spectrograms_per_beat(values, [0, 3, 6])
    np.testing.assert_allclose(actual.transpose(0, 2, 1), expected, rtol=1e-7, equal_nan=True)


@pytest.mark.parametrize("memory", [0.001, 64])
def test_tiled_calculation_matches_direct_reference_and_preserves_source(memory):
    source = np.random.default_rng(182).normal(2, 4, (9, 4, 4, 7))
    source[2, 0, 1, 3] = np.nan
    source[4, 1, 2, 2] = np.inf
    source[:, 2] = np.nan
    original = source.copy()
    result = calculate_cvi(source, [1, 4, 7], [2, 3, -3, -2], sample_count=8, working_memory_mb=memory)
    expected = direct(source, [1, 4, 7], 8)
    for key, reference in zip(("thresholds", "per_beat", "cvi", "valid_counts"), expected):
        np.testing.assert_allclose(result[key], reference, equal_nan=True)
    assert result["flags"][2] == 2
    np.testing.assert_array_equal(source, original)
    assert result["pixel_count"] == disc_mask(4, 7).sum()
    # Spatial in-place medians operate on owned buffers, leaving the CVI source intact.
    reduced = preprocess_spectrograms(source, None, workers=1)
    np.testing.assert_array_equal(source, original)
    from input_output.schema import EyeFlowOutputPaths
    actual = reduced[EyeFlowOutputPaths.active().spectrograms_spatial_median][0]
    median = np.median(source, axis=(-1, -2))
    median[~np.isfinite(source).all(axis=(-1, -2))] = np.nan
    np.testing.assert_allclose(actual, np.maximum(median, 0), equal_nan=True)


def test_invalid_images_excluded_wholly_from_reference_and_fixed_denominator():
    source = np.zeros((5, 1, 1, 5))
    source[0, 0, 0] = [np.nan, 1000, 1000, 1000, 1000]
    source[2, 0, 0] = [0, 0, 10, 10, 10]
    result = calculate_cvi(source, [0, 2, 4], [1], sample_count=2)
    assert np.isnan(result["per_beat"][0, 0, 0])
    assert result["valid_counts"][0, 0] == 1
    assert result["cvi"][0, 0] == .6
    assert result["thresholds"][0] < 10  # invalid image's 1000s cannot affect threshold
    source[:] = np.nan
    result = calculate_cvi(source, [0, 2, 4], [1], sample_count=2)
    assert np.all(np.isnan(result["cvi"]))
    assert result["flags"][0] == 2
    assert np.all(result["valid_counts"] == 0)


def test_nonfinite_outside_disc_is_ignored_and_equality_is_black():
    source = np.full((3, 1, 5, 5), -2.)
    source[:, :, 0, 0] = np.inf
    result = calculate_cvi(source, [0, 2], [-1], sample_count=4)
    assert result["flags"][0] == 1
    assert result["thresholds"][0] == -2
    assert np.all(result["cvi"] == 0)
    assert np.all(result["valid_counts"] == 1)


@pytest.mark.parametrize("boundaries", [[], [0]])
def test_insufficient_complete_beats_are_rejected_without_inventing_timing(boundaries):
    with pytest.raises(ValueError, match="At least two"):
        calculate_cvi(np.zeros((3, 1, 1, 5)), boundaries, [1])


@pytest.mark.parametrize("parameters", [{"histogram_bins": 1}, {"working_memory_mb": 0}])
def test_invalid_calculation_arguments(parameters):
    with pytest.raises(ValueError):
        calculate_cvi(np.ones((3, 1, 1, 5)), [0, 2], [1], **parameters)


def test_empty_geometric_intersection_stays_invalid():
    result = calculate_cvi(np.ones((3, 1, 2, 1)), [0, 2], [1.], sample_count=2)
    assert result["pixel_count"] == 0
    assert result["flags"][0] == 2
    assert np.isnan(result["thresholds"][0])
    assert np.isnan(result["cvi"]).all()
    assert not result["valid_counts"].any()


@pytest.mark.parametrize("memory", [.001, 64])
def test_binary_averages_match_direct_reference_with_equal_phase_frequency_weights(memory):
    source = np.random.default_rng(28).normal(size=(7, 5, 3, 4))
    thresholds = np.array([-.5, 0, .5, 1, np.nan])
    valid = np.ones((4, 2, 5), dtype=bool)
    valid[0, 1, 1] = False
    valid[1, :, 2] = False
    valid[2, 0, 3] = False
    valid[3, :, 1:4] = False
    original = source.copy()
    actual = calculate_binary_occupancy_maps(source, [0, 3, 6], thresholds, valid,
                                             peak_indices=(1, 3), working_memory_mb=memory)
    s = np.stack([align_power(source[:, f].reshape(7, -1), np.array([0, 3, 6]), 4)
                  for f in range(5)], axis=2).reshape(4, 2, 5, 3, 4)
    valid &= np.isfinite(thresholds)[None, None, :]
    binary = np.where(valid[..., None, None], s > thresholds[None, None, :, None, None], np.nan)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        beat_mean = np.nanmean(binary, axis=1)
        phase_video = np.nanmean(beat_mean[:, 1:4], axis=1)
        frequency_video = np.nanmean(beat_mean, axis=0)
    support = disc_mask(3, 4)
    phase_video[:, ~support] = np.nan
    frequency_video[:, ~support] = np.nan
    np.testing.assert_allclose(actual["phase_images"], phase_video, rtol=1e-7, equal_nan=True)
    np.testing.assert_allclose(actual["frequency_images"], frequency_video, rtol=1e-7, equal_nan=True)
    np.testing.assert_array_equal(actual["selected_bins"], [1, 2, 3])
    np.testing.assert_array_equal(actual["phase_bin_counts"], valid.any(axis=1)[:, 1:4].sum(axis=1))
    np.testing.assert_array_equal(actual["frequency_phase_counts"], valid.any(axis=1).sum(axis=0))
    assert np.isnan(actual["phase_images"][3]).all()
    assert np.isnan(actual["frequency_images"][4]).all()
    assert np.any((actual["phase_images"] > 0) & (actual["phase_images"] < 1))
    np.testing.assert_array_equal(source, original)


def test_averaging_binary_values_does_not_rethreshold_or_binarize_mean_power():
    source = np.zeros((3, 2, 1, 5))
    source[:, 1] = 100
    result = calculate_binary_occupancy_maps(source, [0, 2], [1, 1], np.ones((2, 1, 2), bool),
                                            peak_indices=(0, 1))
    np.testing.assert_array_equal(result["phase_images"], .5)
    np.testing.assert_array_equal(result["frequency_images"][0], 0)
    np.testing.assert_array_equal(result["frequency_images"][1], 1)
    missing = calculate_binary_occupancy_maps(source, [0, 2], [1, 1], np.zeros((2, 1, 2), bool))
    assert missing["phase_images"] is None
    assert np.isnan(missing["frequency_images"]).all()


def test_mean_frequency_map_uses_equal_valid_bin_weights_and_bounded_reads():
    frames = np.array([[[0., .25, np.nan], [1., np.nan, np.nan]],
                       [[1., .75, np.nan], [0., .2, np.nan]],
                       [[0., np.nan, np.nan], [.5, .8, np.nan]]])
    original = frames.copy()

    class FrameSource:
        shape = frames.shape
        indexes = []

        def __getitem__(self, index):
            assert isinstance(index, int)
            self.indexes.append(index)
            return frames[index]

    source = FrameSource()
    actual = mean_frequency_occupancy(source)
    expected = [[1/3, .5, np.nan], [.5, .5, np.nan]]
    np.testing.assert_allclose(actual, expected, equal_nan=True)
    assert source.indexes == [0, 1, 2]
    np.testing.assert_array_equal(frames, original)
