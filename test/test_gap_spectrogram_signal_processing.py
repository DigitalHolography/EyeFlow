"""Physical temporal FFT cutoffs and centered SVD across concatenated beats."""

import numpy as np
import pytest

from calculations.math.low_rank import remove_leading_svd_modes
from pipelines.velocity.gaps.signal_processing import (
    band_limit_gap_spectrograms,
    svd_filter_gap_spectrograms,
)


def test_fft_cutoff_uses_original_beat_duration_and_preserves_dc():
    phase = np.arange(128) / 128
    # Harmonic 10 is 20 Hz for a half-second beat, 10 Hz for a one-second beat.
    slow = np.sin(2 * np.pi * 3 * phase)
    tenth = np.cos(2 * np.pi * 10 * phase)
    values = np.broadcast_to((20 + slow + tenth)[:, None, None], (128, 2, 2)).copy()
    result = band_limit_gap_spectrograms(values, np.array([0.5, 1.0]))
    np.testing.assert_allclose(result.filtered[:, 0, 0], 20 + slow, atol=3e-6)
    np.testing.assert_allclose(result.filtered[:, 0, 1], values[:, 0, 1], atol=3e-6)
    np.testing.assert_allclose(result.filtered.mean(axis=0), 20, atol=3e-6)
    np.testing.assert_array_equal(result.highest_harmonic, [7, 15])
    assert result.filtered.dtype == np.float32


def test_fft_keeps_inclusive_15_hz_boundary_and_dc_only_when_cutoff_is_zero():
    phase = np.arange(128) / 128
    values = (5 + np.cos(2 * np.pi * 15 * phase) + np.cos(2 * np.pi * 16 * phase))[:, None, None]
    result = band_limit_gap_spectrograms(values, np.array([1.0]))
    np.testing.assert_allclose(result.filtered[:, 0, 0], 5 + np.cos(2 * np.pi * 15 * phase), atol=1e-6)
    np.testing.assert_allclose(band_limit_gap_spectrograms(values, np.array([1.0]), cutoff_hz=0).filtered, 5, atol=1e-6)


def test_fft_incomplete_series_stays_unchanged_without_affecting_other_bins_or_beats():
    values = np.arange(128 * 2 * 2, dtype=np.float32).reshape(128, 2, 2)
    values[8, 0, 0] = np.nan
    original = values.copy()
    result = band_limit_gap_spectrograms(values, np.ones(2))
    np.testing.assert_array_equal(result.filtered[:, 0, 0], values[:, 0, 0])
    assert np.isfinite(result.filtered[:, 1, 0]).all()
    np.testing.assert_array_equal(result.filtered_series, [[False, True], [True, True]])
    np.testing.assert_array_equal(values, original)


def _three_mode_matrix():
    sample = np.arange(128)
    # Mutually orthogonal centered temporal modes and distinct singular values.
    matrix = np.column_stack([
        20 + 8 * np.sin(2 * np.pi * sample / 128),
        30 + 4 * np.cos(2 * np.pi * sample / 128),
        40 + np.sin(4 * np.pi * sample / 128),
    ])
    expected = matrix.copy()
    expected[:, :2] = [20, 30]
    return matrix, expected


def test_svd_removes_two_leading_modes_restores_baseline_and_retains_chronological_continuous_output():
    matrix, expected = _three_mode_matrix()
    values = matrix.reshape(2, 64, 3).transpose(1, 2, 0)
    original = values.copy()
    filtered, result = svd_filter_gap_spectrograms(values)
    assert filtered.shape == (128, 3)
    np.testing.assert_allclose(filtered, expected, atol=2e-6)
    np.testing.assert_allclose(filtered.mean(axis=0), [20, 30, 40], atol=3e-6)
    assert result.modes_removed == 2
    assert result.applied
    assert filtered.dtype == np.float32
    np.testing.assert_array_equal(values, original)


def test_svd_preserves_invalid_rows_excludes_empty_bins_and_never_fills_nan():
    matrix, _ = _three_mode_matrix()
    values = np.column_stack([matrix, np.full(128, np.nan)])
    values[10, 0] = np.nan
    original = values.copy()
    result = remove_leading_svd_modes(values)
    assert result.applied
    np.testing.assert_array_equal(result.filtered[10], values[10].astype(np.float32))
    assert np.isnan(result.filtered[:, 3]).all()
    assert not result.valid_rows[10]
    np.testing.assert_array_equal(result.valid_columns, [True, True, True, False])
    np.testing.assert_allclose(result.filtered[result.valid_rows, :3].mean(axis=0), values[result.valid_rows, :3].mean(axis=0), atol=5e-6)
    np.testing.assert_array_equal(values, original)


def test_insufficient_svd_dimensions_return_unchanged_instead_of_erasing_all_modes():
    for values in (np.ones((128, 2)), np.full((128, 3), np.nan), np.ones((2, 3))):
        result = remove_leading_svd_modes(values)
        assert not result.applied
        assert result.modes_removed == 0
        np.testing.assert_array_equal(result.filtered, values)


@pytest.mark.parametrize("durations", [[0], [-1], [np.nan], [1, 2]])
def test_invalid_beat_durations_fail_explicitly(durations):
    with pytest.raises(ValueError, match="duration"):
        band_limit_gap_spectrograms(np.ones((128, 3, 1)), durations)


@pytest.mark.parametrize("mode_count", [-1, 1.5, True])
def test_invalid_svd_mode_count_fails_explicitly(mode_count):
    with pytest.raises(ValueError, match="mode_count"):
        remove_leading_svd_modes(np.ones((128, 3)), mode_count=mode_count)
