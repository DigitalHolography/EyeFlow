"""Spectral artifact axes, native mask orientation, and shared output routing."""

import numpy as np
import pytest
from matplotlib.colors import Normalize, PowerNorm

from pipelines.low_freq_spectra.artifacts import (
    FIGURE_ASPECT_RATIO,
    FIGURE_WIDTH_INCHES,
    SPECTROGRAM_ASPECT_RATIO,
    SPECTROGRAM_WIDTH_INCHES,
    _gap_spectrogram_figure,
    _spectrum_curve,
)
from calculations.math.frequency import frequency_coordinates_khz


def test_heatmap_preserves_frame_bin_values_and_nan_gaps():
    frames, dt_seconds = 640, 0.025
    time_ticks, monochrome = np.arange(0, 16, 5), True
    values = np.arange(frames * 32, dtype=float).reshape(frames, 32)
    values[51, 2] = np.nan
    original = values.copy()
    fig = _gap_spectrogram_figure(
        values, frequencies_khz=frequency_coordinates_khz(32, 64000),
        dt_seconds=dt_seconds, monochrome=monochrome,
    )
    ax, colorbar = fig.axes
    image = ax.images[0]

    np.testing.assert_array_equal(image.get_array().filled(np.nan), values.T)
    np.testing.assert_array_equal(values, original)
    np.testing.assert_allclose(ax.get_xticks(), time_ticks)
    np.testing.assert_allclose(image.get_extent(), [
        -0.5 * dt_seconds, (frames - 0.5) * dt_seconds, -0.5, 31.5,
    ])
    # Every pixel's center is its zero-based frame index times the HD interval.
    left, right = image.get_extent()[:2]
    np.testing.assert_allclose(
        left + (np.arange(frames) + 0.5) * (right - left) / frames,
        np.arange(frames) * dt_seconds, atol=1e-14,
    )
    np.testing.assert_array_equal(ax.get_yticks(), [0, 4, 8, 12, 16, 20, 24, 28, 31])
    assert [label.get_text() for label in ax.get_yticklabels()] == ["0", "8", "16", "24", "-32", "-24", "-16", "-8", "-2"]
    assert ax.get_xlabel() == "Time (s)"
    assert ax.get_ylabel() == "Frequency (kHz)"
    assert image.origin == "lower"
    assert image.get_cmap().name == ("gray_r" if monochrome else "viridis")
    if monochrome:
        np.testing.assert_array_equal(image.cmap([0, 1.0])[:, :3], [[1, 1, 1], [0, 0, 0]])
        quarter_power = np.nanmin(values) + 0.25 * np.ptp(values[np.isfinite(values)])
        assert 0.25 < image.norm(quarter_power) < 0.3
        np.testing.assert_allclose(image.norm.inverse(image.norm(quarter_power)), quarter_power)
    assert type(image.norm) is (PowerNorm if monochrome else Normalize)
    assert (image.norm.vmin, image.norm.vmax) == (np.nanmin(values), np.nanmax(values))
    assert colorbar.get_ylabel() == "Power (a.u.)"
    assert image.colorbar.orientation == "vertical"
    np.testing.assert_allclose(fig.get_size_inches(), [
        SPECTROGRAM_WIDTH_INCHES, SPECTROGRAM_WIDTH_INCHES / SPECTROGRAM_ASPECT_RATIO,
    ])
    assert fig.get_figheight() < FIGURE_WIDTH_INCHES / FIGURE_ASPECT_RATIO
    for axis in fig.axes:
        assert all(spine.get_linewidth() == 0.4 for spine in axis.spines.values())
        assert axis.yaxis.label.get_fontsize() == 10
        assert axis.yaxis.label.get_fontfamily() == ["Times New Roman"]
        assert all(label.get_fontsize() == 8 for label in axis.get_xticklabels())
        assert all(label.get_rotation() == 0 for label in axis.get_xticklabels())
        assert all(label.get_fontsize() == 8 for label in axis.get_yticklabels())
        assert all(label.get_rotation() == 0 for label in axis.get_yticklabels())
    fig.set_dpi(320)
    fig.canvas.draw()
    bounds = [label.get_window_extent(fig.canvas.get_renderer()) for label in ax.get_xticklabels()]
    assert all(left.x1 < right.x0 for left, right in zip(bounds, bounds[1:]))
    bounds = [label.get_window_extent(fig.canvas.get_renderer()) for label in ax.get_yticklabels()]
    assert all(bottom.y1 < top.y0 for bottom, top in zip(bounds, bounds[1:]))
    for axis in fig.axes:
        for label in (axis.xaxis.label, axis.yaxis.label):
            if label.get_text():
                box = label.get_window_extent(fig.canvas.get_renderer())
                assert 0 <= box.x0 < box.x1 <= fig.bbox.width
                assert 0 <= box.y0 < box.y1 <= fig.bbox.height


def test_spectra_keep_both_frequency_halves_and_share_stored_peaks_and_hv_window():
    first, second = np.ones(64), np.ones(64)
    first[5], second[7] = 10, 20
    first[40], second[50] = 10000, 20000  # Outside the Hv neighborhood.
    first[12] = np.nan
    first[24:40] = 4
    second[24:40] = 6
    first[29:32] = [2, 4, 100]
    second[29:32] = [3, 6, 100]
    original = first.copy()
    frequencies = frequency_coordinates_khz(64, 64000)
    a = _spectrum_curve(first, np.full(64, 2), frequencies_khz=frequencies, label="Early", linestyle="--", index_max1=5, index_max2=55)
    b = _spectrum_curve(second, np.full(64, 5), frequencies_khz=frequencies, label="Late", linestyle=":", index_max1=5, index_max2=55)
    np.testing.assert_array_equal(a.frequencies_khz, frequencies)
    np.testing.assert_array_equal(b.frequencies_khz, frequencies)
    assert a.peak_bins == b.peak_bins == (5, 55)
    np.testing.assert_array_equal(a.valley_bins, np.arange(24, 40))
    np.testing.assert_allclose(a.power, first / 4)
    np.testing.assert_allclose(b.power, second / 6)
    assert a.power[5] == 10 / 4
    assert b.power[7] == 20 / 6
    assert b.power[0] == 1 / 6
    np.testing.assert_array_equal(a.std, [2 / 4] * 64)
    np.testing.assert_array_equal(b.std, [5 / 6] * 64)
    assert np.isnan(a.power[12])
    np.testing.assert_equal(first, original)
