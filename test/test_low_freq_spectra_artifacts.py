"""Spectral artifact axes, native mask orientation, and shared output routing."""

from io import BytesIO
from types import SimpleNamespace
from unittest.mock import patch

import h5py
import numpy as np
import pytest
from matplotlib.colors import Normalize, PowerNorm
from PIL import Image

from input_output.writers.h5_output import PipelineH5Output
from input_output.runs.layout import HoloRunLayout
from input_output.output_manager import OutputManager
from input_output.schema import EyeFlowOutputPaths
from pipelines.low_freq_spectra.artifacts import (
    FIGURE_ASPECT_RATIO,
    FIGURE_WIDTH_INCHES,
    SPECTROGRAM_ASPECT_RATIO,
    SPECTROGRAM_WIDTH_INCHES,
    SPECTRUM_FONT_SIZE,
    SPECTRUM_HEIGHT_SCALE,
    SPECTRUM_TICK_FONT_SIZE,
    _gap_median_figure,
    _gap_metric_figure,
    _gap_spectrogram_figure,
    _gaps_mask_figure,
    _spectrum_curve,
    _spectrum_baseline,
    _spectra_figure,
    export_spectrogram_figures,
)
from calculations.math.frequency import frequency_coordinates_khz
from test_spectrogram_frequency import _write_frequency_table


@pytest.mark.parametrize("frames, dt_seconds, time_ticks, monochrome", [
    (103, 0.01, np.arange(0, 1.01, 0.2), False),
    (640, 0.025, np.arange(0, 16, 5), True),
])
def test_heatmap_preserves_frame_bin_values_and_nan_gaps(frames, dt_seconds, time_ticks, monochrome):
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


def test_median_labels_signed_khz_without_reordering_or_filling_nan():
    values = np.array([20, 5, np.nan, 10])
    ax = _gap_median_figure(values, frequencies_khz=frequency_coordinates_khz(4, 8000)).axes[0]
    np.testing.assert_array_equal(ax.lines[0].get_xdata(), [0, 1, 2, 3])
    np.testing.assert_array_equal(ax.lines[0].get_ydata(), values)
    zero_line = ax.lines[1]
    np.testing.assert_array_equal(zero_line.get_ydata(), [0, 0])
    assert zero_line.get_linestyle() == ":"
    assert zero_line.get_color() == "black"
    assert zero_line.get_linewidth() == 0.4
    np.testing.assert_array_equal(ax.get_xticks(), [0, 1, 2, 3])
    assert [label.get_text() for label in ax.get_xticklabels()] == ["0", "2", "-4", "-2"]
    assert all(label.get_fontsize() == 8 for label in ax.get_xticklabels())
    assert all(label.get_rotation() == 90 for label in ax.get_xticklabels())
    assert all(label.get_fontsize() == 8 for label in ax.get_yticklabels())
    assert all(label.get_rotation() == 0 for label in ax.get_yticklabels())
    assert ax.get_xlabel() == "Frequency (kHz)"
    assert ax.get_ylabel() == "Power (a.u.)"


def test_mask_preserves_published_xy_orientation_and_square_pixels():
    mask = np.array([[0, 1, 1], [0, 0, 0], [1, 0, 0]], dtype=np.uint8)
    fig = _gaps_mask_figure(mask)
    ax = fig.axes[0]
    image = ax.images[0]
    np.testing.assert_array_equal(image.get_array(), mask.T.astype(bool))
    assert image.origin == "lower"
    assert ax.get_aspect() == 1
    np.testing.assert_array_equal(image.cmap(image.norm([0, 1]))[:, :3],
                                  [[0, 0, 0], [1, 1, 1]])
    assert ax.get_xticks().size == ax.get_yticks().size == 0
    assert ax.get_xlabel() == ax.get_ylabel() == ""
    assert not ax.axison
    np.testing.assert_allclose(ax.get_position().bounds, [0, 0, 1, 1])
    np.testing.assert_allclose(fig.get_size_inches(), [FIGURE_WIDTH_INCHES] * 2)


@pytest.mark.parametrize("metric, unit, label", [
    ("A0", "power_unit", "$A_{0}$ (power_unit)"), ("SE_50", "1", "$SE_{50}$"),
    ("SE", "1", "$SE$"),
    ("A0", "a.u.*kHz", "$A_{0}$ (a.u.)"),
    ("A", "a.u.*kHz", "$A$ (a.u.)"),
    ("A_25", "a.u.*kHz", "$A_{25}$ (a.u.)"),
    ("A_50", "a.u.*kHz", "$A_{50}$ (a.u.)"),
    ("A_75", "a.u.*kHz", "$A_{75}$ (a.u.)"),
    ("SEII_25", "1", "$SEII_{25}$"),
    ("A", "a.u. * kHZ", "$A$ (a.u.)"),
])
def test_metric_curves_preserve_phase_coordinates_nan_gaps_and_units(metric, unit, label):
    phase = np.array([0, 0.1, 0.5, 0.9])
    values = np.array([4, np.nan, 2, 3])
    fig = _gap_metric_figure(values, phase, metric=metric, unit=unit)
    ax = fig.axes[0]
    np.testing.assert_array_equal(ax.lines[0].get_xdata(), phase)
    np.testing.assert_array_equal(ax.lines[0].get_ydata(), values)
    assert ax.get_xlim() == (0, 1)
    assert ax.get_xlabel() == r"Cardiac phase ($t/T$)"
    assert ax.get_ylabel() == label
    for text in (ax.xaxis.label, ax.yaxis.label, *ax.get_xticklabels(), *ax.get_yticklabels()):
        assert text.get_fontfamily() == ["Times New Roman"]
    assert ax.xaxis.label.get_math_fontfamily() == ax.yaxis.label.get_math_fontfamily() == "stix"
    fig.canvas.draw()
    assert all(tick.get_rotation() == 0 for tick in ax.get_xticklabels())
    assert all(tick.get_fontsize() == 8 for tick in ax.get_xticklabels())
    assert all(spine.get_linewidth() == 0.4 for spine in ax.spines.values())
    np.testing.assert_allclose(fig.get_size_inches(), [
        FIGURE_WIDTH_INCHES, FIGURE_WIDTH_INCHES / FIGURE_ASPECT_RATIO,
    ])


def test_unbinned_frequency_helper_uses_sampling_rate_and_actual_bin_domain():
    for count in (32, 64, 128):
        frequencies = frequency_coordinates_khz(count, 64000)
        assert frequencies[0] == 0
        assert frequencies[count // 2] == -32
        assert frequencies[count // 2 - 1] == 32 - 64 / count
        assert frequencies[-1] == -64 / count
        np.testing.assert_allclose(frequencies, np.fft.fftfreq(count, d=1 / 64000) / 1000)


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
    np.testing.assert_array_equal(first, original)


@pytest.mark.parametrize("bin_count", [31, 64])
def test_frequency_helpers_and_full_spectrum_preserve_unbinned_domain(bin_count):
    stop = (bin_count + 1) // 2
    signed = frequency_coordinates_khz(bin_count, 64000)
    coordinates = signed[signed >= 0]
    np.testing.assert_allclose(coordinates, np.arange(stop) * 64 / bin_count)
    assert coordinates[0] == 0
    np.testing.assert_array_equal(coordinates, frequency_coordinates_khz(bin_count, 64000)[:stop])
    values = np.ones(bin_count)
    values[3] = 10
    curve = _spectrum_curve(
        values, np.ones(bin_count), frequencies_khz=frequency_coordinates_khz(bin_count, 64000),
        label="Median", linestyle="-", index_max1=5, index_max2=bin_count - 2,
    )
    np.testing.assert_array_equal(curve.frequencies_khz, frequency_coordinates_khz(bin_count, 64000))
    assert curve.power.size == bin_count
    # The peak does not affect normalization by the central-quarter Hv.
    np.testing.assert_array_equal(curve.power, np.r_[np.ones(3), 10, np.ones(bin_count - 4)])


@pytest.mark.parametrize("bin_count, end, positive_count", [(64, 0.5, 32), (48, 0.25, 16)])
def test_frequency_figures_label_the_full_averaged_signed_domain(bin_count, end, positive_count):
    frequencies = frequency_coordinates_khz(
        bin_count, 37037, fft_bin_count=256, binning_factor=4, frequency_end=end,
    )
    spectrum = np.ones(bin_count)
    median = _gap_median_figure(spectrum, frequencies_khz=frequencies).axes[0]
    heatmap = _gap_spectrogram_figure(
        spectrum[None, :], frequencies_khz=frequencies, dt_seconds=0.01,
    ).axes[0]
    for positions, labels, whole_numbers in (
        (median.get_xticks(), median.get_xticklabels(), False),
        (heatmap.get_yticks(), heatmap.get_yticklabels(), True),
    ):
        ticks = dict(zip(positions, (label.get_text() for label in labels)))
        assert ticks[0] == ("0" if whole_numbers else "0.217")
        assert ticks[positive_count] == ("-18" if whole_numbers else "-18.3")
        for position, label in ticks.items():
            frequency = frequencies[int(position)]
            assert label == (str(round(float(frequency))) if whole_numbers else f"{frequency:.3g}")
    assert median.get_xticklabels()[-1].get_text() == "-0.362"
    assert heatmap.get_yticks()[-1] == bin_count - 1
    assert heatmap.get_yticklabels()[-1].get_text() == "0"
    heatmap.figure.set_dpi(320)
    heatmap.figure.canvas.draw()
    renderer = heatmap.figure.canvas.get_renderer()
    bounds = [label.get_window_extent(renderer) for label in heatmap.get_yticklabels()]
    assert all(bottom.y1 < top.y0 for bottom, top in zip(bounds, bounds[1:]))
    curve = _spectrum_curve(
        spectrum, np.ones(bin_count), frequencies_khz=frequencies,
        label="Median", linestyle="-", index_max1=5, index_max2=bin_count - 2,
    )
    np.testing.assert_array_equal(curve.frequencies_khz, frequencies)
    assert curve.frequencies_khz[0] == pytest.approx(0.217013671875)
    assert frequencies[curve.peak_bins[0]] == pytest.approx(3.110529296875)


@pytest.mark.parametrize("selection", [(0,), (1, 2), (0, 1, 2)])
def test_full_spectrum_styles_hv_shading_peaks_zero_legend_and_double_width(selection, caplog):
    curves = []
    for label, style, peak, scale, deviation in (
        ("Median", "-", 4, 10, 2), ("Early", "--", 5, 20, 3), ("Late", ":", 6, 30, 36),
    ):
        values = np.linspace(4, 1, 64)
        values[peak] = scale
        curves.append(_spectrum_curve(
            values, np.full(64, deviation), frequencies_khz=frequency_coordinates_khz(64, 64000),
            label=label, linestyle=style, index_max1=4, index_max2=59,
        ))
    curves = [curves[index] for index in selection]
    fig = _spectra_figure(curves)
    ax = fig.axes[0]
    centers = {line.get_label(): line for line in ax.lines if line.get_label() in {"median", "early", "late"}}
    assert {name: line.get_linestyle() for name, line in centers.items()} == {
        curve.label.lower(): curve.linestyle for curve in curves
    }
    assert all(line.get_color() == "black" for line in centers.values())
    zero = next(line for line in ax.lines if np.all(np.asarray(line.get_ydata()) == 0))
    assert zero.get_linestyle() == ":"
    assert len(ax.lines) == len(curves) + 1  # Signals and zero; peaks use ordinary axis ticks.
    region = next(patch for patch in ax.patches if patch.get_alpha() is not None)
    assert region.get_x() == 23.5 and region.get_width() == 16
    assert 0 < region.get_alpha() <= 0.15
    assert ax.get_xlabel() == "Frequency (kHz)" and ax.get_ylabel() == "Power (a.u.)"
    assert ax.get_yscale() == "linear"
    assert ax.get_xlim() == (0, 63)
    frequency_labels = dict(zip(ax.get_xticks(), (tick.get_text() for tick in ax.get_xticklabels())))
    assert frequency_labels[4] == "4" and frequency_labels[59] == "-5"
    assert all("." not in label for label in frequency_labels.values())
    assert all("." not in tick.get_text() for tick in ax.get_yticklabels())
    np.testing.assert_allclose(fig.get_size_inches(), [
        2 * FIGURE_WIDTH_INCHES, FIGURE_WIDTH_INCHES / FIGURE_ASPECT_RATIO * SPECTRUM_HEIGHT_SCALE,
    ])
    assert ax.xaxis.label.get_fontsize() == ax.yaxis.label.get_fontsize() == SPECTRUM_FONT_SIZE
    assert all(tick.get_fontsize() == SPECTRUM_TICK_FONT_SIZE for tick in ax.get_xticklabels())
    assert all(tick.get_fontsize() == SPECTRUM_TICK_FONT_SIZE for tick in ax.get_yticklabels())
    assert ax.get_ylim()[0] < min(0, *(np.nanmin(curve.power - curve.std) for curve in curves))
    assert ax.get_ylim()[1] > max(0, *(np.nanmax(curve.power + curve.std) for curve in curves))
    assert len(ax.collections) == len(curves)
    region_rgba = np.asarray(region.get_facecolor())
    region_gray = np.mean(region_rgba[:3] * region_rgba[3] + 1 - region_rgba[3])
    for band in ax.collections:
        rgba = band.get_facecolors()[0]
        band_gray = np.mean(rgba[:3] * rgba[3] + 1 - rgba[3])
        assert 0.55 < band_gray < region_gray - 0.15
        assert band.get_zorder() < min(line.get_zorder() for line in centers.values())
        assert not band.get_hatch()
        assert band.get_edgecolors().size == 0
    legend = ax.get_legend()
    assert [text.get_text() for text in legend.get_texts()] == [curve.label.lower() for curve in curves]
    assert all(text.get_fontsize() == ax.xaxis.label.get_fontsize() for text in legend.get_texts())
    assert all(label.get_rotation() == 0 for label in ax.get_xticklabels())
    eps = BytesIO()
    with caplog.at_level("WARNING", logger="matplotlib.backends.backend_ps"):
        fig.savefig(eps, format="eps", dpi=320)
    assert b"colorimage" in eps.getvalue()
    assert "does not support transparency" not in caplog.text


@pytest.mark.parametrize("values", [np.zeros(64), np.full(64, np.nan)])
def test_unusable_normalization_does_not_fabricate_a_spectrum(values):
    assert _spectrum_curve(values, np.zeros(64), frequencies_khz=frequency_coordinates_khz(64, 64000), label="Median", linestyle="-", index_max1=0, index_max2=63) is None


@pytest.mark.parametrize("bin_count, start, stop", [
    (2, 0, 2), (31, 11, 19), (64, 24, 40),
])
def test_normalization_uses_central_quarter_hv_of_the_full_spectrum(bin_count, start, stop):
    values = np.full(bin_count, 10000.0)
    values[start:stop] = np.arange(1, stop - start + 1)
    assert _spectrum_baseline(values) == (stop - start + 1) / 2


def test_normalization_median_ignores_nan_without_expanding_the_neighborhood():
    values = np.full(64, 10000.0)
    values[24:40] = np.nan
    values[29:32] = [4, np.nan, 12]
    assert _spectrum_baseline(values) == 8
    values[24:40] = np.nan
    assert _spectrum_baseline(values) is None


@pytest.mark.parametrize("products", ["all", "empty_mask", "one_metric", "one_spectrum", "none"])
def test_exports_only_available_products_in_shared_namespace(tmp_path, products):
    manager = OutputManager(HoloRunLayout.from_holo(tmp_path / "scan.holo"))
    schema = EyeFlowOutputPaths.active()
    expected = set()
    with h5py.File(tmp_path / "out.h5", "w") as h5:
        h5.create_dataset(schema.sampling_frequency_hz, data=8000)
        _write_frequency_table(h5, 32 if products in {"all", "one_spectrum"} else 4, 8000)
        if products in {"all", "raw"}:
            h5.create_dataset(schema.spectrograms_masked, data=np.arange(192).reshape(6, 32) if products == "all" else np.arange(24).reshape(6, 4))
            expected.update({"spectrograms_gaps", "spectrograms_gaps_monochrome"})
        if products in {"all", "one_spectrum"}:
            values = np.ones(32)
            values[2] = 8
            h5.create_dataset(schema.spectrograms_masked_time_aggregated_median,
                              data=values)
            h5.create_dataset(schema.spectrograms_masked_time_aggregated_median_beat,
                              data=np.full((128, 32), 100))
            expected.add("spectrogram_gaps_median")
            h5.create_dataset(schema.spectrograms_masked_time_aggregated_std, data=np.ones(32))
            h5.create_dataset(f"{schema.spectrograms_masked_metrics_root}/index_max1", data=2)
            h5.create_dataset(f"{schema.spectrograms_masked_metrics_root}/index_max2", data=29)
            expected.add("gaps_positive_freq_median_phase_beat")
        if products == "all":
            early, late = np.ones(32), np.ones(32)
            early[4], late[3] = 9, 12
            h5.create_dataset(schema.spectrograms_masked_time_aggregated_median_early, data=early)
            h5.create_dataset(schema.spectrograms_masked_time_aggregated_median_late, data=late)
            h5.create_dataset(schema.spectrograms_masked_time_aggregated_std_early, data=np.ones(32))
            h5.create_dataset(schema.spectrograms_masked_time_aggregated_std_late, data=np.ones(32))
            expected.update({
                "gaps_positive_freq_median_phase_early_late", "gaps_positive_freq_median_phase_all",
            })
        if products in {"all", "mask", "empty_mask"}:
            mask = np.eye(4, dtype=bool) if products != "empty_mask" else np.zeros((4, 4), bool)
            h5.create_dataset(schema.segmentation.gaps, data=mask.astype(np.uint8))
            expected.add("gaps_mask")
        if products in {"all", "metrics", "one_metric"}:
            names = ("A0_t",) if products == "one_metric" else (
                "A0_t", "A_t", "SE", "A_25_t", "A_50_t", "A_75_t", "SE_25", "SE_50", "SE_75",
                "SEII_25_t", "SEII_50_t", "SEII_75_t",
            )
            for name in names:
                metric = h5.create_dataset(f"{schema.spectrograms_masked_metrics_root}/{name}", data=[2, np.nan, 1, 3])
                metric.attrs["phase_samples"] = [0, 0.25, 0.5, 0.75]
                metric.attrs["unit"] = "1" if name.startswith(("R", "SE")) else "a.u."
                expected.add(f"gaps_{name.removesuffix('_t')}")
        paths = export_spectrogram_figures(SimpleNamespace(
            manager=manager, h5=PipelineH5Output(h5),
        ), dt_seconds=0.01)

    assert {path.relative_to(manager.layout.ef_dir).as_posix() for path in paths} == {
        f"{extension}/spectrograms/scan_{name}.{extension}"
        for extension in ("png", "eps") for name in expected
    }
    for path in paths:
        if path.suffix == ".eps":
            assert path.read_bytes().startswith(b"%!PS-Adobe")
        else:
            with Image.open(path) as image:
                assert image.format == "PNG"
                ratio = 1 if "gaps_mask" in path.name else FIGURE_ASPECT_RATIO
                if "spectrograms_gaps" in path.name:
                    ratio = SPECTROGRAM_ASPECT_RATIO
                elif "positive_freq" in path.name:
                    ratio = 2 * FIGURE_ASPECT_RATIO / SPECTRUM_HEIGHT_SCALE
                assert np.isclose(image.width / image.height, ratio, rtol=0.01)
                image.verify()
            if "monochrome" in path.name:
                with Image.open(path) as image:
                    rgb = np.asarray(image.convert("RGB"))
                    np.testing.assert_array_equal(rgb[..., 0], rgb[..., 1])
                    np.testing.assert_array_equal(rgb[..., 1], rgb[..., 2])
            if products == "empty_mask" and "gaps_mask" in path.name:
                with Image.open(path) as image:
                    assert not np.any(np.asarray(image.convert("RGB")))


@pytest.mark.parametrize("dt_seconds", [None, 0, np.nan])
def test_invalid_frame_interval_skips_heatmaps_but_preserves_other_figures(tmp_path, dt_seconds):
    schema = EyeFlowOutputPaths.active()
    with h5py.File(tmp_path / "out.h5", "w") as h5:
        _write_frequency_table(h5, 4, 8000)
        h5.create_dataset(schema.spectrograms_masked, data=np.ones((6, 4)))
        h5.create_dataset(schema.spectrograms_masked_time_aggregated_median, data=np.ones(4))
        h5.create_dataset(schema.segmentation.gaps, data=np.eye(4))
        output = SimpleNamespace(
            manager=OutputManager(HoloRunLayout.from_holo(tmp_path / "scan.holo")),
            h5=PipelineH5Output(h5),
        )
        with patch("pipelines.low_freq_spectra.artifacts._save_figure_pair", return_value=[]) as save:
            export_spectrogram_figures(output, dt_seconds=dt_seconds)
        assert {call.args[1] for call in save.call_args_list} == {"spectrogram_gaps_median", "gaps_mask"}
        assert schema.spectrograms_masked in h5


@pytest.mark.parametrize("frequency", [None, np.nan])
def test_missing_calibration_preserves_mask_exports_without_mislabeled_spectral_axes(tmp_path, frequency):
    manager = OutputManager(HoloRunLayout.from_holo(tmp_path / "scan.holo"))
    schema = EyeFlowOutputPaths.active()
    with h5py.File(tmp_path / "out.h5", "w") as h5:
        h5.create_dataset(schema.spectrograms_masked, data=np.ones((2, 64)))
        h5.create_dataset(schema.segmentation.gaps, data=np.eye(4))
        if frequency is not None:
            h5.create_dataset(schema.sampling_frequency_hz, data=frequency)
        paths = export_spectrogram_figures(SimpleNamespace(manager=manager, h5=PipelineH5Output(h5)))
        assert {path.name for path in paths} == {"scan_gaps_mask.png", "scan_gaps_mask.eps"}
        assert schema.spectrograms_masked in h5


def test_export_reads_shared_peak_without_mutating_hdf5(tmp_path):
    manager = OutputManager(HoloRunLayout.from_holo(tmp_path / "scan.holo"))
    schema = EyeFlowOutputPaths.active()
    with h5py.File(tmp_path / "out.h5", "w") as h5:
        h5.create_dataset(schema.sampling_frequency_hz, data=64000)
        _write_frequency_table(h5, 64, 64000)
        peak = h5.create_dataset(f"{schema.spectrograms_masked_metrics_root}/index_max1", data=5)
        peak.attrs["valley_neighborhood_length"] = 4
        h5.create_dataset(f"{schema.spectrograms_masked_metrics_root}/index_max2", data=55)
        spectra = {}
        for label, center_path, std_path, peak, maximum, deviation, baseline in (
            ("Median", schema.spectrograms_masked_time_aggregated_median, schema.spectrograms_masked_time_aggregated_std, 5, 1000, 1, 2),
            ("Early", schema.spectrograms_masked_time_aggregated_median_early, schema.spectrograms_masked_time_aggregated_std_early, 7, 20, 2, 4),
            ("Late", schema.spectrograms_masked_time_aggregated_median_late, schema.spectrograms_masked_time_aggregated_std_late, 2, 40, 4, 8),
        ):
            values = np.ones(64)
            values[peak] = maximum
            values[12] = np.nan
            values[40] = 10000  # Outside the central-quarter Hv neighborhood.
            values[24:40] = 100  # Distinguish the default quarter from the recorded 4-bin window.
            values[30:34] = baseline
            values[29:32] = [baseline, np.nan, baseline]
            std = np.full(64, deviation)
            spectra[label] = (values, std, baseline)
            h5.create_dataset(center_path, data=values)
            h5.create_dataset(std_path, data=std)
        original = {}
        h5.visititems(lambda name, obj: original.update({name: obj[()]}) if isinstance(obj, h5py.Dataset) else None)
        with patch("pipelines.low_freq_spectra.artifacts._spectra_figure", wraps=_spectra_figure) as draw:
            paths = export_spectrogram_figures(SimpleNamespace(manager=manager, h5=PipelineH5Output(h5)))
        assert draw.call_count == 3  # Build each figure once for its PNG/EPS pair.
        for call in draw.call_args_list:
            for curve in call.args[0]:
                np.testing.assert_array_equal(curve.frequencies_khz, frequency_coordinates_khz(
                    64, 64000, fft_bin_count=256, binning_factor=4,
                ))
                assert curve.peak_bins == (5, 55)
                np.testing.assert_array_equal(curve.valley_bins, np.arange(30, 34))
                values, std, own_median = spectra[curve.label]
                np.testing.assert_allclose(curve.power, values / own_median)
                np.testing.assert_allclose(curve.std, std / own_median)
        current = {}
        h5.visititems(lambda name, obj: current.update({name: obj[()]}) if isinstance(obj, h5py.Dataset) else None)
        assert original.keys() == current.keys()
        for name in original:
            np.testing.assert_array_equal(original[name], current[name])
        assert len([path for path in paths if "positive_freq" in path.name]) == 6


@pytest.mark.parametrize("half, peak", [
    (1, None), (2, None), (1, np.nan), (2, -1), (1, 32), (2, 16.5),
])
def test_spectrum_export_does_not_infer_missing_or_invalid_stored_peaks(tmp_path, half, peak):
    manager = OutputManager(HoloRunLayout.from_holo(tmp_path / "scan.holo"))
    schema = EyeFlowOutputPaths.active()
    with h5py.File(tmp_path / "out.h5", "w") as h5:
        h5.create_dataset(schema.sampling_frequency_hz, data=64000)
        _write_frequency_table(h5, 32, 64000)
        h5.create_dataset(schema.spectrograms_masked_time_aggregated_median, data=np.ones(32))
        h5.create_dataset(schema.spectrograms_masked_time_aggregated_std, data=np.ones(32))
        for side, valid in ((1, 5), (2, 25)):
            value = peak if side == half else valid
            if value is not None:
                h5.create_dataset(f"{schema.spectrograms_masked_metrics_root}/index_max{side}", data=value)
        paths = export_spectrogram_figures(SimpleNamespace(manager=manager, h5=PipelineH5Output(h5)))
        assert not any("positive_freq" in path.name for path in paths)
