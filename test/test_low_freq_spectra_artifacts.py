"""Spectral artifact axes, native mask orientation, and shared output routing."""

from types import SimpleNamespace
from unittest.mock import patch

import h5py
import numpy as np
import pytest
from matplotlib.colors import Normalize
from PIL import Image

from input_output.h5_access import PipelineH5Output
from input_output.holo_run_layout import HoloRunLayout
from input_output.output_manager import OutputManager
from input_output.schema import EyeFlowOutputPaths
from pipelines.low_freq_spectra.artifacts import (
    FIGURE_ASPECT_RATIO,
    FIGURE_WIDTH_INCHES,
    _gap_median_figure,
    _gap_metric_figure,
    _gap_spectrogram_figure,
    _gaps_mask_figure,
    _positive_spectrum_curve,
    _positive_spectra_figure,
    export_spectrogram_figures,
)
from pipelines.low_freq_spectra.frequency import frequency_coordinates_khz, positive_frequency_coordinates_khz


@pytest.mark.parametrize("frames, frame_ticks", [(103, [0, 50, 100]), (640, list(range(0, 640, 100)))])
def test_heatmap_preserves_frame_bin_values_and_nan_gaps(frames, frame_ticks):
    values = np.arange(frames * 32, dtype=float).reshape(frames, 32)
    values[51, 2] = np.nan
    original = values.copy()
    fig = _gap_spectrogram_figure(values, sampling_frequency_hz=64000)
    ax, colorbar = fig.axes
    image = ax.images[0]

    np.testing.assert_array_equal(image.get_array().filled(np.nan), values.T)
    np.testing.assert_array_equal(values, original)
    np.testing.assert_array_equal(ax.get_xticks(), frame_ticks)
    np.testing.assert_array_equal(ax.get_yticks(), [0, 4, 8, 12, 16, 20, 24, 28, 31])
    assert [label.get_text() for label in ax.get_yticklabels()] == ["0", "4", "8", "12", "-16", "-12", "-8", "-4", "-1"]
    assert ax.get_xlabel() == "Time (frames)"
    assert ax.get_ylabel() == "Frequency (kHz)"
    assert image.origin == "lower"
    assert image.get_cmap().name == "viridis"
    assert type(image.norm) is Normalize
    assert colorbar.get_ylabel() == "Power (a.u.)"
    assert image.colorbar.orientation == "vertical"
    np.testing.assert_allclose(fig.get_size_inches(), [
        FIGURE_WIDTH_INCHES, FIGURE_WIDTH_INCHES / FIGURE_ASPECT_RATIO,
    ])
    for axis in fig.axes:
        assert all(spine.get_linewidth() == 0.4 for spine in axis.spines.values())
        assert axis.yaxis.label.get_fontsize() == 10
        assert axis.yaxis.label.get_fontfamily() == ["Times New Roman"]
        assert all(label.get_fontsize() == 8 for label in axis.get_xticklabels())
        assert all(label.get_rotation() == 0 for label in axis.get_xticklabels())
        assert all(label.get_fontsize() == 8 for label in axis.get_yticklabels())
        assert all(label.get_rotation() == 0 for label in axis.get_yticklabels())
    fig.canvas.draw()
    bounds = [label.get_window_extent(fig.canvas.get_renderer()) for label in ax.get_xticklabels()]
    assert all(left.x1 < right.x0 for left, right in zip(bounds, bounds[1:]))


def test_median_labels_signed_khz_without_reordering_or_filling_nan():
    values = np.array([20, 5, np.nan, 10])
    ax = _gap_median_figure(values, sampling_frequency_hz=8000).axes[0]
    np.testing.assert_array_equal(ax.lines[0].get_xdata(), [0, 1, 2, 3])
    np.testing.assert_array_equal(ax.lines[0].get_ydata(), values)
    np.testing.assert_array_equal(ax.get_xticks(), [0, 1, 2, 3])
    assert [label.get_text() for label in ax.get_xticklabels()] == ["0", "1", "-2", "-1"]
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


@pytest.mark.parametrize("metric, unit, label", [("A0", "power_unit", "A0 (power_unit)"), ("R1", "1", "R1")])
def test_metric_curves_preserve_phase_coordinates_nan_gaps_and_units(metric, unit, label):
    phase = np.array([0, 0.1, 0.5, 0.9])
    values = np.array([4, np.nan, 2, 3])
    fig = _gap_metric_figure(values, phase, metric=metric, unit=unit)
    ax = fig.axes[0]
    np.testing.assert_array_equal(ax.lines[0].get_xdata(), phase)
    np.testing.assert_array_equal(ax.lines[0].get_ydata(), values)
    assert ax.get_xlim() == (0, 1)
    assert ax.get_xlabel() == "Cardiac phase (t/T)"
    assert ax.get_ylabel() == label
    assert all(tick.get_rotation() == 0 for tick in ax.get_xticklabels())
    assert all(tick.get_fontsize() == 8 for tick in ax.get_xticklabels())
    assert all(spine.get_linewidth() == 0.4 for spine in ax.spines.values())
    np.testing.assert_allclose(fig.get_size_inches(), [
        FIGURE_WIDTH_INCHES, FIGURE_WIDTH_INCHES / FIGURE_ASPECT_RATIO,
    ])


def test_frequency_calibration_uses_hd_sampling_rate_and_actual_bin_domain():
    for count in (32, 64, 128):
        frequencies = frequency_coordinates_khz(count, 64000)
        assert frequencies[0] == 0
        assert frequencies[count // 2] == -16
        assert frequencies[count // 2 - 1] == 16 - 32 / count
        assert frequencies[-1] == -32 / count


def test_positive_spectra_share_the_stored_peak_cutoff_and_normalize_their_own_std():
    first, second = np.ones(64), np.ones(64)
    first[5], second[7] = 10, 20
    first[40], second[50] = 10000, 20000  # Negative peaks cannot affect positive normalization.
    first[12] = np.nan
    original = first.copy()
    a = _positive_spectrum_curve(first, np.full(64, 2), sampling_frequency_hz=64000, label="Early", linestyle="--", index_max1=5)
    b = _positive_spectrum_curve(second, np.full(64, 5), sampling_frequency_hz=64000, label="Late", linestyle=":", index_max1=5)
    np.testing.assert_array_equal(a.frequencies_khz, np.arange(5, 32) / 2)
    np.testing.assert_array_equal(b.frequencies_khz, np.arange(5, 32) / 2)
    assert a.power[0] == b.power[2] == 1
    assert b.power[0] == 1 / 20
    np.testing.assert_array_equal(a.std, [0.2] * 27)
    np.testing.assert_array_equal(b.std, [0.25] * 27)
    assert np.isnan(a.power[7])
    np.testing.assert_array_equal(first, original)


@pytest.mark.parametrize("bin_count", [32, 64, 128])
def test_positive_frequency_calibration_always_anchors_bin_32_at_quarter_sampling_rate(bin_count):
    coordinates = positive_frequency_coordinates_khz(64000)
    np.testing.assert_array_equal(coordinates, np.linspace(0, 16, 33))
    assert coordinates[0] == 0
    assert coordinates[32] == 16
    values = np.ones(bin_count)
    values[3] = 10
    curve = _positive_spectrum_curve(
        values, np.ones(bin_count), sampling_frequency_hz=64000,
        label="Median", linestyle="-", index_max1=5,
    )
    np.testing.assert_array_equal(curve.frequencies_khz, np.arange(5, 32) / 2)
    assert curve.frequencies_khz[-1] == 15.5
    assert curve.power.size == 27
    # Normalize by the curve's own maximum even when it precedes the shared cutoff.
    np.testing.assert_array_equal(curve.power, np.full(27, 0.1))


def test_positive_plot_styles_envelopes_zero_reference_and_data_margins():
    curves = []
    for label, style, peak, scale, deviation in (
        ("Median", "-", 4, 10, 2), ("Early", "--", 5, 20, 3), ("Late", ":", 6, 30, 36),
    ):
        values = np.linspace(4, 1, 64)
        values[peak] = scale
        curves.append(_positive_spectrum_curve(
            values, np.full(64, deviation), sampling_frequency_hz=64000,
            label=label, linestyle=style, index_max1=4,
        ))
    fig = _positive_spectra_figure(curves)
    ax = fig.axes[0]
    centers = {line.get_label(): line for line in ax.lines if line.get_label() in {"Median", "Early", "Late"}}
    assert {name: line.get_linestyle() for name, line in centers.items()} == {"Median": "-", "Early": "--", "Late": ":"}
    assert all(line.get_color() == "black" for line in centers.values())
    assert any(line.get_linestyle() == ":" and np.all(np.asarray(line.get_ydata()) == 0) for line in ax.lines)
    assert ax.get_xlabel() == "Frequency (kHz)" and ax.get_ylabel() == "Power (a.u.)"
    assert ax.get_yscale() == "linear"
    assert ax.get_xlim()[0] < 2 and ax.get_xlim()[1] > 15.5
    assert ax.get_ylim()[0] < min(np.nanmin(curve.power - curve.std) for curve in curves)
    assert ax.get_ylim()[1] > max(np.nanmax(curve.power + curve.std) for curve in curves)
    # Three filled bands plus dashed and dotted clipped interior patterns.
    assert len(ax.collections) == 5
    assert all(label.get_rotation() == 0 for label in ax.get_xticklabels())


@pytest.mark.parametrize("values", [np.zeros(64), np.full(64, np.nan)])
def test_unusable_normalization_does_not_fabricate_a_positive_spectrum(values):
    assert _positive_spectrum_curve(values, np.zeros(64), sampling_frequency_hz=64000, label="Median", linestyle="-", index_max1=0) is None


@pytest.mark.parametrize("products", ["all", "raw", "mask", "empty_mask", "metrics", "one_metric", "one_spectrum", "none"])
def test_exports_only_available_products_in_shared_namespace(tmp_path, products):
    manager = OutputManager(HoloRunLayout.from_holo(tmp_path / "scan.holo"))
    schema = EyeFlowOutputPaths.active()
    expected = set()
    with h5py.File(tmp_path / "out.h5", "w") as h5:
        h5.create_dataset(schema.sampling_frequency_hz, data=8000)
        if products in {"all", "raw"}:
            h5.create_dataset(schema.spectrograms_masked, data=np.arange(24).reshape(6, 4))
            expected.add("spectrograms_gaps")
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
            expected.add("gaps_positive_freq_median_phase_beat")
        if products == "all":
            early, late = np.ones(32), np.ones(32)
            early[4], late[3] = 9, 12
            h5.create_dataset(schema.spectrograms_masked_time_aggregated_median_early, data=early)
            h5.create_dataset(schema.spectrograms_masked_time_aggregated_median_late, data=late)
            h5.create_dataset(schema.spectrograms_masked_time_aggregated_std_early, data=np.ones(32))
            h5.create_dataset(schema.spectrograms_masked_time_aggregated_std_late, data=np.ones(32))
            expected.update({"gaps_positive_freq_median_phase_early_late", "gaps_positive_freq_median_phase_all"})
        if products in {"all", "mask", "empty_mask"}:
            mask = np.eye(4, dtype=bool) if products != "empty_mask" else np.zeros((4, 4), bool)
            h5.create_dataset(schema.segmentation.gaps, data=mask.astype(np.uint8))
            expected.add("gaps_mask")
        if products in {"all", "metrics", "one_metric"}:
            names = ("A0",) if products == "one_metric" else ("A0", "A1", "A2", "R1", "R2")
            for name in names:
                metric = h5.create_dataset(f"{schema.spectrograms_masked_metrics_root}/{name}_t", data=[2, np.nan, 1, 3])
                metric.attrs["phase_samples"] = [0, 0.25, 0.5, 0.75]
                metric.attrs["unit"] = "1" if name.startswith("R") else "a.u."
                expected.add(f"gaps_{name}")
        paths = export_spectrogram_figures(SimpleNamespace(
            manager=manager, h5=PipelineH5Output(h5),
        ))

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
                assert np.isclose(image.width / image.height, ratio, rtol=0.01)
                image.verify()
            if products == "empty_mask" and "gaps_mask" in path.name:
                with Image.open(path) as image:
                    assert not np.any(np.asarray(image.convert("RGB")))


@pytest.mark.parametrize("frequency", [None, 0, np.nan])
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


def test_export_reads_shared_peak_from_h5_and_does_not_save_frequency_coordinates(tmp_path):
    manager = OutputManager(HoloRunLayout.from_holo(tmp_path / "scan.holo"))
    schema = EyeFlowOutputPaths.active()
    with h5py.File(tmp_path / "out.h5", "w") as h5:
        h5.create_dataset(schema.sampling_frequency_hz, data=64000)
        h5.create_dataset(f"{schema.spectrograms_masked_metrics_root}/index_max1", data=5)
        for center_path, std_path, peak in (
            (schema.spectrograms_masked_time_aggregated_median, schema.spectrograms_masked_time_aggregated_std, 5),
            (schema.spectrograms_masked_time_aggregated_median_early, schema.spectrograms_masked_time_aggregated_std_early, 7),
            (schema.spectrograms_masked_time_aggregated_median_late, schema.spectrograms_masked_time_aggregated_std_late, 2),
        ):
            values = np.ones(32)
            values[peak] = 10
            h5.create_dataset(center_path, data=values)
            h5.create_dataset(std_path, data=np.ones(32))
        original = {}
        h5.visititems(lambda name, obj: original.update({name: obj[()]}) if isinstance(obj, h5py.Dataset) else None)
        with patch("pipelines.low_freq_spectra.artifacts._positive_spectra_figure", wraps=_positive_spectra_figure) as draw:
            paths = export_spectrogram_figures(SimpleNamespace(manager=manager, h5=PipelineH5Output(h5)))
        assert draw.call_count == 6  # Three figures, PNG and EPS.
        for call in draw.call_args_list:
            for curve in call.args[0]:
                np.testing.assert_array_equal(curve.frequencies_khz, np.arange(5, 32) / 2)
        current = {}
        h5.visititems(lambda name, obj: current.update({name: obj[()]}) if isinstance(obj, h5py.Dataset) else None)
        assert original.keys() == current.keys()
        for name in original:
            np.testing.assert_array_equal(original[name], current[name])
        assert len([path for path in paths if "positive_freq" in path.name]) == 6


@pytest.mark.parametrize("peak", [None, np.nan, -1, 32, 5.5])
def test_positive_export_does_not_infer_a_peak_if_the_stored_index_is_missing_or_invalid(tmp_path, peak):
    manager = OutputManager(HoloRunLayout.from_holo(tmp_path / "scan.holo"))
    schema = EyeFlowOutputPaths.active()
    with h5py.File(tmp_path / "out.h5", "w") as h5:
        h5.create_dataset(schema.sampling_frequency_hz, data=64000)
        h5.create_dataset(schema.spectrograms_masked_time_aggregated_median, data=np.ones(32))
        h5.create_dataset(schema.spectrograms_masked_time_aggregated_std, data=np.ones(32))
        if peak is not None:
            h5.create_dataset(f"{schema.spectrograms_masked_metrics_root}/index_max1", data=peak)
        paths = export_spectrogram_figures(SimpleNamespace(manager=manager, h5=PipelineH5Output(h5)))
        assert not any("positive_freq" in path.name for path in paths)
