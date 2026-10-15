"""Spectral artifact axes, native mask orientation, and shared output routing."""

from types import SimpleNamespace

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
    _gap_spectrogram_figure,
    _gaps_mask_figure,
    export_spectrogram_figures,
)


def test_heatmap_preserves_frame_bin_values_and_nan_gaps():
    values = np.arange(103 * 32, dtype=float).reshape(103, 32)
    values[51, 2] = np.nan
    original = values.copy()
    fig = _gap_spectrogram_figure(values)
    ax, colorbar = fig.axes
    image = ax.images[0]

    np.testing.assert_array_equal(image.get_array().filled(np.nan), values.T)
    np.testing.assert_array_equal(values, original)
    np.testing.assert_array_equal(ax.get_xticks(), [0, 50, 100])
    np.testing.assert_array_equal(ax.get_yticks(), np.arange(0, 32, 4))
    assert ax.get_xlabel() == "Time (frames)"
    assert ax.get_ylabel() == "Frequency band"
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
        assert axis.yaxis.label.get_fontsize() == 8
        assert axis.yaxis.label.get_fontfamily() == ["Times New Roman"]
        assert all(label.get_fontsize() == 5 for label in axis.get_xticklabels())
        assert all(label.get_rotation() == 0 for label in axis.get_xticklabels())
        assert all(label.get_fontsize() == 5 for label in axis.get_yticklabels())
        assert all(label.get_rotation() == 0 for label in axis.get_yticklabels())


def test_median_uses_integer_bins_without_reordering_or_filling_nan():
    values = np.array([20, 5, np.nan, 10])
    ax = _gap_median_figure(values).axes[0]
    np.testing.assert_array_equal(ax.lines[0].get_xdata(), [0, 1, 2, 3])
    np.testing.assert_array_equal(ax.lines[0].get_ydata(), values)
    np.testing.assert_array_equal(ax.get_xticks(), [0, 2])
    assert all(label.get_fontsize() == 5 for label in ax.get_xticklabels())
    assert all(label.get_rotation() == 90 for label in ax.get_xticklabels())
    assert all(label.get_fontsize() == 5 for label in ax.get_yticklabels())
    assert all(label.get_rotation() == 0 for label in ax.get_yticklabels())
    assert ax.get_xlabel() == "Frequency band"
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


@pytest.mark.parametrize("products", ["all", "raw", "mask", "empty_mask", "none"])
def test_exports_only_available_products_in_shared_namespace(tmp_path, products):
    manager = OutputManager(HoloRunLayout.from_holo(tmp_path / "scan.holo"))
    schema = EyeFlowOutputPaths.active()
    expected = set()
    with h5py.File(tmp_path / "out.h5", "w") as h5:
        if products in {"all", "raw"}:
            h5.create_dataset(schema.spectrograms_masked, data=np.arange(24).reshape(6, 4))
            expected.add("spectrograms_gaps")
        if products == "all":
            h5.create_dataset(schema.spectrograms_masked_time_aggregated_median,
                              data=[8, 3, np.nan, 5])
            expected.add("spectrogram_gaps_median")
        if products in {"all", "mask", "empty_mask"}:
            mask = np.eye(4, dtype=bool) if products != "empty_mask" else np.zeros((4, 4), bool)
            h5.create_dataset(schema.segmentation.gaps, data=mask.astype(np.uint8))
            expected.add("gaps_mask")
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
