"""Shared article formatting preserves scientific artists and plot geometry."""

from dataclasses import replace
from io import BytesIO
from pathlib import Path

import matplotlib as mpl
import numpy as np
import pytest

from input_output.article_formatting import (
    ARTICLE_STYLE,
    apply_article_style,
    create_article_figure,
    prepare_article_fill_layer,
)
from pipelines.velocity_analysis.artifacts.figures.plotting import (
    _colorbar,
    _image_map,
    _line_plot,
)


def test_article_style_formats_text_regions_without_changing_scientific_artists():
    settings = {key: mpl.rcParams[key] for key in (
        "font.family", "font.size", "mathtext.fontset", "axes.linewidth", "axes.xmargin",
    )}
    fig, ax = create_article_figure()
    line, = ax.plot([0, 1, 2], [1, np.nan, 3], "r--", label="signal")
    image = ax.imshow([[1, 2], [3, 4]], extent=(0, 2, 0, 4), aspect="auto")
    colorbar = fig.colorbar(image, ax=ax, label="Power")
    ax.set(xlabel=r"Phase ($t/T$)", ylabel=r"$A_1$", title="Metric")
    annotation = ax.text(0.5, 0.5, "peak")
    legend = ax.legend()
    fig.suptitle("Recording")
    ax.set(xlim=(0, 2), ylim=(-1, 5))
    before = (fig.get_size_inches().copy(), ax.get_xlim(), ax.get_ylim(), image.get_clim())
    style = replace(ARTICLE_STYLE, label_size=12, tick_size=9, legend_size=11,
                    title_size=13, annotation_size=8, label_padding=7, tick_padding=3)

    apply_article_style(fig, style)
    fig.canvas.draw()

    np.testing.assert_array_equal(line.get_ydata(), [1, np.nan, 3])
    assert line.get_color() == "r" and line.get_linestyle() == "--"
    np.testing.assert_array_equal(fig.get_size_inches(), before[0])
    assert (ax.get_xlim(), ax.get_ylim(), image.get_clim()) == before[1:]
    assert ax.xaxis.label.get_fontfamily() == ["Times New Roman"]
    assert ax.yaxis.label.get_math_fontfamily() == "stix"
    assert ax.xaxis.label.get_fontsize() == colorbar.ax.yaxis.label.get_fontsize() == 12
    assert ax.title.get_fontsize() == 13 and annotation.get_fontsize() == 8
    assert all(text.get_fontsize() == 11 for text in legend.get_texts())
    assert all(text.get_fontsize() == 9 for text in colorbar.ax.get_yticklabels())
    assert ax.xaxis.labelpad == 7
    assert all(mpl.rcParams[key] == value for key, value in settings.items())


def test_selected_axes_preserve_sibling_styling_and_hidden_borders():
    fig, ax = create_article_figure()
    sibling = ax.twinx()
    sibling.set_ylabel("Acceleration", fontsize=17, fontfamily="DejaVu Sans")
    ax.spines["top"].set_visible(False)
    sibling_size = sibling.yaxis.label.get_fontsize()
    apply_article_style(fig, axes=[ax])
    assert not ax.spines["top"].get_visible()
    assert sibling.yaxis.label.get_fontsize() == sibling_size
    assert sibling.yaxis.label.get_fontfamily() == ["DejaVu Sans"]


def test_explicit_layout_padding_data_margins_and_full_bleed():
    fig, ax = create_article_figure(width=6, aspect_ratio=2, padding_inches=0.12)
    ax.plot([0, 1], [1, 2])
    apply_article_style(fig, replace(ARTICLE_STYLE, x_margin=0, y_margin=0.1))
    fig.canvas.draw()
    np.testing.assert_allclose(fig.get_size_inches(), [6, 3])
    np.testing.assert_allclose(ax.get_xlim(), [0, 1])
    np.testing.assert_allclose(ax.get_ylim(), [0.9, 2.1])
    assert fig.get_layout_engine().get()["w_pad"] == 0.12
    assert fig.get_layout_engine().get()["h_pad"] == 0.12
    square, mask_ax = create_article_figure(width=2, square=True)
    mask_ax.set_axis_off()
    apply_article_style(square)
    assert not mask_ax.axison and square.get_layout_engine() is None
    np.testing.assert_allclose(mask_ax.get_position().bounds, [0, 0, 1, 1])
    np.testing.assert_allclose(square.get_size_inches(), [2, 2])
    with pytest.raises(ValueError, match="padding"):
        create_article_figure(square=True, padding_inches=0.1)


@pytest.mark.parametrize("kind", ["signal", "map", "colorbar"])
def test_velocity_diagnostic_exporters_apply_the_shared_article_style(kind):
    class CaptureWriter:
        def savefig(self, fig, suffix, **kwargs):
            self.figure = fig
            fig.canvas.draw()
            return Path(suffix)

    writer = CaptureWriter()
    if kind == "signal":
        _line_plot(writer, "signal.png", np.array([0, 1]),
                   [(np.array([1, 2]), "--", "red", "artery")],
                   xlabel="Time (s)", ylabel="Velocity (mm/s)")
    elif kind == "map":
        _image_map(writer, np.ones((2, 2)), "map.png", cmap="viridis",
                   colorbar=True, label="Velocity (mm/s)")
        assert not writer.figure.axes[0].axison
        assert writer.figure.axes[0].get_aspect() == 1
    else:
        _colorbar(writer, "colorbar.png", cmap="viridis", vmin=0, vmax=1,
                  label="Velocity (mm/s)")
    for ax in writer.figure.axes:
        for text in (ax.xaxis.label, ax.yaxis.label, *ax.get_xticklabels()):
            assert text.get_fontfamily() == [ARTICLE_STYLE.font_family]
        assert ax.xaxis.label.get_fontsize() == ARTICLE_STYLE.label_size
    from matplotlib import pyplot as plt
    plt.close(writer.figure)


@pytest.mark.parametrize("horizontal, vertical", [(True, True), (True, False),
                                                  (False, True), (False, False)])
def test_signal_padding_controls_data_margins_independently(horizontal, vertical):
    with mpl.rc_context({"axes.xmargin": 0.1, "axes.ymargin": 0.2}):
        fig, ax = create_article_figure(
            padding_inches=0.12, signal_padding=(horizontal, vertical),
        )
        ax.plot([2, 8], [10, 30])
        apply_article_style(fig)
        fig.canvas.draw()
        np.testing.assert_allclose(ax.get_xlim(), [1.4, 8.6] if horizontal else [2, 8])
        np.testing.assert_allclose(ax.get_ylim(), [6, 34] if vertical else [10, 30])
        assert fig.get_layout_engine().get()["w_pad"] == 0.12
        assert fig.get_layout_engine().get()["h_pad"] == 0.12


@pytest.mark.parametrize("invalid", [False, (True,), (True, "False")])
def test_signal_padding_requires_two_boolean_values(invalid):
    with pytest.raises(ValueError, match="tuple of booleans"):
        create_article_figure(signal_padding=invalid)


@pytest.mark.parametrize("kind", ["spectrum", "velocity", "flow"])
def test_line_and_shading_overrides_reach_consumers_and_composite_for_eps(kind, monkeypatch, caplog):
    overrides = dict(signal_line_width=1.7, reference_line_width=0.6,
                     envelope_color="0.2", envelope_alpha=0.3,
                     region_color="0.7", region_alpha=0.12)
    style = replace(ARTICLE_STYLE, **overrides)
    beat_values = np.array([[1, 2], [2, 4], [3, 4]], dtype=float)
    if kind == "spectrum":
        from pipelines.low_freq_spectra import artifacts as plots
        style = replace(plots.SPECTRUM_STYLE, **overrides)
        monkeypatch.setattr(plots, "SPECTRUM_STYLE", style)
        curve = plots._spectrum_curve(
            np.array([1, 4, 3, 2, 2, 3, 4, 1]), np.full(8, 0.1),
            frequencies_khz=np.fft.fftfreq(8), label="Median", linestyle="--",
            index_max1=1, index_max2=6,
        )
        fig = plots._spectra_figure([curve])
        expected = curve.power.copy()
        signal, reference = fig.axes[0].lines
        region = next(patch for patch in fig.axes[0].patches if patch.get_alpha() is not None)
        assert region.get_alpha() == style.region_alpha
    elif kind == "velocity":
        from pipelines.velocity_analysis.artifacts import velocity_signals as plots
        monkeypatch.setattr(plots, "COMPACT_ARTICLE_STYLE", style)
        fig = plots._velocity_figure(beat_values[:, :, None, None])
        expected = np.median(beat_values, axis=1)
        reference, signal = fig.axes[0].lines
    else:
        from pipelines.blood_volume_rate import outputs as plots
        monkeypatch.setattr(plots, "COMPACT_ARTICLE_STYLE", style)
        fig = plots._blood_volume_rate_figure(beat_values)
        expected = np.median(beat_values, axis=1)
        reference, signal = fig.axes[0].lines
    ax = fig.axes[0]
    assert signal.get_linewidth() == 1.7 and reference.get_linewidth() == 0.6
    np.testing.assert_allclose(signal.get_ydata(), expected)
    np.testing.assert_allclose(ax.collections[0].get_facecolor()[0], [0.2, 0.2, 0.2, 0.3])
    assert not ax.collections[0].get_edgecolors().size
    prepare_article_fill_layer(ax)  # Multiple envelopes share one composite background.
    assert sum(patch.get_gid() == "article_fill_background" for patch in ax.patches) == 1
    assert ax.get_rasterization_zorder() == 1
    assert ax.collections[0].get_zorder() < 1 <= signal.get_zorder()
    eps = BytesIO()
    with caplog.at_level("WARNING", logger="matplotlib.backends.backend_ps"):
        fig.savefig(eps, format="eps", dpi=120)
    assert b"colorimage" in eps.getvalue()
    assert "does not support transparency" not in caplog.text


def test_histogram_and_fit_use_shared_settings_without_changing_bin_counts(monkeypatch):
    from pipelines.blood_volume_rate import outputs as plots
    style = replace(ARTICLE_STYLE, histogram_color="0.6", histogram_edge_color="0.1",
                    histogram_edge_width=0.5, fit_line_width=1.2)
    monkeypatch.setattr(plots, "COMPACT_ARTICLE_STYLE", style)
    ax = plots._lumen_diameter_distribution_figure(np.array([10, 20, 30, 40, 50])).axes[0]
    assert sum(patch.get_height() for patch in ax.patches) == 5
    for patch in ax.patches:
        np.testing.assert_allclose(patch.get_facecolor()[:3], [0.6] * 3)
        np.testing.assert_allclose(patch.get_edgecolor()[:3], [0.1] * 3)
        assert patch.get_linewidth() == 0.5
    assert ax.lines[0].get_linewidth() == 1.2
