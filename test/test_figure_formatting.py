"""General figure controls preserve mappings, themes, and scientific artists."""

import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from matplotlib.cm import ScalarMappable
import numpy as np
import pytest

from input_output.figure_formatting import (
    ColorbarLayout, LegendLayout, TickFormat,
    add_colorbar, add_legend, apply_tick_format, format_tick_labels,
)


@pytest.fixture(autouse=True)
def headless_figures():
    with mpl.rc_context({"backend": "Agg"}):
        yield
    plt.close("all")


@pytest.mark.parametrize("style, expected", [
    (TickFormat(decimals=0), ["0", "14", "-18", "0"]),
    (TickFormat(decimals=1), ["0.0", "14.1", "-18.3", "0.0"]),
    (TickFormat(significant_digits=3), ["0", "14.1", "-18.3", "-0.001"]),
])
def test_mapped_tick_precision_uses_physical_values_and_preserves_bin_order(style, expected):
    fig, ax = plt.subplots()
    ax.plot([0, 7], [1, 2])
    ax.set_xlim(0, 7)
    apply_tick_format(ax.xaxis, style, positions=[0, 2, 5, 7], values=[0, 14.1, -18.3, -0.001])
    assert [label.get_text() for label in ax.get_xticklabels()] == expected
    np.testing.assert_array_equal(ax.get_xticks(), [0, 2, 5, 7])
    assert ax.get_xlim() == (0, 7)
    np.testing.assert_array_equal(ax.lines[0].get_ydata(), [1, 2])


def test_default_format_preserves_phase_fractions_and_logarithmic_locators():
    fig, (phase, spectrum) = plt.subplots(1, 2)
    phase.set_xticks([0, 1 / 3, 2 / 3, 1], labels=["0", "1/3", "2/3", "1"])
    spectrum.set_yscale("log")
    locators = [phase.xaxis.get_major_locator(), spectrum.yaxis.get_major_locator()]
    formatters = [phase.xaxis.get_major_formatter(), spectrum.yaxis.get_major_formatter()]
    for axis in (phase.xaxis, spectrum.yaxis):
        apply_tick_format(axis)
    assert phase.xaxis.get_major_locator() is locators[0]
    assert spectrum.yaxis.get_major_locator() is locators[1]
    assert phase.xaxis.get_major_formatter() is formatters[0]
    assert spectrum.yaxis.get_major_formatter() is formatters[1]
    assert [label.get_text() for label in phase.get_xticklabels()] == ["0", "1/3", "2/3", "1"]
    assert spectrum.get_yscale() == "log"


def test_explicit_ticks_do_not_expand_limits_or_change_autoscaling():
    fig, ax = plt.subplots()
    ax.plot([0, 1], [2, 4])
    limits = ax.get_xlim()
    assert ax.get_autoscalex_on()
    apply_tick_format(ax.xaxis, positions=[-20, 0, 1, 20])
    assert ax.get_xlim() == limits
    assert ax.get_autoscalex_on()


def test_numeric_tick_density_precision_and_major_minor_appearance():
    fig, ax = plt.subplots()
    ax.plot([0, 20], [-1, 1])
    ax.minorticks_on()
    limits = ax.get_xlim()
    apply_tick_format(ax.xaxis, TickFormat(
        nbins=4, integer=True, decimals=0, direction="in", length=6,
        minor_length=2, width=0.7, padding=5,
    ))
    fig.canvas.draw()
    assert ax.get_xlim() == limits
    assert all(value == round(value) for value in ax.get_xticks())
    assert ax.xaxis.get_major_formatter()(-0.001) == "0"
    assert all(tick.tick1line.get_markersize() == 6 for tick in ax.xaxis.get_major_ticks())
    assert all(tick.tick1line.get_markersize() == 2 for tick in ax.xaxis.get_minor_ticks())
    assert all(tick.get_pad() == 5 for tick in ax.xaxis.get_major_ticks())


def test_scientific_notation_and_offset_are_opt_in():
    fig, ax = plt.subplots()
    ax.plot([1e6, 2e6], [1, 2])
    limits = ax.get_xlim()
    apply_tick_format(ax.xaxis, TickFormat(nbins=4, scientific=True,
                                         scientific_limits=(0, 0), use_offset=False))
    fig.canvas.draw()
    assert "6" in ax.xaxis.get_offset_text().get_text()
    assert not ax.xaxis.get_major_formatter().get_useOffset()
    assert ax.get_xlim() == limits
    apply_tick_format(ax.xaxis, TickFormat(scientific=False, use_offset=False))
    fig.canvas.draw()
    assert ax.xaxis.get_offset_text().get_text() == ""


@pytest.mark.parametrize("options", [
    {"decimals": -1}, {"decimals": 1, "significant_digits": 3},
    {"decimals": 1, "scientific": True}, {"nbins": 0},
])
def test_conflicting_or_invalid_tick_options_fail_early(options):
    with pytest.raises(ValueError):
        TickFormat(**options)


def test_automatic_density_cannot_replace_mapped_or_logarithmic_ticks():
    fig, ax = plt.subplots()
    with pytest.raises(ValueError, match="explicit tick positions"):
        apply_tick_format(ax.xaxis, TickFormat(nbins=4), positions=[0, 1])
    with pytest.raises(ValueError, match="same length"):
        apply_tick_format(ax.xaxis, positions=[0, 1], values=[1])
    ax.set_xscale("log")
    before = ax.xaxis.get_major_locator()
    with pytest.raises(ValueError, match="linear axis"):
        apply_tick_format(ax.xaxis, TickFormat(nbins=4))
    assert ax.xaxis.get_major_locator() is before


def test_general_legend_layout_preserves_theme_labels_and_signal_styles():
    with mpl.rc_context({"font.family": "monospace", "legend.fontsize": 13,
                         "legend.frameon": False}):
        fig, ax = plt.subplots()
        first, = ax.plot([0, 1], [1, 2], "r--", label="early")
        second, = ax.plot([0, 1], [2, 3], "b:", label="late")
        legend = add_legend(ax, LegendLayout(
            loc="upper left", bbox_to_anchor=(1.02, 1), ncols=2, labelspacing=0.3,
            columnspacing=1.2, handlelength=1.5, frameon=False,
        ), frameon=True)
        fig.canvas.draw()
        assert [text.get_text() for text in legend.get_texts()] == ["early", "late"]
        assert all(text.get_fontfamily() == ["monospace"] for text in legend.get_texts())
        assert all(text.get_fontsize() == 13 for text in legend.get_texts())
        assert legend.get_frame_on()  # Explicit overrides win, including boolean values.
        assert first.get_linestyle() == "--" and second.get_color() == "b"
        np.testing.assert_allclose(
            legend.get_bbox_to_anchor().transformed(ax.transAxes.inverted()).bounds[:2], [1.02, 1],
        )
        figure_legend = add_legend(fig, LegendLayout(loc="upper center", ncols=2),
                                   handles=[first, second])
        assert len(figure_legend.get_texts()) == 2


def test_colorbar_with_explicit_axes_keeps_geometry_and_tick_format_on_update():
    fig, ax = plt.subplots()
    ax.plot([0, 1], [1, 2])
    cax = fig.add_axes((0.2, 0.1, 0.6, 0.05))
    main_bounds, colorbar_bounds = ax.get_position().bounds, cax.get_position().bounds
    mappable = ScalarMappable(norm=Normalize(0, 2), cmap="viridis")
    norm, cmap = mappable.norm, mappable.cmap
    colorbar = add_colorbar(
        fig, mappable, cax=cax, label="Velocity (mm/s)",
        layout=ColorbarLayout(orientation="horizontal", fraction=0.03, pad=0.02,
                              label_padding=7, label_size=11, tick_size=9, outline_visible=False),
        ticks=TickFormat(decimals=2),
    )
    assert colorbar.orientation == "horizontal"
    assert cax.get_xlabel() == "Velocity (mm/s)" and cax.xaxis.labelpad == 7
    assert cax.xaxis.label.get_fontsize() == 11
    assert all(tick.get_fontsize() == 9 for tick in cax.get_xticklabels())
    assert not colorbar.outline.get_visible()
    np.testing.assert_allclose(ax.get_position().bounds, main_bounds)
    np.testing.assert_allclose(cax.get_position().bounds, colorbar_bounds)
    assert mappable.norm is norm and mappable.cmap is cmap
    mappable.set_clim(0, 3)
    colorbar.update_ticks()
    assert colorbar.formatter(0.25) == "0.25"
    assert colorbar.vmin == 0 and colorbar.vmax == 3


def test_attached_colorbar_preserves_log_normalization_and_theme():
    with mpl.rc_context({"font.family": "monospace", "axes.labelsize": 14}):
        fig, ax = plt.subplots()
        norm = LogNorm(1, 100)
        image = ax.imshow([[1, 10], [10, 100]], norm=norm, cmap="magma")
        colorbar = add_colorbar(fig, image, ax=ax, label="Power",
                               layout=ColorbarLayout(pad=0.04, fraction=0.05, shrink=0.8))
        fig.canvas.draw()
        assert colorbar.mappable is image and colorbar.norm is norm
        assert colorbar.ax.get_yscale() == "log"
        assert colorbar.ax.yaxis.label.get_fontfamily() == ["monospace"]
        assert colorbar.ax.yaxis.label.get_fontsize() == 14
        assert image.get_clim() == (1, 100)
