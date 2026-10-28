"""Shared, opt-in Matplotlib formatting for article figure artifacts.

Build the plot first, then apply a style before saving PNG/EPS. Writers remain
format-neutral; scientific labels, limits, colors, and overlays belong to callers.
Matplotlib imports stay lazy and no process-wide rcParams are changed.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from math import isfinite

from .formatting import ColorbarLayout, LegendLayout, TickFormat, apply_tick_format

ARTICLE_DPI = 320
ARTICLE_WIDTH_INCHES = 8.0 / 2.54
ARTICLE_ASPECT_RATIO = 1.618


@dataclass(frozen=True)
class ArticleStyle:
    """Text/line sizes in points, fill colors/opacity, and optional margins.

    ``None`` preserves a plot's own label/tick padding and data margins.
    Pass line widths and fill keyword dictionaries when creating the artists;
    applying text/border styling never guesses which artist is a signal, a
    reference, a map overlay, or an envelope.
    """

    font_family: str = "Times New Roman"
    math_font_family: str = "stix"
    label_size: float = 10
    tick_size: float = 8
    legend_size: float = 10
    title_size: float = 10
    annotation_size: float = 9
    border_width: float = 0.4
    border_color: str = "black"
    signal_line_width: float = 0.8
    secondary_line_width: float = 1.5
    reference_line_width: float = 0.4
    event_line_width: float = 1.0
    fit_line_width: float = 0.625
    envelope_color: str = "0.8"
    envelope_alpha: float = 1.0
    region_color: str = "0.5"
    region_alpha: float = 0.1
    histogram_color: str = "0.8"
    histogram_edge_color: str = "0.25"
    histogram_edge_width: float = 0.25
    x_tick_rotation: float = 0
    y_tick_rotation: float = 0
    label_padding: float | None = None
    tick_padding: float | None = None
    x_margin: float | None = None
    y_margin: float | None = None
    x_ticks: TickFormat = TickFormat()
    y_ticks: TickFormat = TickFormat()
    legend_layout: LegendLayout = LegendLayout()
    colorbar_layout: ColorbarLayout = ColorbarLayout()

    def envelope_kwargs(self) -> dict:
        """Solid dispersion fills; overlap opacity is chosen by the preset."""
        return {"facecolor": self.envelope_color, "alpha": self.envelope_alpha,
                "edgecolor": "none", "linewidth": 0.0}

    def region_kwargs(self) -> dict:
        """A softer background for an analysis region, distinct from envelopes."""
        return {"facecolor": self.region_color, "alpha": self.region_alpha,
                "edgecolor": "none"}

    def histogram_kwargs(self) -> dict:
        return {"color": self.histogram_color, "edgecolor": self.histogram_edge_color,
                "linewidth": self.histogram_edge_width}

    def scaled_fonts(self, factor: float) -> ArticleStyle:
        """Match apparent text size across canvases displayed at equal width."""
        if not isfinite(factor) or factor <= 0:
            raise ValueError("Font scale must be positive and finite.")
        return replace(self, **{
            name: getattr(self, name) * factor
            for name in ("label_size", "tick_size", "legend_size", "title_size", "annotation_size")
        })


ARTICLE_STYLE = ArticleStyle()
COMPACT_ARTICLE_STYLE = replace(
    ARTICLE_STYLE, label_size=8, legend_size=8, title_size=8, annotation_size=8,
    legend_layout=LegendLayout(loc="upper right", frameon=False, handlelength=2),
)
SPECTRUM_ARTICLE_STYLE = replace(
    ARTICLE_STYLE, envelope_color="0.35", envelope_alpha=0.45,
    y_ticks=TickFormat(decimals=0, nbins=5, integer=True),
    legend_layout=LegendLayout(loc="upper center", frameon=False, labelspacing=0.25, borderpad=0.1),
)
DIAGNOSTIC_ARTICLE_STYLE = replace(
    ARTICLE_STYLE, signal_line_width=2, reference_line_width=1.5,
    legend_layout=LegendLayout(loc="best"),
    colorbar_layout=ColorbarLayout(fraction=0.046, pad=0.04),
)
SPECTROGRAM_ARTICLE_STYLE = replace(
    ARTICLE_STYLE, colorbar_layout=ColorbarLayout(orientation="vertical", pad=0.02),
)


def prepare_article_fill_layer(ax) -> None:
    """Composite translucent fills for EPS, leaving signals/text vector.

    Call when using transparent fills. Put backgrounds/regions/envelopes below
    zorder 1 and signals/references above it. Repeated calls share one opaque
    background, so overlapping envelopes blend before EPS rasterization.
    """
    from matplotlib.patches import Rectangle

    gid = "article_fill_background"
    if not any(patch.get_gid() == gid for patch in ax.patches):
        ax.add_patch(Rectangle(
            (0, 0), 1, 1, transform=ax.transAxes, facecolor=ax.get_facecolor(),
            edgecolor="none", zorder=-2, gid=gid,
        ))
    ax.set_rasterization_zorder(1)


def create_article_figure(
    *, width: float = ARTICLE_WIDTH_INCHES,
    aspect_ratio: float = ARTICLE_ASPECT_RATIO,
    square: bool = False,
    padding_inches: float | None = None,
    signal_padding: tuple[bool, bool] = (True, True),
):
    """Create an Agg figure and axes with fixed export dimensions.

    Normal figures use constrained layout; ``padding_inches`` controls its
    outer padding. Square figures are full bleed, for masks and image maps.
    ``signal_padding`` is (horizontal, vertical): False removes the automatic
    data margin on that axis; True keeps Matplotlib's default margin. Explicit
    axis limits or ArticleStyle margins applied later take precedence.
    Existing multi-panel figures can use :func:`apply_article_style` directly.
    """
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    if not isfinite(width) or width <= 0 or not isfinite(aspect_ratio) or aspect_ratio <= 0:
        raise ValueError("Figure width and aspect ratio must be positive and finite.")
    if padding_inches is not None and (not isfinite(padding_inches) or padding_inches < 0):
        raise ValueError("Figure padding must be nonnegative and finite.")
    if square and padding_inches not in (None, 0):
        raise ValueError("Full-bleed square figures cannot have outer padding.")
    if (not isinstance(signal_padding, tuple) or len(signal_padding) != 2
            or any(not isinstance(value, bool) for value in signal_padding)):
        raise ValueError("signal_padding must be a (horizontal, vertical) tuple of booleans.")
    height = width / (1.0 if square else aspect_ratio)
    fig = Figure(figsize=(width, height), constrained_layout=not square)
    FigureCanvasAgg(fig)
    if not square and padding_inches is not None:
        fig.get_layout_engine().set(w_pad=padding_inches, h_pad=padding_inches)
    ax = fig.add_axes((0, 0, 1, 1)) if square else fig.subplots()
    margins = {axis: 0 for axis, enabled in zip(("x", "y"), signal_padding) if not enabled}
    if margins:
        ax.margins(**margins)
    return fig, ax


def apply_article_style(fig, style: ArticleStyle = ARTICLE_STYLE, *, axes=None) -> None:
    """Format an existing figure, or selected axes, without resizing it.

    Apply after labels, legends, annotations, and colorbars are created. Colorbar
    axes are included by default. No data, units, limits, or line styles change
    unless the caller explicitly requests data margins in the style.
    """
    from matplotlib.text import Text

    selected = list(fig.axes if axes is None else axes)
    containers = [fig] if axes is None else selected
    for container in containers:
        for text in container.findobj(match=Text):
            text.set_fontfamily(style.font_family)
            text.set_math_fontfamily(style.math_font_family)
            text.set_fontsize(style.annotation_size)
    for ax in selected:
        tick_options = {"pad": style.tick_padding} if style.tick_padding is not None else {}
        ax.tick_params(axis="both", which="both", labelsize=style.tick_size,
                       width=style.border_width, **tick_options)
        ax.tick_params(axis="x", labelrotation=style.x_tick_rotation)
        ax.tick_params(axis="y", labelrotation=style.y_tick_rotation)
        for axis in (ax.xaxis, ax.yaxis):
            axis.label.set_fontsize(style.label_size)
            axis.get_offset_text().set_fontsize(style.tick_size)
            if style.label_padding is not None:
                axis.labelpad = style.label_padding
        ax.title.set_fontsize(style.title_size)
        for spine in ax.spines.values():
            spine.set_color(style.border_color)
            spine.set_linewidth(style.border_width)
        legend = ax.get_legend()
        if legend is not None:
            for text in [*legend.get_texts(), legend.get_title()]:
                text.set_fontsize(style.legend_size)
        margins = {name: value for name, value in (("x", style.x_margin), ("y", style.y_margin))
                   if value is not None}
        if margins:
            ax.margins(**margins)
        apply_tick_format(ax.xaxis, style.x_ticks)
        apply_tick_format(ax.yaxis, style.y_ticks)
    if axes is None:
        for legend in fig.legends:
            for text in [*legend.get_texts(), legend.get_title()]:
                text.set_fontsize(style.legend_size)
        for text in fig.texts:
            text.set_fontsize(style.title_size)
