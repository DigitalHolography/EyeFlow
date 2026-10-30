"""Opt-in tick, legend, and colorbar controls for any Matplotlib figure.

No article fonts, export formats, backend selection, or global rcParams are
imposed. Unspecified options preserve Matplotlib defaults or existing ticks.
"""

from __future__ import annotations

from dataclasses import dataclass, fields


@dataclass(frozen=True)
class TickFormat:
    """Numeric precision/density and tick appearance; None preserves settings.

    ``nbins`` is the maximum number of major intervals, as in MaxNLocator.
    Explicit formatting is opt-in so phase fractions and logarithmic axes keep
    their own locators/formatters by default. Dimensions are in points.
    """

    decimals: int | None = None
    significant_digits: int | None = None
    nbins: int | None = None
    integer: bool = False
    scientific: bool | None = None
    scientific_limits: tuple[int, int] | None = None
    use_offset: bool | None = None
    direction: str | None = None
    length: float | None = None
    width: float | None = None
    padding: float | None = None
    minor_length: float | None = None

    def __post_init__(self):
        for name, minimum in (("decimals", 0), ("significant_digits", 1), ("nbins", 1)):
            value = getattr(self, name)
            if value is not None and (not isinstance(value, int) or value < minimum):
                raise ValueError(f"{name} must be an integer >= {minimum}.")
        if self.decimals is not None and self.significant_digits is not None:
            raise ValueError("Choose decimal precision or significant digits, not both.")
        if self.has_precision and (self.scientific is True or self.scientific_limits is not None
                                   or self.use_offset is not None):
            raise ValueError("Fixed precision cannot be combined with scientific/offset settings.")

    @property
    def has_precision(self) -> bool:
        return self.decimals is not None or self.significant_digits is not None


def format_tick_labels(values, style: TickFormat = TickFormat()) -> list[str]:
    """Format supplied physical values, including mapped bins; suppress -0."""
    if style.scientific is True or style.scientific_limits is not None or style.use_offset is not None:
        raise ValueError("Mapped labels use precision settings, not scientific offsets.")
    labels = []
    for value in values:
        value = float(value)
        if style.decimals is not None:
            label = f"{value:.{style.decimals}f}"
        elif style.significant_digits is not None:
            label = f"{value:.{style.significant_digits}g}"
        else:
            label = str(value)
        labels.append(label.lstrip("-") if float(label) == 0 else label)
    return labels


def apply_tick_format(axis, style: TickFormat = TickFormat(), *, positions=None, values=None) -> None:
    """Configure one X/Y Axis without changing its scale, limits, or data.

    Supply positions and physical values for bin-index axes. Their labels are
    formatted from those values rather than the index. A default TickFormat
    preserves custom labels and locators. Density options replace only linear
    locators; explicit positions and automatic density are mutually exclusive.
    """
    from matplotlib.ticker import FuncFormatter, MaxNLocator, ScalarFormatter

    automatic = style.nbins is not None or style.integer
    if values is not None and positions is None:
        raise ValueError("Mapped tick values require positions.")
    if automatic and positions is not None:
        raise ValueError("Choose explicit tick positions or automatic density, not both.")
    if automatic and axis.get_scale() != "linear":
        raise ValueError("Automatic numeric tick density requires a linear axis.")
    if positions is not None:
        positions = list(positions)
        labels = None if values is None else format_tick_labels(values, style)
        if labels is not None and len(labels) != len(positions):
            raise ValueError("Tick positions and mapped values must have the same length.")
        ax = axis.axes
        horizontal = axis.axis_name == "x"
        limits = ax.get_xlim() if horizontal else ax.get_ylim()
        autoscale = ax.get_autoscalex_on() if horizontal else ax.get_autoscaley_on()
        axis.set_ticks(positions, labels=labels)
        if horizontal:
            ax.set_xlim(limits, auto=autoscale)
        else:
            ax.set_ylim(limits, auto=autoscale)
    elif automatic:
        options = {"integer": style.integer}
        if style.nbins is not None:
            options["nbins"] = style.nbins
        axis.set_major_locator(MaxNLocator(**options))
    if values is None:
        if style.has_precision:
            axis.set_major_formatter(FuncFormatter(
                lambda value, position: format_tick_labels([value], style)[0],
            ))
        elif (style.scientific is not None or style.scientific_limits is not None
              or style.use_offset is not None):
            formatter = ScalarFormatter()
            if style.scientific is not None:
                formatter.set_scientific(style.scientific)
            if style.scientific_limits is not None:
                formatter.set_powerlimits(style.scientific_limits)
            if style.use_offset is not None:
                formatter.set_useOffset(style.use_offset)
            axis.set_major_formatter(formatter)
    options = {key: value for key, value in (
        ("direction", style.direction), ("width", style.width), ("pad", style.padding),
    ) if value is not None}
    if options:
        axis.axes.tick_params(axis=axis.axis_name, which="both", **options)
    if style.length is not None:
        axis.axes.tick_params(axis=axis.axis_name, which="major", length=style.length)
    if style.minor_length is not None:
        axis.axes.tick_params(axis=axis.axis_name, which="minor", length=style.minor_length)


@dataclass(frozen=True)
class LegendLayout:
    """Reusable legend placement/spacing; None delegates to Matplotlib."""

    loc: str | int | None = None
    ncols: int | None = None
    frameon: bool | None = None
    bbox_to_anchor: tuple[float, ...] | None = None
    fontsize: float | str | None = None
    borderpad: float | None = None
    labelspacing: float | None = None
    columnspacing: float | None = None
    handlelength: float | None = None
    handletextpad: float | None = None
    borderaxespad: float | None = None


def _configured(layout) -> dict:
    return {field.name: getattr(layout, field.name) for field in fields(layout)
            if getattr(layout, field.name) is not None}


def add_legend(ax, layout: LegendLayout = LegendLayout(), *, handles=None, labels=None, **overrides):
    """Create a legend without imposing article styling or changing handles.

    Explicit keyword overrides win over the layout. Pass a Figure instead of
    Axes to create a figure-level legend using its supported locations.
    """
    options = _configured(layout)
    if handles is not None:
        options["handles"] = handles
    if labels is not None:
        options["labels"] = labels
    return ax.legend(**(options | overrides))


@dataclass(frozen=True)
class ColorbarLayout:
    """Colorbar geometry plus optional label/tick/outline appearance.

    pad/fraction/shrink/anchors follow Matplotlib's relative geometry; label
    padding/font sizes/outline width use points. An explicit cax retains its
    supplied geometry. No colormap, normalization, or units are chosen here.
    """

    orientation: str | None = None
    location: str | None = None
    fraction: float | None = None
    pad: float | None = None
    shrink: float | None = None
    aspect: float | None = None
    anchor: tuple[float, float] | None = None
    panchor: tuple[float, float] | bool | None = None
    label_padding: float | None = None
    label_size: float | None = None
    tick_size: float | None = None
    outline_visible: bool | None = None
    outline_width: float | None = None


def add_colorbar(
    fig, mappable, *, ax=None, cax=None, label: str | None = None,
    layout: ColorbarLayout = ColorbarLayout(), ticks: TickFormat | None = None, **overrides,
):
    """Add an attached or standalone colorbar; retain calibration and colors."""
    options = _configured(layout) | overrides
    appearance = {name: options.pop(name) for name in (
        "label_padding", "label_size", "tick_size", "outline_visible", "outline_width",
    ) if name in options}
    colorbar = fig.colorbar(mappable, ax=ax, cax=cax, **options)
    label_options = {}
    if "label_padding" in appearance:
        label_options["labelpad"] = appearance["label_padding"]
    if "label_size" in appearance:
        label_options["fontsize"] = appearance["label_size"]
    if label is not None:
        colorbar.set_label(label, **label_options)
    if "tick_size" in appearance:
        colorbar.ax.tick_params(labelsize=appearance["tick_size"])
    if "outline_visible" in appearance:
        colorbar.outline.set_visible(appearance["outline_visible"])
    if "outline_width" in appearance:
        colorbar.outline.set_linewidth(appearance["outline_width"])
    if ticks is not None:
        axis = colorbar.ax.xaxis if colorbar.orientation == "horizontal" else colorbar.ax.yaxis
        apply_tick_format(axis, ticks)
        colorbar.locator = axis.get_major_locator()
        colorbar.formatter = axis.get_major_formatter()
    return colorbar
