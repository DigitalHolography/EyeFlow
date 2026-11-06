"""One-page displacement and spatial-gradient PDF presentation."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure

from .width_data import PhaseSummary, WidthReportData


def build_width_page(report: WidthReportData, folder_name: str) -> Figure:
    """Show phase curves beside concise summaries of the same values."""
    fig = plt.figure(figsize=(21 / 2.54, 29.7 / 2.54), dpi=150)
    fig.patch.set_facecolor("white")
    fig.suptitle(f"{folder_name} — Motion and lumen width", fontsize=15, y=0.965)
    fig.text(0.5, 0.925, "Descriptive, pixel-based results; no clinical thresholds",
             ha="center", fontsize=9)
    grid = fig.add_gridspec(2, 2, left=0.08, right=0.94, bottom=0.11,
                            top=0.88, hspace=0.35, wspace=0.28)

    _draw_section(
        fig, grid[0, 0], grid[0, 1],
        title="Displacement: radial wall movement",
        series=report.displacement,
        ylabel="Radial movement (pixels)",
        empty_message="Displacement profiles not available",
    )
    _draw_section(
        fig, grid[1, 0], grid[1, 1],
        title="Spatial gradient: lumen width",
        series=report.gradient,
        ylabel="Lumen width (pixels)",
        empty_message="Spatial-gradient metrics not available",
    )
    fig.text(
        0.08, 0.055,
        "Lines show each phase's median over finite samples. The reported value is "
        "the median of that line.\n"
        "Displacement uses absolute local radial motion in wall bands relative to "
        "the temporal-mean image.\n"
        "Gradient width is masked transverse peak-to-peak spacing, reduced over "
        "beats and radii before this summary.",
        fontsize=7.3, va="bottom",
    )
    return fig


def _draw_section(
    fig: Figure, text_cell, plot_cell, *, title: str,
    series: tuple[PhaseSummary, ...], ylabel: str,
    empty_message: str,
) -> None:
    """Draw one analysis family as matching text and phase traces."""
    text_ax = fig.add_subplot(text_cell)
    plot_ax = fig.add_subplot(plot_cell)
    text_ax.axis("off")
    text_ax.set_title(title, loc="left", fontsize=10, fontweight="bold")
    if not series:
        text_ax.text(0, 0.91, empty_message, fontsize=9,
                     transform=text_ax.transAxes)
    else:
        lines = [_summary_line(item) for item in series]
        text_ax.text(0, 0.91, "\n\n".join(lines),
                     fontsize=8, va="top", linespacing=1.3,
                     transform=text_ax.transAxes)
        for item in series:
            if np.isfinite(item.curve).any():
                phase = np.arange(item.curve.size) / item.curve.size
                plot_ax.plot(phase, item.curve, label=item.label, linewidth=1.2)
    if plot_ax.lines:
        plot_ax.legend(fontsize=7, loc="best")
    else:
        plot_ax.text(0.5, 0.5, "No finite phase curve", ha="center",
                     va="center", transform=plot_ax.transAxes, fontsize=8)
    plot_ax.set_xlim(0, 1)
    plot_ax.set_xlabel("Cardiac phase")
    plot_ax.set_ylabel(ylabel)


def _summary_line(item: PhaseSummary) -> str:
    """Format the plotted median and its sample support."""
    value = "unavailable" if item.median is None else f"{item.median:.2f} {item.unit}"
    beats = f", {item.beat_count} beats" if item.beat_count is not None else ""
    return (
        f"{item.label}: {value}\n"
        f"{item.valid_fraction * 100:.1f}% finite; "
        f"{item.branch_count} branches{beats}"
    )
