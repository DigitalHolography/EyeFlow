"""Systole-index figure exporters."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from calculations.blood_flow_velocity.signal_analysis.waveform.cycles import (
    cycle_extrema as _cycle_extrema,
)
from plotting.article import (
    ARTICLE_ASPECT_RATIO,
    DIAGNOSTIC_ARTICLE_STYLE as STYLE,
)
from plotting.units import display_unit
from input_output.writers.png import FigureArtifactWriter as FigureWriter

from .common import (
    PulseFigureContext,
    _log,
    _plt,
    _vector,
    display_velocity as _display_velocity,
)
from .plotting import _style_axes


def _export_systole_plots(writer: FigureWriter, ctx: PulseFigureContext) -> list[Path]:
    peaks = ctx.cycle_boundary_indexes
    if peaks.size == 0:
        _log(ctx, "Skipping systole index PNGs; no beat indices are available.")
        return []
    artery = _display_velocity(_vector(ctx.velocity.continuous("artery")))
    vein = _display_velocity(_vector(ctx.velocity.continuous("vein")))
    maxima, minima = _cycle_extrema(artery, peaks)
    artifact_groups = [
        _systole_plot(
            writer,
            "find_systoles_indices_artery.png",
            ctx.time,
            artery,
            peaks,
            maxima,
            minima,
            ctx.velocity_axis_label,
            "Artery systole indices",
        ),
        _systole_plot(
            writer,
            "find_systoles_indices_vein.png",
            ctx.time,
            vein,
            peaks,
            maxima,
            minima,
            ctx.velocity_axis_label,
            "Vein systole indices",
        ),
    ]
    for vessel, artifact_paths in zip(("artery", "vein"), artifact_groups):
        writer.register_artifact(
            "systole", vessel,
            next(path for path in artifact_paths if path.suffix.lower() == ".png"),
        )
    return [path for artifact_paths in artifact_groups for path in artifact_paths]

def _systole_plot(
    writer: FigureWriter,
    suffix: str,
    time: np.ndarray,
    filtered: np.ndarray,
    peaks: np.ndarray,
    maxima: np.ndarray,
    minima: np.ndarray,
    velocity_label: str = f"Velocity ({display_unit('mm/s')})",
    title: str = "Systole indices",
) -> list[Path]:
    fig, ax = _plt().subplots(figsize=(7.5, 7.5 / ARTICLE_ASPECT_RATIO))
    ax.plot(
        time[: filtered.size], filtered, color="black",
        linewidth=STYLE.secondary_line_width,
    )
    for indexes, marker_color in (
        (peaks, "k"), (maxima, "tab:red"), (minima, "tab:blue")
    ):
        valid = indexes[(indexes >= 0) & (indexes < filtered.size)]
        if valid.size:
            ax.scatter(
                time[valid], filtered[valid], color=marker_color, s=28, zorder=3
            )
            for idx in valid:
                ax.axvline(
                    time[idx], color=marker_color, linestyle="--",
                    linewidth=STYLE.event_line_width, alpha=0.75,
                )
    ax.set_xlabel("Time (s)")
    ax.set_ylabel(velocity_label)
    _style_axes(ax)
    ax.set_title(
        title, fontsize=10, fontfamily="sans-serif", fontweight="normal"
    )
    return writer.save_figure_pair(fig, suffix)
