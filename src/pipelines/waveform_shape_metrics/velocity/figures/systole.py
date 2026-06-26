"""Systole-index PNG exporters."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from calculations.blood_flow_velocity.signal_analysis.waveform import (
    cycle_extrema as _cycle_extrema,
)
from input_output.writers.png import PngArtifactWriter as FigureWriter

from .common import (
    PulseFigureContext,
    _log,
    _plt,
    _safe_indexes,
    _vector,
    display_velocity as _display_velocity,
)
from .plotting import _style_axes


def _export_systole_plots(writer: FigureWriter, ctx: PulseFigureContext) -> list[Path]:
    peaks = _safe_indexes(ctx.analysis.get("beat_indices"))
    if peaks.size == 0:
        _log(ctx, "Skipping systole index PNGs; no beat indices are available.")
        return []
    artery = _display_velocity(
        _vector(ctx.analysis.get("retinal_artery_velocity_signal_filtered"))
    )
    artery_deriv = _display_velocity(
        _vector(ctx.analysis.get("retinal_artery_velocity_signal_derivative"))
    )
    vein = _display_velocity(_vector(ctx.analysis.get("retinal_vein_velocity_signal_filtered")))
    vein_deriv = _display_velocity(
        _vector(ctx.analysis.get("retinal_vein_velocity_signal_derivative"))
    )
    maxima, minima = _cycle_extrema(artery, peaks)
    return [
        _systole_plot(
            writer,
            "find_systoles_indices_artery.png",
            ctx.time,
            artery,
            artery_deriv,
            peaks,
            maxima,
            minima,
        ),
        _systole_plot(
            writer,
            "find_systoles_indices_vein.png",
            ctx.time,
            vein,
            vein_deriv,
            peaks,
            maxima,
            minima,
        ),
    ]

def _systole_plot(
    writer: FigureWriter,
    suffix: str,
    time: np.ndarray,
    filtered: np.ndarray,
    derivative: np.ndarray,
    peaks: np.ndarray,
    maxima: np.ndarray,
    minima: np.ndarray,
) -> Path:
    fig, ax = _plt().subplots(figsize=(7.5, 4.0))
    ax.plot(time[: derivative.size], derivative, color="0.55", linewidth=1.5)
    ax.plot(time[: filtered.size], filtered, color="k", linewidth=1.5)
    for indexes, color in ((peaks, "k"), (maxima, "tab:red"), (minima, "tab:blue")):
        valid = indexes[(indexes >= 0) & (indexes < filtered.size)]
        if valid.size:
            ax.scatter(time[valid], filtered[valid], color=color, s=28, zorder=3)
            for idx in valid:
                ax.axvline(time[idx], color=color, linestyle="--", linewidth=1.0, alpha=0.75)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Velocity (mm/s)")
    _style_axes(ax)
    return writer.savefig(fig, suffix)
