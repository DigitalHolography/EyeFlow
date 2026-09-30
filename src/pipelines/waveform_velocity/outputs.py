"""Export waveform-velocity figure artifacts."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from calculations.math import nanmedian
from input_output.writers.eps import EpsArtifactWriter
from input_output.writers.png import PngArtifactWriter
from pipelines.retinal_velocity.semantics import (
    resolve_velocity_semantics,
    velocity_unit_from_payload,
)

VELOCITY_FIGURE_DPI = 320
VELOCITY_FIGURE_WIDTH_INCHES = 8.0 / 2.54
VELOCITY_FIGURE_ASPECT_RATIO = 1.618
VELOCITY_FONT_FAMILY = "Times New Roman"
VELOCITY_FONT_SIZE = 8
VELOCITY_ENVELOPE_GRAY = "0.8"
VELOCITY_BOX_LINE_WIDTH = 0.4
VELOCITY_SIGNAL_LINE_WIDTH = 0.8


def export_velocity_signals(
    output,
    artery_safe_segment_velocity,
    vein_safe_segment_velocity,
) -> list[Path]:
    """Export median safe segment velocity and its beatwise SD."""

    paths: list[Path] = []
    for vessel, values in (
        ("artery", artery_safe_segment_velocity),
        ("vein", vein_safe_segment_velocity),
    ):
        velocity = _metric_data(values)
        semantics = resolve_velocity_semantics(
            unit=velocity_unit_from_payload(values)
        )
        stem = f"velocity/{vessel}"
        paths.append(
            PngArtifactWriter(output, stem).save_figure(
                _velocity_figure(velocity, semantics.axis_label),
                "velocity.png",
                dpi=VELOCITY_FIGURE_DPI,
                bbox_inches=None,
            )
        )
        paths.append(
            EpsArtifactWriter(output, stem).save_figure(
                _velocity_figure(velocity, semantics.axis_label),
                "velocity.eps",
                dpi=VELOCITY_FIGURE_DPI,
            )
        )
    return paths


def _velocity_figure(
    safe_segment_velocity: np.ndarray,
    velocity_label: str = r"$v(t)$ (mm/s)",
):
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    from matplotlib.text import Text

    values = np.asarray(safe_segment_velocity, dtype=np.float64)
    if values.ndim != 4:
        raise ValueError(
            "safe segment velocity must have exactly four dimensions "
            "(time, beat, branch, radius)."
        )

    spatial_median = nanmedian(values, axis=(2, 3), dtype=np.float64)
    signal = nanmedian(spatial_median, axis=1, dtype=np.float64)
    standard_deviation = _nanstd(spatial_median, axis=1)
    cardiac_phase = np.linspace(0.0, 1.0, values.shape[0])

    figure_height = VELOCITY_FIGURE_WIDTH_INCHES / VELOCITY_FIGURE_ASPECT_RATIO
    fig = Figure(
        figsize=(VELOCITY_FIGURE_WIDTH_INCHES, figure_height),
        constrained_layout=True,
    )
    FigureCanvasAgg(fig)
    ax = fig.subplots()
    ax.fill_between(
        cardiac_phase,
        signal - standard_deviation,
        signal + standard_deviation,
        color=VELOCITY_ENVELOPE_GRAY,
        linewidth=0.0,
    )
    ax.axhline(
        0.0,
        color="black",
        linestyle=":",
        linewidth=VELOCITY_BOX_LINE_WIDTH,
    )
    ax.plot(
        cardiac_phase,
        signal,
        color="black",
        linestyle="-",
        linewidth=VELOCITY_SIGNAL_LINE_WIDTH,
    )

    ax.set_xlim(0.0, 1.0)
    ax.set_xlabel(r"Cardiac Phase $t/T$", fontsize=VELOCITY_FONT_SIZE)
    ax.set_ylabel(velocity_label, fontsize=VELOCITY_FONT_SIZE)
    ax.tick_params(axis="both", labelsize=VELOCITY_FONT_SIZE)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(VELOCITY_BOX_LINE_WIDTH)
    for text in fig.findobj(match=Text):
        text.set_fontfamily(VELOCITY_FONT_FAMILY)
    return fig


def _nanstd(values: np.ndarray, *, axis: int) -> np.ndarray:
    finite = np.isfinite(values)
    finite_count = np.sum(finite, axis=axis)
    mean = np.divide(
        np.sum(np.where(finite, values, 0.0), axis=axis),
        finite_count,
        out=np.full(values.shape[0], np.nan, dtype=np.float64),
        where=finite_count > 0,
    )
    squared_deviation = np.where(
        finite,
        (values - np.expand_dims(mean, axis=axis)) ** 2,
        0.0,
    )
    variance = np.divide(
        np.sum(squared_deviation, axis=axis),
        finite_count,
        out=np.full(values.shape[0], np.nan, dtype=np.float64),
        where=finite_count > 0,
    )
    return np.sqrt(variance)


def _metric_data(value) -> np.ndarray:
    if isinstance(value, tuple) and len(value) == 2 and isinstance(value[1], dict):
        value = value[0]
    return np.asarray(value)


__all__ = ["export_velocity_signals"]
