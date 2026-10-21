"""Export velocity-analysis figure artifacts."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from calculations.math import nanmedian
from input_output.article_formatting import (
    ARTICLE_ASPECT_RATIO as VELOCITY_FIGURE_ASPECT_RATIO,
    ARTICLE_DPI as VELOCITY_FIGURE_DPI,
    ARTICLE_WIDTH_INCHES as VELOCITY_FIGURE_WIDTH_INCHES,
    COMPACT_ARTICLE_STYLE,
    apply_article_style,
    create_article_figure,
    prepare_article_fill_layer,
)
from input_output.writers.article import save_article_figure_pair
from pipelines.velocity.semantics import (
    resolve_velocity_semantics,
    velocity_unit_from_payload,
)

VELOCITY_ENVELOPE_GRAY = COMPACT_ARTICLE_STYLE.envelope_color


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
        paths.extend(save_article_figure_pair(
            output, _velocity_figure(velocity, semantics.axis_label), "velocity",
            stem=stem, dpi=VELOCITY_FIGURE_DPI,
        ))
    return paths


def _velocity_figure(
    safe_segment_velocity: np.ndarray,
    velocity_label: str = r"$v(t)$ (mm/s)",
):
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

    fig, ax = create_article_figure(
        width=VELOCITY_FIGURE_WIDTH_INCHES, aspect_ratio=VELOCITY_FIGURE_ASPECT_RATIO,
    )
    if COMPACT_ARTICLE_STYLE.envelope_alpha < 1:
        prepare_article_fill_layer(ax)
    ax.fill_between(
        cardiac_phase,
        signal - standard_deviation,
        signal + standard_deviation,
        **COMPACT_ARTICLE_STYLE.envelope_kwargs(),
        zorder=0,
    )
    ax.axhline(
        0.0,
        color="black",
        linestyle=":",
        linewidth=COMPACT_ARTICLE_STYLE.reference_line_width,
    )
    ax.plot(
        cardiac_phase,
        signal,
        color="black",
        linestyle="-",
        linewidth=COMPACT_ARTICLE_STYLE.signal_line_width,
    )

    ax.set_xlim(0.0, 1.0)
    ax.set_xlabel(r"Cardiac Phase $t/T$")
    ax.set_ylabel(velocity_label)
    apply_article_style(fig, COMPACT_ARTICLE_STYLE)
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
