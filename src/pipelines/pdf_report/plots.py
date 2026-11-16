"""Reusable plots shared by the moments PDF and HTML presentations."""

from __future__ import annotations

import numpy as np

from plotting.article import ARTICLE_ASPECT_RATIO, prepare_article_fill_layer
from plotting.colors import opaque_on_background
from plotting.units import display_unit

from .data import SystoleReportData

GRAPH_TITLE_STYLE = {
    "fontsize": 10,
    "fontfamily": "sans-serif",
    "fontweight": "normal",
}


def plot_systole_indices(
    ax, data: SystoleReportData | None, title: str,
) -> None:
    """Render systole diagnostics at the presentation's final dimensions."""
    ax.set_title(title, **GRAPH_TITLE_STYLE)
    if data is None:
        ax.axis("off")
        ax.text(
            0.5, 0.5, "Systole indices unavailable",
            ha="center", va="center", transform=ax.transAxes,
        )
        return

    values = np.asarray(data.velocity, dtype=float)
    if data.sample_interval_seconds is None:
        horizontal = np.arange(values.size, dtype=float)
        xlabel = "Frame"
    else:
        horizontal = np.arange(values.size, dtype=float) * data.sample_interval_seconds
        xlabel = "Time (s)"

    ax.plot(horizontal, values, color="black", linewidth=1.2, zorder=2)
    _draw_systole_markers(
        ax, horizontal, values, data.systolic_boundaries,
        marker="v", size=18, color="black", label="Systolic boundary",
    )
    _draw_systole_markers(
        ax, horizontal, values, data.maxima,
        marker="o", size=14, color="tab:red", label="Systolic maximum",
    )
    _draw_systole_markers(
        ax, horizontal, values, data.minima,
        marker="o", size=16, color="tab:blue", label="Diastolic minimum",
    )
    ax.set_xlabel(xlabel, fontsize=8)
    ax.set_ylabel(f"Velocity ({display_unit(data.unit)})", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.set_box_aspect(1 / ARTICLE_ASPECT_RATIO)
    ax.margins(x=0.01)


def _draw_systole_markers(
    ax, horizontal, values, indexes, *, marker, size, color, label,
) -> None:
    valid = _valid_signal_indexes(indexes, values.size)
    if not valid.size:
        return
    ax.scatter(
        horizontal[valid], values[valid], marker=marker, s=size,
        facecolors=color, edgecolors=color, linewidths=0.7,
        label=label, zorder=3,
    )
    for index in valid:
        ax.axvline(
            horizontal[index],
            color=opaque_on_background(color, opacity=0.65),
            linestyle=":", linewidth=0.7, zorder=0,
        )


def _valid_signal_indexes(indexes, size: int) -> np.ndarray:
    indexes = np.asarray(indexes, dtype=np.int64).reshape(-1)
    return indexes[(indexes >= 0) & (indexes < size)]


def plot_filtered_velocity(
    ax, values: np.ndarray | None, unit: str, title: str, *,
    unfiltered_values: np.ndarray | None = None,
) -> None:
    """Plot filtered beats over a robust view of the unfiltered samples."""
    ax.set_title(title, **GRAPH_TITLE_STYLE)
    if values is None:
        ax.axis("off")
        ax.text(0.5, 0.5, "Filtered per-beat velocity unavailable",
                ha="center", va="center", transform=ax.transAxes, wrap=True)
        return

    values = np.asarray(values, dtype=float)
    finite = values[np.isfinite(values)]
    if values.ndim != 2 or not finite.size:
        ax.axis("off")
        ax.text(0.5, 0.5, "No finite filtered per-beat velocity samples",
                ha="center", va="center", transform=ax.transAxes, wrap=True)
        return

    phase = np.arange(values.shape[1], dtype=float) / values.shape[1]
    samples = np.where(np.isfinite(values), values, np.nan)
    counts = np.sum(np.isfinite(samples), axis=0)
    mean_beat = np.divide(
        np.nansum(samples, axis=0), counts,
        out=np.full(values.shape[1], np.nan), where=counts > 0,
    )
    prepare_article_fill_layer(ax)
    _plot_unfiltered_velocity(ax, unfiltered_values)
    ax.plot(phase, mean_beat, color="black", linewidth=1.5, label="Mean beat")

    ax.set_xlim(0, 1)
    ax.set_xlabel("Cardiac phase", fontsize=8)
    ax.set_ylabel(f"Velocity ({display_unit(unit)})", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.set_box_aspect(1 / ARTICLE_ASPECT_RATIO)


def _plot_unfiltered_velocity(ax, values: np.ndarray | None) -> None:
    """Add a robust raw-signal envelope and mark phase-wise IQR outliers."""
    if values is None:
        return
    samples = np.asarray(values, dtype=float)
    if samples.ndim != 2 or not all(samples.shape):
        return
    if not np.any(np.isfinite(samples)):
        return
    phase = np.arange(samples.shape[1], dtype=float) / samples.shape[1]
    quantiles = _phase_quantiles(samples, (10.0, 25.0, 50.0, 75.0, 90.0))
    lower, first_quartile, median, third_quartile, upper = quantiles
    envelope_valid = np.isfinite(lower) & np.isfinite(upper)
    if np.any(envelope_valid):
        ax.fill_between(
            phase, lower, upper, where=envelope_valid,
            color="black", alpha=0.18, linewidth=0,
            label="Unfiltered 10th-90th percentile envelope", zorder=-2,
        )
    ax.plot(
        phase, median, color="black", linestyle=(0, (3, 2)),
        linewidth=0.9, label="Unfiltered median", zorder=1,
    )

    interquartile_range = third_quartile - first_quartile
    low_fence = first_quartile - 1.5 * interquartile_range
    high_fence = third_quartile + 1.5 * interquartile_range
    finite = np.isfinite(samples)
    outlier = finite & (
        (samples < low_fence[None, :]) | (samples > high_fence[None, :])
    )
    beat_indexes, phase_indexes = np.nonzero(outlier)
    if beat_indexes.size:
        ax.scatter(
            phase[phase_indexes], samples[beat_indexes, phase_indexes],
            s=10, facecolors="none", edgecolors="black", linewidths=0.55,
            label="Unfiltered IQR outlier", zorder=2,
        )


def _phase_quantiles(samples: np.ndarray, quantiles) -> np.ndarray:
    """Compute phase-wise finite quantiles without all-NaN warnings."""
    output = np.full((len(quantiles), samples.shape[1]), np.nan, dtype=float)
    for phase_index in range(samples.shape[1]):
        finite = samples[:, phase_index]
        finite = finite[np.isfinite(finite)]
        if finite.size:
            output[:, phase_index] = np.percentile(finite, quantiles)
    return output


__all__ = [
    "GRAPH_TITLE_STYLE",
    "plot_filtered_velocity",
    "plot_systole_indices",
]
