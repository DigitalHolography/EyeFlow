"""A4 page layout for one velocity workflow and its optional results."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure

from input_output.reports.pdf_report import load_first_image
from input_output.writers.png import ReportImagePaths

from .data import (
    ReportData,
    SpectralReportData,
    VolumeRateReportData,
)
from .palette import VESSEL_COLORS

GRAPH_TITLE_STYLE = {
    "fontsize": 10,
    "fontfamily": "sans-serif",
    "fontweight": "normal",
}

VOLUME_RATE_METHODS = (
    ("gradient_edges", "Spatial-Gradient Edge (Static)", "local"),
    ("masked_edges", "Mask-Derived Geometry", "total"),
)


def build_report_pages(
    report: ReportData, folder_name: str, png_dir: Path | None,
    report_images: ReportImagePaths,
) -> list[Figure]:
    """Build ordered velocity, BVR, and spectral report pages."""
    pages = [_velocity_page(report, folder_name, png_dir, report_images)]
    if report.volume_rates:
        pages.append(_volume_rate_page(report.volume_rates, folder_name))
    if report.spectral is not None:
        pages.append(_spectral_page(report.spectral, folder_name))
    return pages


def _velocity_page(report: ReportData, folder_name: str, png_dir: Path | None,
                   images: ReportImagePaths):
    fig = plt.figure(figsize=(21 / 2.54, 29.7 / 2.54), dpi=150)
    fig.patch.set_facecolor("white")
    _add_page_header(fig, folder_name, "Velocity")
    for col, vessel in enumerate(("artery", "vein")):
        left = 0.08 + col * 0.48
        _add_report_image(
            fig, (left, 0.68, 0.39, 0.21), png_dir, folder_name,
            "vessel_map", vessel, images,
            f"{vessel.title()} RMS-frequency map",
        )
        _add_report_image(
            fig, (left, 0.47, 0.39, 0.18), png_dir, folder_name,
            "systole", vessel, images, f"{vessel.title()} systole indices",
        )
        velocity_ax = fig.add_axes((left, 0.27, 0.39, 0.17))
        _plot_filtered_velocity(
            velocity_ax,
            report.velocity_per_beat.get(vessel),
            report.velocity_units.get(vessel, "mm/s"),
            f"{vessel.title()} filtered per-beat velocity",
        )

    params_ax = fig.add_axes((0.08, 0.025, 0.84, 0.19))
    params_ax.axis("off")
    _add_parameters_section(params_ax, report.parameters)
    return fig


def _add_report_image(
    fig: Figure,
    bounds: tuple[float, float, float, float],
    png_dir: Path | None,
    folder_name: str,
    kind: str,
    vessel: str,
    report_images: ReportImagePaths,
    title: str,
) -> None:
    """Place one existing artifact without rebuilding its calculation."""
    ax = fig.add_axes(bounds)
    ax.axis("off")
    image = _try_load_report_image(
        png_dir, folder_name, kind, vessel, report_images
    )
    if image is None:
        ax.text(0.5, 0.5, f"{title} unavailable", ha="center", va="center",
                transform=ax.transAxes, wrap=True)
    else:
        if kind == "systole":
            image = _without_embedded_title(image)
        _show_report_image(ax, image)
    if kind in {"vessel_map", "systole"}:
        ax.set_title(title, **GRAPH_TITLE_STYLE)


def _plot_filtered_velocity(
    ax, values: np.ndarray | None, unit: str, title: str,
) -> None:
    """Plot the exact filtered per-beat samples used by report statistics."""
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
    _plot_beat_whiskers(
        ax, phase, samples.T, mean_beat, color="0.45", label="Beat min/max",
    )
    ax.plot(phase, mean_beat, color="black", linewidth=1.5, label="Mean beat")

    for index, level in enumerate((np.min(finite), np.max(finite))):
        ax.axhline(level, color="0.5", linestyle=":", linewidth=0.8,
                   label="All-beat min/max" if index == 0 else None)
    finite_mean = mean_beat[np.isfinite(mean_beat)]
    if finite_mean.size:
        for index, level in enumerate((np.min(finite_mean), np.max(finite_mean))):
            ax.axhline(level, color="0.25", linestyle="--", linewidth=0.8,
                       label="Mean-beat extrema" if index == 0 else None)

    ax.set_xlim(0, 1)
    ax.set_xlabel("Cardiac phase", fontsize=8)
    ax.set_ylabel(f"Velocity ({unit})", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.legend(loc="lower right", fontsize=6)


def _finite_parameter(parameter: dict | None) -> float | None:
    if not isinstance(parameter, dict):
        return None
    try:
        value = float(parameter["value"])
    except (KeyError, TypeError, ValueError):
        return None
    return value if np.isfinite(value) else None


def _volume_rate_page(
    methods: dict[str, dict[str, VolumeRateReportData]], folder_name: str,
) -> Figure:
    fig = plt.figure(figsize=(21 / 2.54, 29.7 / 2.54), dpi=150)
    fig.patch.set_facecolor("white")
    _add_page_header(fig, folder_name, "Blood-Volume Rate")
    selected = [method for method in VOLUME_RATE_METHODS if method[0] in methods]
    layouts = (
        ((0.54, 0.27, 0.84), (0.13, 0.27, 0.43))
        if len(selected) == 2 else ((0.49, 0.31, 0.83),)
    )
    for (key, subtitle, scope), (bottom, height, title_y) in zip(selected, layouts):
        fig.text(0.08, title_y, subtitle, fontsize=13, fontweight="bold", va="top")
        text_ax = fig.add_axes((0.08, bottom, 0.32, height))
        plot_ax = fig.add_axes((0.48, bottom, 0.44, height))
        text_ax.axis("off")
        _draw_volume_method(text_ax, plot_ax, methods[key], scope)
    return fig


def _spectral_page(data: SpectralReportData, folder_name: str) -> Figure:
    fig = plt.figure(figsize=(21 / 2.54, 29.7 / 2.54), dpi=150)
    fig.patch.set_facecolor("white")
    _add_page_header(fig, folder_name, "Low-Frequency Spectra")
    fig.text(
        0.08, 0.84, "Low-Frequency Gap Spectra (Shared Across Workflows)",
        fontsize=13, fontweight="bold", va="top",
    )
    text_ax = fig.add_axes((0.08, 0.43, 0.38, 0.36))
    plot_ax = fig.add_axes((0.54, 0.43, 0.38, 0.36))
    text_ax.axis("off")
    _draw_spectral_section(text_ax, plot_ax, data)
    return fig


def _draw_spectral_section(text_ax, plot_ax, data: SpectralReportData) -> None:
    def number(value, suffix=""):
        return "unavailable" if value is None else f"{value:.2f}{suffix}"

    lines = [
        f"Median spectral excess (SEI): {number(data.sei)}",
        f"Early / late SEI: {number(data.sei_early)} / {number(data.sei_late)}",
        "",
        "Descriptive quality indicators (no pass/fail threshold):",
        f"Finite acquisition spectra: {number(_percent(data.acquisition_valid_fraction), '%')}",
        f"Valid cardiac phases for SE: {number(_percent(data.valid_se_phase_fraction), '%')}",
        f"Retained gap-mask field: {number(_percent(data.retained_field_fraction), '%')}",
        "",
        "SEI is the median of the phase-resolved SE curve.",
        "Mask coverage is not an anatomical quality score.",
    ]
    text_ax.text(0, 0.95, "\n".join(lines), va="top", fontsize=9,
                 linespacing=1.5, transform=text_ax.transAxes)
    if data.se_curve is None or not data.se_curve.size:
        plot_ax.text(0.5, 0.5, "Phase-resolved SE unavailable",
                     ha="center", va="center", transform=plot_ax.transAxes)
    else:
        phase = np.arange(data.se_curve.size) / data.se_curve.size
        plot_ax.plot(phase, data.se_curve, color="black", linewidth=1.2)
        if data.sei is not None:
            plot_ax.axhline(data.sei, color="tab:blue", linestyle="--",
                            linewidth=0.9, label="Median SEI")
            plot_ax.legend(fontsize=8)
    plot_ax.set_xlim(0, 1)
    plot_ax.set_xlabel("Cardiac phase")
    plot_ax.set_ylabel("Spectral excess SE (unitless)")


def _draw_volume_method(
    text_ax, plot_ax, rates: dict[str, VolumeRateReportData], scope: str,
) -> None:
    lines = []
    unit = "mm^3/s"
    for vessel in ("artery", "vein"):
        data = rates.get(vessel)
        if data is None:
            continue
        unit = data.unit
        finite = data.phase_curve[np.isfinite(data.phase_curve)]
        minimum = float(np.min(finite)) if finite.size else None
        maximum = float(np.max(finite)) if finite.size else None
        pulsatile_range = None if not finite.size else maximum - minimum
        lines.extend((
            vessel.title(),
            f"  Median {scope} BVR: {_rate_value(data.overall_value, data.unit)}",
            f"  Cycle minimum: {_rate_value(minimum, data.unit)}",
            f"  Cycle maximum: {_rate_value(maximum, data.unit)}",
            f"  Pulsatile range: {_rate_value(pulsatile_range, data.unit)}",
            "",
        ))
        phase = np.arange(data.phase_curve.size) / data.phase_curve.size
        _plot_beat_whiskers(
            plot_ax, phase, data.beat_curves, data.phase_curve,
            color=VESSEL_COLORS[vessel],
        )
        plot_ax.plot(
            phase, data.phase_curve, label=f"{vessel.title()} median",
            color=VESSEL_COLORS[vessel], linewidth=1.6,
        )
    text_ax.text(0, 0.98, "\n".join(lines), va="top", fontsize=8.5,
                 linespacing=1.35, transform=text_ax.transAxes)
    plot_ax.set_title(
        f"Median {scope} BVR with beat min/max whiskers", **GRAPH_TITLE_STYLE,
    )
    plot_ax.set_xlim(0, 1)
    plot_ax.set_xlabel("Cardiac phase")
    plot_ax.set_ylabel(f"BVR ({unit})")
    plot_ax.axhline(0, color="0.5", linestyle=":", linewidth=0.8)
    plot_ax.legend(fontsize=8)


def _rate_value(value: float | None, unit: str) -> str:
    return "Unavailable" if value is None else f"{value:.2f} {unit}"


def _plot_beat_whiskers(
    ax,
    phase: np.ndarray,
    beat_curves: np.ndarray,
    center: np.ndarray,
    *,
    color: str,
    label: str | None = None,
    maximum_count: int = 24,
) -> None:
    """Draw sparse capped min/max ranges around a phase-resolved center."""
    values = np.asarray(beat_curves, dtype=float)
    center = np.asarray(center, dtype=float)
    if values.ndim != 2 or values.shape[0] != center.size:
        return
    masked = np.ma.masked_invalid(values)
    minimum = np.asarray(masked.min(axis=1).filled(np.nan), dtype=float)
    maximum = np.asarray(masked.max(axis=1).filled(np.nan), dtype=float)
    step = max(1, int(np.ceil(center.size / maximum_count)))
    indexes = np.arange(0, center.size, step)
    valid = (
        np.isfinite(center[indexes])
        & np.isfinite(minimum[indexes])
        & np.isfinite(maximum[indexes])
    )
    indexes = indexes[valid]
    if not indexes.size:
        return
    ax.errorbar(
        phase[indexes], center[indexes],
        yerr=np.vstack((
            center[indexes] - minimum[indexes],
            maximum[indexes] - center[indexes],
        )),
        fmt="none", ecolor=color, alpha=0.32,
        elinewidth=0.7, capsize=1.8, capthick=0.7,
        label=label, zorder=1,
    )


def _percent(fraction: float | None) -> float | None:
    return None if fraction is None else 100 * fraction


def _add_page_header(fig: Figure, folder_name: str, section: str) -> None:
    """Identify both the acquisition and the subject of a report page."""
    section = section[:1].upper() + section[1:]
    fig.text(0.5, 0.978, folder_name, fontsize=10, color="0.35",
             ha="center", va="top")
    fig.text(0.5, 0.952, section, fontsize=16, fontweight="bold",
             ha="center", va="top")


def _show_report_image(ax, image: np.ndarray) -> None:
    """Draw an image without distorting its pixel aspect ratio."""
    height, width = np.asarray(image).shape[:2]
    ax.imshow(image, aspect="equal")
    ax.set_anchor("C")
    ax.set_adjustable("box")
    ax.margins(0)

    if width > 0 and height > 0:
        zoom = 1.04
        x_center = (width - 1) / 2
        y_center = (height - 1) / 2
        half_width = width / (2 * zoom)
        half_height = height / (2 * zoom)
        ax.set_xlim(x_center - half_width, x_center + half_width)
        ax.set_ylim(y_center + half_height, y_center - half_height)


def _without_embedded_title(image: np.ndarray) -> np.ndarray:
    """Remove the controlled PNG title margin before adding the PDF title."""
    image = np.asarray(image)
    crop = int(round(image.shape[0] * 0.06))
    return image[crop:] if 0 < crop < image.shape[0] else image


def _add_parameters_section(ax, parameters: dict[str, Any]) -> None:
    """Show paired vessel measurements and one shared cardiac value."""
    ax.text(0, 0.98, "Summary", fontsize=12, fontweight="bold", va="top")
    ax.text(
        1, 0.98,
        f"Heart rate  |  {_parameter_value(parameters, 'heart_beat', decimals=1)}",
        fontsize=10, fontweight="bold", ha="right", va="top",
    )
    _draw_vessel_column(ax, parameters, "artery", 0.00)
    _draw_vessel_column(ax, parameters, "vein", 0.52)


def _draw_vessel_column(
    ax, parameters: dict[str, Any], vessel: str, left: float,
) -> None:
    arterial = vessel == "artery"
    prefix = "Arterial" if arterial else "Venous"
    metric_prefix = "A" if arterial else "V"
    ax.text(
        left + 0.025, 0.68, "Artery" if arterial else "Vein",
        color="black", fontsize=10, fontweight="bold", va="center",
    )
    rows = (
        ("Minimum velocity", f"Min_{prefix}_Velocity"),
        ("Mean velocity", f"Average_{prefix}_Velocity"),
        ("Maximum velocity", f"Max_{prefix}_Velocity"),
        ("Resistivity index", f"{metric_prefix}RI"),
        ("Pulsatility index", f"{metric_prefix}PI"),
    )
    for row, (label, key) in enumerate(rows):
        y = 0.56 - row * 0.115
        ax.text(left + 0.025, y, label, fontsize=8.2, color="0.25", va="center")
        ax.text(
            left + 0.455, y, _parameter_value(parameters, key),
            fontsize=8.5, fontweight="bold", ha="right", va="center",
        )


def _parameter_value(
    parameters: dict[str, Any], name: str, *, decimals: int = 2,
) -> str:
    parameter = parameters.get(name)
    value = _finite_parameter(parameter)
    if value is None:
        return "Unavailable"
    unit = str(parameter.get("unit", "")).strip()
    return f"{value:.{decimals}f}{f' {unit}' if unit else ''}"


def _try_load_report_image(
    png_dir: Path | None, folder_name: str, kind: str, vessel_type: str,
    report_images: ReportImagePaths,
):
    """Prefer the current variant's registered artifact, then its legacy path."""
    possible_paths = []
    registered = report_images.get((kind, vessel_type))
    if registered is not None:
        possible_paths.append(Path(registered))
    if png_dir:
        suffixes = {
            "vessel_map": f"{vessel_type}_seg_map_bkg.png",
            "systole": f"find_systoles_indices_{vessel_type}.png",
        }
        possible_paths.append(png_dir / f"{folder_name}_{suffixes[kind]}")
    return load_first_image(possible_paths)
