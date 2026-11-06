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


def build_report_pages(
    report: ReportData, folder_name: str, png_dir: Path | None,
    report_images: ReportImagePaths,
) -> list[Figure]:
    """Build the velocity page and an optional-results page when applicable."""
    pages = [_velocity_page(report, folder_name, png_dir, report_images)]
    if report.spectral is not None or report.volume_rates:
        pages.append(_optional_results_page(report, folder_name))
    return pages


def _velocity_page(report: ReportData, folder_name: str, png_dir: Path | None,
                   images: ReportImagePaths):
    fig = plt.figure(figsize=(21 / 2.54, 29.7 / 2.54), dpi=150)
    fig.patch.set_facecolor("white")
    _add_page_header(fig, folder_name, "Velocity")
    for col, vessel in enumerate(("artery", "vein")):
        left = 0.08 + col * 0.48
        map_ax = fig.add_axes((left, 0.635, 0.39, 0.26))
        map_ax.axis("off")
        image = _try_load_vessel_image(png_dir, folder_name, vessel, images)
        if image is None:
            map_ax.text(0.5, 0.5, "RMS-frequency map unavailable",
                        ha="center", va="center", transform=map_ax.transAxes)
        else:
            _show_report_image(map_ax, image)
        map_ax.set_title(f"{vessel.title()} RMS-frequency map", fontsize=10)

        velocity_ax = fig.add_axes((left, 0.35, 0.39, 0.235))
        _plot_filtered_velocity(
            velocity_ax, report.velocity_per_beat.get(vessel),
            report.velocity_units.get(vessel, "mm/s"),
            report.parameters.get(f"{vessel[0].upper()}RI"),
        )
        velocity_ax.set_title(f"{vessel.title()} filtered per-beat velocity", fontsize=10)

    fig.text(
        0.08, 0.295,
        "Statistics use finite per-beat samples derived from filtered continuous "
        "traces (signs removed).\n"
        "RI and PI summarize the mean beat; plotted all-beat extrema are separate.",
        fontsize=8, va="top", linespacing=1.4,
    )
    params_ax = fig.add_axes((0.08, 0.045, 0.84, 0.19))
    params_ax.axis("off")
    _add_parameters_section(params_ax, report.parameters)
    return fig


def _plot_filtered_velocity(ax, values, unit: str, ri: dict | None) -> None:
    if values is None:
        ax.text(0.5, 0.5, "Filtered per-beat velocity unavailable",
                ha="center", va="center", transform=ax.transAxes, wrap=True)
        return
    finite = values[np.isfinite(values)]
    if not finite.size:
        ax.text(0.5, 0.5, "No finite filtered per-beat velocity samples",
                ha="center", va="center", transform=ax.transAxes, wrap=True)
        return
    phase = np.arange(values.shape[1], dtype=float) / values.shape[1]
    for beat in values:
        ax.plot(phase, beat, color="0.5", alpha=0.10, linewidth=0.6)
    mean_beat = np.mean(values, axis=0)
    ax.plot(phase, mean_beat, color="black", linewidth=1.5, label="Mean beat")
    for index, level in enumerate((np.min(finite), np.max(finite))):
        ax.axhline(level, color="0.4", linestyle=":", linewidth=0.7,
                   label="All-beat min/max" if index == 0 else None)
    if ri is not None and np.all(np.isfinite(mean_beat)):
        for index, level in enumerate((np.min(mean_beat), np.max(mean_beat))):
            ax.axhline(level, color="tab:blue", linestyle="--", linewidth=0.7,
                       label="Mean-beat RI extrema" if index == 0 else None)
    if ri is not None:
        ax.text(0.97, 0.97, f"Mean-beat RI = {ri['value']:.2f}",
                ha="right", va="top", transform=ax.transAxes, fontsize=8,
                backgroundcolor="white")
    ax.set_xlim(0, 1)
    ax.set_xlabel("Cardiac phase")
    ax.set_ylabel(f"Velocity ({unit})")
    ax.legend(loc="lower right", fontsize=7)


def _optional_results_page(report: ReportData, folder_name: str):
    fig = plt.figure(figsize=(21 / 2.54, 29.7 / 2.54), dpi=150)
    fig.patch.set_facecolor("white")
    sections = []
    if report.spectral is not None:
        sections.append(("Low-frequency gap spectra (shared across workflows)",
                         report.spectral, _draw_spectral_section))
    if report.volume_rates:
        sections.append(("Blood-volume rate (this workflow)",
                         report.volume_rates, _draw_volume_section))
    heading = " and ".join(
        title for title, present in (
            ("Low-frequency spectra", report.spectral is not None),
            ("blood-volume rate", bool(report.volume_rates)),
        ) if present
    )
    _add_page_header(fig, folder_name, heading)
    gs = fig.add_gridspec(len(sections), 2, left=0.08, right=0.94,
                          bottom=0.06, top=0.88, hspace=0.32, wspace=0.23)
    for row, (title, data, draw) in enumerate(sections):
        text_ax = fig.add_subplot(gs[row, 0])
        plot_ax = fig.add_subplot(gs[row, 1])
        text_ax.axis("off")
        text_ax.set_title(title, loc="left", fontsize=11, fontweight="bold")
        draw(text_ax, plot_ax, data)
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


def _draw_volume_section(text_ax, plot_ax,
                         rates: dict[str, VolumeRateReportData]) -> None:
    lines = ["Signed blood-volume rate (method shown per vessel):"]
    for vessel in ("artery", "vein"):
        data = rates.get(vessel)
        if data is None:
            continue
        value = "unavailable" if data.overall_value is None else f"{data.overall_value:.2f} {data.unit}"
        lines.extend((
            f"{vessel.title()}: {value}",
            f"  {data.method}; {data.beat_count} beats",
            f"  {data.valid_fraction * 100:.1f}% finite samples",
        ))
        phase = np.arange(data.phase_curve.size) / data.phase_curve.size
        plot_ax.plot(phase, data.phase_curve, label=vessel.title(), linewidth=1.2)
    lines.extend(("", "Total: beat median per phase, then phase median.",
                  "Gradient: local mean per phase, then phase mean",
                  "(not total flow). Coverage is descriptive, not pass/fail."))
    text_ax.text(0, 0.95, "\n".join(lines), va="top", fontsize=9,
                 linespacing=1.5, transform=text_ax.transAxes)
    plot_ax.set_xlim(0, 1)
    plot_ax.set_xlabel("Cardiac phase")
    plot_ax.set_ylabel("Flow rate (mm$^3$/s)")
    plot_ax.axhline(0, color="0.5", linestyle=":", linewidth=0.8)
    plot_ax.legend(fontsize=8)


def _percent(fraction: float | None) -> float | None:
    return None if fraction is None else 100 * fraction


def _add_page_header(fig: Figure, folder_name: str, section: str) -> None:
    """Identify both the acquisition and the subject of a report page."""
    fig.text(0.5, 0.965, folder_name, fontsize=16, fontweight="bold",
             ha="center", va="top")
    fig.text(0.5, 0.935, section, fontsize=11, ha="center", va="top")


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


def _add_parameters_section(ax, parameters: dict[str, Any]) -> None:
    """Lay out the formatted measurements in one or two columns."""
    param_texts = format_parameters_for_display(parameters)

    if not param_texts:
        ax.text(0.02, 0.5, "No parameters extracted from H5 files",
                fontsize=10, va='center')
        return

    ax.text(0.02, 0.98, 'Velocity statistics',
            fontsize=14, fontweight='bold', va='top')
    if len(param_texts) > 6:
        midpoint = len(param_texts) // 2
        columns = ((0.02, param_texts[:midpoint]), (0.52, param_texts[midpoint:]))
        y_step = 0.055
    else:
        columns = ((0.02, param_texts),)
        y_step = 0.08
    for x, texts in columns:
        for index, text in enumerate(texts):
            y = 0.86 - index * y_step
            if y > 0.05:
                ax.text(x, y, text, fontsize=9, va='top')


def _try_load_vessel_image(
    png_dir: Path | None, folder_name: str, vessel_type: str,
    report_images: ReportImagePaths,
):
    """Load only a real RMS-frequency map, never an unrelated M0 fallback."""
    possible_paths = []
    registered = report_images.get(("vessel_map", vessel_type))
    if registered is not None:
        possible_paths.append(Path(registered))
    if png_dir:
        possible_paths.append(png_dir / f"{folder_name}_{vessel_type}_seg_map_bkg.png")
    return load_first_image(possible_paths)


_FORMATS = {
    "UnixTimestampFirst": "First Unix Timestamp = %d",
    "UnixTimestampLast": "Last Unix Timestamp = %d",
    "heart_beat": "Spectral heart-rate estimate = %.1f",
    "Average_Arterial_Velocity": "Mean arterial velocity = %.2f",
    "Max_Arterial_Velocity": "Max arterial velocity = %.2f",
    "Min_Arterial_Velocity": "Min arterial velocity = %.2f",
    "Average_Venous_Velocity": "Mean venous velocity = %.2f",
    "Max_Venous_Velocity": "Max venous velocity = %.2f",
    "Min_Venous_Velocity": "Min venous velocity = %.2f",
    "ARI": "Arterial Resistivity Index = %.2f",
    "VRI": "Venous Resistivity Index = %.2f",
    "API": "Arterial Pulsatility Index = %.2f",
    "VPI": "Venous Pulsatility Index = %.2f",
}
_DISPLAY_ORDER = (
    "UnixTimestampFirst", "UnixTimestampLast", "heart_beat",
    "Average_Arterial_Velocity", "Max_Arterial_Velocity", "Min_Arterial_Velocity",
    "Average_Venous_Velocity", "Max_Venous_Velocity", "Min_Venous_Velocity",
    "ARI", "VRI", "API", "VPI", "SystoleDuration", "DiastoleDuration",
)


def format_parameters_for_display(parameters: dict[str, Any]) -> list[str]:
    """Format known parameters first, leaving other extracted values visible."""
    ordered = [name for name in _DISPLAY_ORDER if name in parameters]
    ordered.extend(name for name in parameters if name not in _DISPLAY_ORDER)
    return [_format_value(parameters[name], _FORMATS.get(name, "%s = %g")) for name in ordered]


def _format_value(value: dict | Any, fmt: str) -> str:
    if not isinstance(value, dict):
        return str(value)
    val = value.get("value", 0)
    unit = value.get("unit", "")
    try:
        if isinstance(val, (int, float)):
            formatted = fmt % val
        else:
            formatted = fmt.replace("%s", str(val)).replace("%g", str(val))
    except (TypeError, ValueError):
        formatted = f"{val}"
    return f"{formatted} {unit}".strip()
