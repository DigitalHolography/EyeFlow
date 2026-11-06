"""A4 page layout for one velocity workflow and its optional results."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure
from matplotlib.patches import FancyBboxPatch

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
        _add_report_image(
            fig, (left, 0.68, 0.39, 0.21), png_dir, folder_name,
            "vessel_map", vessel, images,
            f"{vessel.title()} RMS-frequency map",
        )
        velocity_ax = fig.add_axes((left, 0.47, 0.39, 0.18))
        _plot_filtered_velocity(
            velocity_ax,
            report.velocity_per_beat.get(vessel),
            report.velocity_units.get(vessel, "mm/s"),
            report.parameters.get(f"{vessel[0].upper()}RI"),
            f"{vessel.title()} filtered per-beat velocity",
        )
        _add_report_image(
            fig, (left, 0.235, 0.39, 0.18), png_dir, folder_name,
            "systole", vessel, images, f"{vessel.title()} systole indices",
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
    ax, values: np.ndarray | None, unit: str, ri: dict | None, title: str,
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
    for beat in samples:
        ax.plot(phase, beat, color="0.8", alpha=0.25, linewidth=0.6)
    counts = np.sum(np.isfinite(samples), axis=0)
    mean_beat = np.divide(
        np.nansum(samples, axis=0), counts,
        out=np.full(values.shape[1], np.nan), where=counts > 0,
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

    ri_value = _finite_parameter(ri)
    if ri_value is not None:
        ax.text(
            0.97, 0.97, f"Mean per-beat RI = {ri_value:.2f}",
            ha="right", va="top", transform=ax.transAxes, fontsize=7,
            backgroundcolor="white",
        )
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


def _optional_results_page(report: ReportData, folder_name: str):
    fig = plt.figure(figsize=(21 / 2.54, 29.7 / 2.54), dpi=150)
    fig.patch.set_facecolor("white")
    sections = []
    if report.spectral is not None:
        sections.append(("Low-Frequency Gap Spectra (Shared Across Workflows)",
                         report.spectral, _draw_spectral_section))
    if report.volume_rates:
        sections.append(("Blood-Volume Rate (This Workflow)",
                         report.volume_rates, _draw_volume_section))
    heading = " and ".join(
        title for title, present in (
            ("Low-Frequency Spectra", report.spectral is not None),
            ("Blood-Volume Rate", bool(report.volume_rates)),
        ) if present
    )
    _add_page_header(fig, folder_name, heading)
    gs = fig.add_gridspec(len(sections), 2, left=0.08, right=0.94,
                          bottom=0.06, top=0.88, hspace=0.32, wspace=0.23)
    for row, (title, data, draw) in enumerate(sections):
        text_ax = fig.add_subplot(gs[row, 0])
        plot_ax = fig.add_subplot(gs[row, 1])
        text_ax.axis("off")
        text_ax.set_title(title, loc="left", fontsize=13, fontweight="bold")
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
        plot_ax.plot(
            phase, data.phase_curve, label=vessel.title(),
            color=VESSEL_COLORS[vessel], linewidth=1.2,
        )
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
    ax.text(0, 0.98, "Hemodynamic summary", fontsize=12, fontweight="bold", va="top")
    ax.text(
        1, 0.98,
        f"Heart rate  |  {_parameter_value(parameters, 'heart_beat', decimals=1)}",
        fontsize=10, fontweight="bold", ha="right", va="top",
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "#f5f5f5",
              "edgecolor": "0.75", "linewidth": 0.8},
    )
    _draw_vessel_card(ax, parameters, "artery", 0.00)
    _draw_vessel_card(ax, parameters, "vein", 0.52)


def _draw_vessel_card(
    ax, parameters: dict[str, Any], vessel: str, left: float,
) -> None:
    color = VESSEL_COLORS[vessel]
    arterial = vessel == "artery"
    prefix = "Arterial" if arterial else "Venous"
    metric_prefix = "A" if arterial else "V"
    ax.add_patch(FancyBboxPatch(
        (left, 0.03), 0.48, 0.72,
        boxstyle="round,pad=0.012,rounding_size=0.025",
        facecolor="#fff7f7" if arterial else "#f5f8ff",
        edgecolor=color, linewidth=1.1,
    ))
    ax.text(
        left + 0.025, 0.68, "Artery" if arterial else "Vein",
        color=color, fontsize=10, fontweight="bold", va="center",
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
