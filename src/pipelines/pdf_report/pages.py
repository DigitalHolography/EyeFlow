"""A4 page layout for the two-page moments report."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure

from input_output.reports.pdf_report import load_first_image
from input_output.writers.artifact_names import prefixed_artifact_path
from input_output.writers.png import ReportImagePaths
from plotting.article import ARTICLE_ASPECT_RATIO, prepare_article_fill_layer
from plotting.units import display_unit

from .data import (
    ReportData,
    VolumeRateReportData,
)

GRAPH_TITLE_STYLE = {
    "fontsize": 10,
    "fontfamily": "sans-serif",
    "fontweight": "normal",
}
VELOCITY_GRAPH_WIDTH = 0.39
VELOCITY_GRAPH_HEIGHT = (
    VELOCITY_GRAPH_WIDTH * 21 / (29.7 * ARTICLE_ASPECT_RATIO)
)

SECOND_PAGE_ARTIFACTS = {
    "histogram": lambda vessel: Path(f"histogramVelocity{vessel}.png"),
    "branch_identity": lambda vessel: (
        Path("branch_identity")
        / f"{vessel}_13_substack_boxes_on_labels_with_rings.png"
    ),
    "lumen_diameter": lambda vessel: (
        Path("lumen_diameter") / f"{vessel}_lumen_diameter_distribution.png"
    ),
    "blood_volume_rate": lambda vessel: (
        Path("blood_volume_rate") / f"{vessel}_blood_volume_rate.png"
    ),
}


def build_velocity_bvr_pages(
    report: ReportData, folder_name: str, png_dir: Path | None,
    report_images: ReportImagePaths,
) -> list[Figure]:
    """Build the two mandatory pages of the moments report."""
    return [
        _velocity_page(report, folder_name, png_dir, report_images),
        _vessel_analysis_page(report, folder_name, png_dir),
    ]


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
            fig,
            (left, 0.47, VELOCITY_GRAPH_WIDTH, VELOCITY_GRAPH_HEIGHT),
            png_dir, folder_name,
            "systole", vessel, images, f"{vessel.title()} systole indices",
        )
        velocity_ax = fig.add_axes(
            (left, 0.27, VELOCITY_GRAPH_WIDTH, VELOCITY_GRAPH_HEIGHT)
        )
        plot_filtered_velocity(
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
        _show_report_image(ax, image, fill_bounds=(kind == "systole"))
    if kind in {"vessel_map", "systole"}:
        ax.set_title(title, **GRAPH_TITLE_STYLE)


def plot_filtered_velocity(
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
    prepare_article_fill_layer(ax)
    _plot_beat_envelope(
        ax, phase, samples.T, color="0.55", label="Beat min/max envelope",
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
    ax.set_ylabel(f"Velocity ({_display_unit(unit)})", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.set_box_aspect(1 / ARTICLE_ASPECT_RATIO)


def _finite_parameter(parameter: dict | None) -> float | None:
    if not isinstance(parameter, dict):
        return None
    try:
        value = float(parameter["value"])
    except (KeyError, TypeError, ValueError):
        return None
    return value if np.isfinite(value) else None


def _vessel_analysis_page(
    report: ReportData, folder_name: str, png_dir: Path | None,
) -> Figure:
    fig = plt.figure(figsize=(21 / 2.54, 29.7 / 2.54), dpi=150)
    fig.patch.set_facecolor("white")
    _add_page_header(fig, folder_name, "Vessel Geometry and Blood-Volume Rate")
    for column, vessel in enumerate(("artery", "vein")):
        left = 0.055 + column * 0.48
        fig.text(
            left + 0.205, 0.905, vessel.title(), ha="center", va="center",
            fontsize=13, fontweight="bold",
        )
        _add_artifact_image(
            fig, (left, 0.68, 0.41, 0.18), png_dir, folder_name,
            "histogram", vessel, "Velocity histogram",
        )
        _add_artifact_image(
            fig, (left, 0.43, 0.195, 0.19), png_dir, folder_name,
            "branch_identity", vessel, "Branch identity",
        )
        _add_artifact_image(
            fig, (left + 0.215, 0.43, 0.195, 0.19), png_dir, folder_name,
            "lumen_diameter", vessel, "Lumen diameter distribution",
        )
        _add_artifact_image(
            fig, (left, 0.14, 0.255, 0.21), png_dir, folder_name,
            "blood_volume_rate", vessel, "Mask-derived geometry",
        )
        value_ax = fig.add_axes((left + 0.275, 0.14, 0.135, 0.21))
        value_ax.axis("off")
        _draw_volume_values(value_ax, report.mask_volume_rates.get(vessel))
    return fig


def _add_artifact_image(
    fig: Figure,
    bounds: tuple[float, float, float, float],
    png_dir: Path | None,
    folder_name: str,
    kind: str,
    vessel: str,
    title: str,
) -> None:
    ax = fig.add_axes(bounds)
    ax.axis("off")
    path = _artifact_path(png_dir, folder_name, kind, vessel)
    image = load_first_image(() if path is None else (path,))
    if image is None:
        ax.text(
            0.5, 0.5, f"{title} unavailable", ha="center", va="center",
            transform=ax.transAxes, wrap=True,
        )
    else:
        _show_report_image(ax, image)
    ax.set_title(title, **GRAPH_TITLE_STYLE)


def _draw_volume_values(ax, data: VolumeRateReportData | None) -> None:
    ax.text(0, 0.98, "Values", fontsize=9, fontweight="bold", va="top")
    if data is None:
        ax.text(0, 0.78, "Unavailable", fontsize=8, va="top")
        return
    finite = data.phase_curve[np.isfinite(data.phase_curve)]
    minimum = float(np.min(finite)) if finite.size else None
    maximum = float(np.max(finite)) if finite.size else None
    pulsatile_range = None if not finite.size else maximum - minimum
    lines = (
        f"Median\n{_rate_value(data.overall_value, data.unit)}",
        f"Minimum\n{_rate_value(minimum, data.unit)}",
        f"Maximum\n{_rate_value(maximum, data.unit)}",
        f"Range\n{_rate_value(pulsatile_range, data.unit)}",
    )
    ax.text(
        0, 0.82, "\n\n".join(lines), fontsize=7.5, va="top",
        linespacing=1.25, transform=ax.transAxes,
    )


def _rate_value(value: float | None, unit: str) -> str:
    return "Unavailable" if value is None else f"{value:.2f} {_display_unit(unit)}"


def _plot_beat_envelope(
    ax,
    phase: np.ndarray,
    beat_curves: np.ndarray,
    *,
    color: str,
    label: str | None = None,
) -> None:
    """Fill the phase-wise range of finite beat samples."""
    values = np.asarray(beat_curves, dtype=float)
    if values.ndim != 2 or values.shape[0] != phase.size:
        return
    masked = np.ma.masked_invalid(values)
    minimum = np.asarray(masked.min(axis=1).filled(np.nan), dtype=float)
    maximum = np.asarray(masked.max(axis=1).filled(np.nan), dtype=float)
    valid = np.isfinite(minimum) & np.isfinite(maximum)
    if not np.any(valid):
        return
    ax.fill_between(
        phase, minimum, maximum, where=valid,
        color=color, alpha=0.18, linewidth=0,
        label=label, zorder=0,
    )


def _add_page_header(fig: Figure, folder_name: str, section: str) -> None:
    """Identify both the acquisition and the subject of a report page."""
    section = section[:1].upper() + section[1:]
    fig.text(0.5, 0.978, folder_name, fontsize=10, color="0.35",
             ha="center", va="top")
    fig.text(0.5, 0.952, section, fontsize=16, fontweight="bold",
             ha="center", va="top")


def _show_report_image(
    ax, image: np.ndarray, *, fill_bounds: bool = False,
) -> None:
    """Draw an image, optionally filling a matching-aspect report box."""
    height, width = np.asarray(image).shape[:2]
    ax.imshow(image, aspect="auto" if fill_bounds else "equal")
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
    return f"{value:.{decimals}f}{f' {_display_unit(unit)}' if unit else ''}"


def _display_unit(unit: str) -> str:
    """Format report units without changing their stored representation."""
    return display_unit(unit)


def _artifact_path(
    png_dir: Path | None,
    folder_name: str,
    kind: str,
    vessel: str,
) -> Path | None:
    if png_dir is None:
        return None
    suffix = SECOND_PAGE_ARTIFACTS[kind](vessel)
    return Path(png_dir) / prefixed_artifact_path(suffix, folder_name)


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
