"""A4 page layout for the two-page moments report."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure

from input_output.writers.png import ReportImagePaths
from plotting.article import ARTICLE_ASPECT_RATIO
from plotting.units import display_unit

from .assets import (
    ANALYSIS_GRAPH_TITLES,
    load_analysis_artifact,
    load_velocity_artifact,
)
from .data import (
    ReportData,
)
from .plots import GRAPH_TITLE_STYLE, plot_filtered_velocity, plot_systole_indices
from .summaries import (
    SummaryValue,
    analysis_summary_rows,
    heart_rate_summary,
    velocity_summary_rows,
)
SUMMARY_TITLE_SIZE = 12
SUMMARY_HEADING_SIZE = 10
SUMMARY_LABEL_SIZE = 8.2
SUMMARY_VALUE_SIZE = 8.5
SUMMARY_LABEL_COLOR = "0.25"
VELOCITY_GRAPH_WIDTH = 0.39
VELOCITY_GRAPH_HEIGHT = (
    VELOCITY_GRAPH_WIDTH * 21 / (29.7 * ARTICLE_ASPECT_RATIO)
)
SECOND_PAGE_GRAPH_WIDTH = 0.41
SECOND_PAGE_GOLDEN_HEIGHT = (
    SECOND_PAGE_GRAPH_WIDTH * 21 / (29.7 * ARTICLE_ASPECT_RATIO)
)

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
            f"{vessel.title()} RMS-Frequency Map",
        )
        systole_ax = fig.add_axes(
            (left, 0.485, VELOCITY_GRAPH_WIDTH, VELOCITY_GRAPH_HEIGHT)
        )
        plot_systole_indices(
            systole_ax,
            report.systole_signals.get(vessel),
            f"{vessel.title()} Systole Indices",
        )
        velocity_ax = fig.add_axes(
            (left, 0.26, VELOCITY_GRAPH_WIDTH, VELOCITY_GRAPH_HEIGHT)
        )
        plot_filtered_velocity(
            velocity_ax,
            report.velocity_per_beat.get(vessel),
            report.velocity_units.get(vessel, "mm/s"),
            f"{vessel.title()} Per-Beat Velocity",
            unfiltered_values=report.velocity_per_beat_unfiltered.get(vessel),
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
    image = load_velocity_artifact(
        png_dir, folder_name, kind, vessel, report_images
    )
    if image is None:
        ax.text(0.5, 0.5, f"{title} unavailable", ha="center", va="center",
                transform=ax.transAxes, wrap=True)
    else:
        _show_report_image(ax, image)
    if kind == "vessel_map":
        ax.set_title(title, **GRAPH_TITLE_STYLE)


def _vessel_analysis_page(
    report: ReportData, folder_name: str, png_dir: Path | None,
) -> Figure:
    fig = plt.figure(figsize=(21 / 2.54, 29.7 / 2.54), dpi=150)
    fig.patch.set_facecolor("white")
    _add_page_header(fig, folder_name, "Vessel Geometry and Blood-Volume Rate")
    for column, vessel in enumerate(("artery", "vein")):
        left = 0.055 + column * 0.48
        graph_specs = (
            (0.750, 0.130, "histogram"),
            (0.590, 0.130, "branch_identity"),
            (0.380, SECOND_PAGE_GOLDEN_HEIGHT, "lumen_diameter"),
            (0.170, SECOND_PAGE_GOLDEN_HEIGHT, "blood_volume_rate"),
        )
        for bottom, height, kind in graph_specs:
            _add_artifact_image(
                fig, (left, bottom, SECOND_PAGE_GRAPH_WIDTH, height),
                png_dir, folder_name, kind, vessel,
                f"{vessel.title()} {ANALYSIS_GRAPH_TITLES[kind]}",
            )
        value_ax = fig.add_axes((left, 0.003, SECOND_PAGE_GRAPH_WIDTH, 0.150))
        value_ax.axis("off")
        _draw_second_page_summary(
            value_ax, report, vessel, show_title=(vessel == "artery")
        )
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
    image = load_analysis_artifact(png_dir, folder_name, kind, vessel)
    if image is None:
        ax.text(
            0.5, 0.5, f"{title} unavailable", ha="center", va="center",
            transform=ax.transAxes, wrap=True,
        )
    else:
        _show_report_image(ax, image, zoom=1.0)
    ax.set_title(title, **GRAPH_TITLE_STYLE)


def _draw_second_page_summary(
    ax, report: ReportData, vessel: str, *, show_title: bool,
) -> None:
    """Summarize every page-two product with the page-one value styling."""
    if show_title:
        ax.text(
            0.02, 0.98, "Summary", fontsize=SUMMARY_TITLE_SIZE,
            fontweight="bold", va="top",
        )
    for row, item in enumerate(analysis_summary_rows(report, vessel)):
        y = 0.80 - row * 0.105
        ax.text(
            0.02, y, item.label, fontsize=SUMMARY_LABEL_SIZE,
            color=SUMMARY_LABEL_COLOR, va="center",
        )
        ax.text(
            0.98, y, _pdf_value(item.value), fontsize=SUMMARY_VALUE_SIZE,
            fontweight="normal" if item.value.regular_unit else "bold",
            ha="right", va="center",
        )


def _add_page_header(fig: Figure, folder_name: str, section: str) -> None:
    """Identify both the acquisition and the subject of a report page."""
    section = section[:1].upper() + section[1:]
    fig.text(0.5, 0.978, folder_name, fontsize=10, color="0.35",
             ha="center", va="top")
    fig.text(0.5, 0.952, section, fontsize=16, fontweight="bold",
             ha="center", va="top")


def _show_report_image(
    ax, image: np.ndarray, *, zoom: float = 1.04,
) -> None:
    """Draw an image without distorting its aspect ratio."""
    height, width = np.asarray(image).shape[:2]
    ax.imshow(image, aspect="equal")
    ax.set_anchor("C")
    ax.set_adjustable("box")
    ax.margins(0)

    if width > 0 and height > 0:
        x_center = (width - 1) / 2
        y_center = (height - 1) / 2
        half_width = width / (2 * zoom)
        half_height = height / (2 * zoom)
        ax.set_xlim(x_center - half_width, x_center + half_width)
        ax.set_ylim(y_center + half_height, y_center - half_height)


def _add_parameters_section(ax, parameters: dict) -> None:
    """Show paired vessel measurements and one shared cardiac value."""
    ax.text(
        0, 0.98, "Summary", fontsize=SUMMARY_TITLE_SIZE,
        fontweight="bold", va="top",
    )
    ax.text(
        1, 0.98,
        f"Heart rate  |  {_pdf_value(heart_rate_summary(parameters))}",
        fontsize=SUMMARY_HEADING_SIZE, fontweight="bold", ha="right", va="top",
    )
    _draw_vessel_column(ax, parameters, "artery", 0.00)
    _draw_vessel_column(ax, parameters, "vein", 0.52)


def _draw_vessel_column(
    ax, parameters: dict, vessel: str, left: float,
) -> None:
    arterial = vessel == "artery"
    ax.text(
        left + 0.025, 0.68, "Artery" if arterial else "Vein",
        color="black", fontsize=SUMMARY_HEADING_SIZE,
        fontweight="bold", va="center",
    )
    for row, item in enumerate(velocity_summary_rows(parameters, vessel)):
        y = 0.56 - row * 0.115
        ax.text(
            left + 0.025, y, item.label, fontsize=SUMMARY_LABEL_SIZE,
            color=SUMMARY_LABEL_COLOR, va="center",
        )
        ax.text(
            left + 0.455, y, _pdf_value(item.value),
            fontsize=SUMMARY_VALUE_SIZE, fontweight="bold",
            ha="right", va="center",
        )


def _pdf_value(value: SummaryValue) -> str:
    if value.value is None:
        return "Unavailable"
    number = f"{value.value:.{value.decimals}f}"
    if value.regular_unit and value.unit == "µm":
        return rf"$\mathbf{{{number}}}\ \mathrm{{\mu m}}$"
    return f"{number}{f' {_display_unit(value.unit)}' if value.unit else ''}"


def _display_unit(unit: str) -> str:
    """Format report units without changing their stored representation."""
    return display_unit(unit)


__all__ = ["build_velocity_bvr_pages"]
