"""Standalone figure artifacts owned by the moments report pipeline."""

from __future__ import annotations

from pathlib import Path

from input_output.writers.article import save_article_figure_pair
from plotting.article import ARTICLE_ASPECT_RATIO, ARTICLE_DPI
from plotting.headless import create_agg_figure

from .data import ReportData
from .plots import plot_filtered_velocity


def export_moments_report_figures(output, report: ReportData) -> list[Path]:
    """Export report-native graphs as matching PNG and EPS artifacts."""

    paths: list[Path] = []
    for vessel in ("artery", "vein"):
        fig = create_agg_figure(figsize=(7.5, 7.5 / ARTICLE_ASPECT_RATIO))
        ax = fig.subplots()
        fig.patch.set_facecolor("white")
        plot_filtered_velocity(
            ax,
            report.velocity_per_beat.get(vessel),
            report.velocity_units.get(vessel, "mm/s"),
            f"{vessel.title()} Per-Beat Velocity",
            unfiltered_values=report.velocity_per_beat_unfiltered.get(vessel),
        )
        paths.extend(
            save_article_figure_pair(
                output,
                fig,
                f"pdf_report/{vessel}_filtered_per_beat_velocity",
                dpi=ARTICLE_DPI,
            )
        )
    return paths


__all__ = ["export_moments_report_figures"]
