"""Standalone figure artifacts owned by the moments PDF pipeline."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt

from input_output.writers.article import save_article_figure_pair
from plotting.article import ARTICLE_ASPECT_RATIO, ARTICLE_DPI

from .data import ReportData
from .pages import plot_filtered_velocity


def export_moments_report_figures(output, report: ReportData) -> list[Path]:
    """Export report-native graphs as matching PNG and EPS artifacts."""

    paths: list[Path] = []
    for vessel in ("artery", "vein"):
        fig, ax = plt.subplots(figsize=(7.5, 7.5 / ARTICLE_ASPECT_RATIO))
        fig.patch.set_facecolor("white")
        plot_filtered_velocity(
            ax,
            report.velocity_per_beat.get(vessel),
            report.velocity_units.get(vessel, "mm/s"),
            f"{vessel.title()} filtered per-beat velocity",
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
