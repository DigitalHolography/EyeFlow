"""Semantic, self-contained HTML presentation of the moments report data."""

from __future__ import annotations

from html import escape
from pathlib import Path

import matplotlib.pyplot as plt

from input_output.reports.html_report import (
    figure_png_data_uri,
    image_png_data_uri,
)
from input_output.writers.png import ReportImagePaths
from plotting.article import ARTICLE_ASPECT_RATIO

from .assets import (
    ANALYSIS_GRAPH_TITLES,
    load_analysis_artifact,
    load_velocity_artifact,
)
from .data import ReportData
from .plots import plot_filtered_velocity, plot_systole_indices
from .summaries import (
    analysis_summary_rows,
    format_plain_value,
    heart_rate_summary,
    velocity_summary_rows,
)


def build_velocity_bvr_html(
    report: ReportData,
    folder_name: str,
    png_dir: Path | None,
    report_images: ReportImagePaths,
) -> str:
    """Build an expandable HTML report independently from the PDF pages."""
    velocity_columns = "".join(
        _velocity_column(report, folder_name, png_dir, report_images, vessel)
        for vessel in ("artery", "vein")
    )
    analysis_columns = "".join(
        _analysis_column(report, folder_name, png_dir, vessel)
        for vessel in ("artery", "vein")
    )
    heart_rate = escape(format_plain_value(heart_rate_summary(report.parameters)))
    velocity_summary = "".join(
        _summary_block(
            vessel.title(), velocity_summary_rows(report.parameters, vessel)
        )
        for vessel in ("artery", "vein")
    )
    title = f"{folder_name} — Velocity and Blood-Volume Rate"
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{escape(title)}</title>
  <style>{_STYLES}</style>
</head>
<body>
  <header class="report-header">
    <p class="acquisition">{escape(folder_name)}</p>
    <h1>Velocity and Blood-Volume Rate</h1>
  </header>
  <main>
    <section id="velocity">
      <h2>Velocity</h2>
      <div class="vessel-grid">{velocity_columns}</div>
      <section class="summary-card wide-summary">
        <div class="summary-heading">
          <h3>Summary</h3>
          <p><strong>Heart rate</strong> <span>{heart_rate}</span></p>
        </div>
        <div class="summary-grid">{velocity_summary}</div>
      </section>
    </section>
    <section id="vessel-analysis">
      <h2>Vessel Geometry and Blood-Volume Rate</h2>
      <div class="vessel-grid">{analysis_columns}</div>
    </section>
  </main>
</body>
</html>
"""


def _velocity_column(
    report: ReportData,
    folder_name: str,
    png_dir: Path | None,
    report_images: ReportImagePaths,
    vessel: str,
) -> str:
    vessel_title = vessel.title()
    vessel_map = load_velocity_artifact(
        png_dir, folder_name, "vessel_map", vessel, report_images
    )
    systole_uri = _systole_plot_uri(report, vessel)
    velocity_uri = _velocity_plot_uri(report, vessel)
    return f"""
<article class="vessel-column">
  <h3>{vessel_title}</h3>
  {_artifact_figure(f"{vessel_title} RMS-Frequency Map", vessel_map)}
  {_uri_figure(f"{vessel_title} Systole Indices", systole_uri)}
  {_uri_figure(f"{vessel_title} Per-Beat Velocity", velocity_uri)}
</article>
"""


def _analysis_column(
    report: ReportData,
    folder_name: str,
    png_dir: Path | None,
    vessel: str,
) -> str:
    vessel_title = vessel.title()
    figures = tuple(
        _artifact_figure(
            f"{vessel_title} {title}",
            load_analysis_artifact(png_dir, folder_name, kind, vessel),
        )
        for kind, title in ANALYSIS_GRAPH_TITLES.items()
    )
    paired_figures = "".join(figures[:2])
    remaining_figures = "".join(figures[2:])
    summary = _summary_block(
        f"{vessel_title} Summary", analysis_summary_rows(report, vessel)
    )
    return f"""
<article class="vessel-column">
  <h3>{vessel_title}</h3>
  <div class="paired-graphs">{paired_figures}</div>
  {remaining_figures}
  <section class="summary-card">{summary}</section>
</article>
"""


def _systole_plot_uri(report: ReportData, vessel: str) -> str:
    fig, ax = _report_figure()
    plot_systole_indices(
        ax,
        report.systole_signals.get(vessel),
        f"{vessel.title()} Systole Indices",
    )
    return figure_png_data_uri(fig)


def _velocity_plot_uri(report: ReportData, vessel: str) -> str:
    fig, ax = _report_figure()
    plot_filtered_velocity(
        ax,
        report.velocity_per_beat.get(vessel),
        report.velocity_units.get(vessel, "mm/s"),
        f"{vessel.title()} Per-Beat Velocity",
        unfiltered_values=report.velocity_per_beat_unfiltered.get(vessel),
    )
    return figure_png_data_uri(fig)


def _report_figure():
    fig, ax = plt.subplots(figsize=(7.5, 7.5 / ARTICLE_ASPECT_RATIO))
    fig.patch.set_facecolor("white")
    return fig, ax


def _artifact_figure(title: str, image) -> str:
    if image is None:
        return _unavailable_figure(title)
    return _uri_figure(title, image_png_data_uri(image))


def _uri_figure(title: str, uri: str) -> str:
    safe_title = escape(title)
    return f"""
<figure class="graph-card">
  <figcaption>{safe_title}</figcaption>
  <img src="{uri}" alt="{safe_title}">
</figure>
"""


def _unavailable_figure(title: str) -> str:
    safe_title = escape(title)
    return f"""
<figure class="graph-card unavailable">
  <figcaption>{safe_title}</figcaption>
  <div role="img" aria-label="{safe_title} unavailable">Unavailable</div>
</figure>
"""


def _summary_block(title: str, rows) -> str:
    body = "".join(
        "<tr>"
        f"<th scope=\"row\">{escape(row.label)}</th>"
        f"<td>{escape(format_plain_value(row.value))}</td>"
        "</tr>"
        for row in rows
    )
    return f"""
<div class="summary-block">
  <h4>{escape(title)}</h4>
  <table><tbody>{body}</tbody></table>
</div>
"""


_STYLES = """
:root {
  color-scheme: light;
  font-family: Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont,
    "Segoe UI", sans-serif;
  color: #171717;
  background: white;
}
* { box-sizing: border-box; }
body { margin: 0; }
.report-header, main {
  width: min(1500px, calc(100% - 32px));
  margin-inline: auto;
}
.report-header {
  padding: 36px 8px 20px;
  text-align: center;
}
.acquisition { margin: 0 0 6px; color: #666; font-size: 0.9rem; }
h1, h2, h3, h4 { margin-top: 0; }
h1 { margin-bottom: 0; font-size: clamp(1.6rem, 3vw, 2.5rem); }
h2 { font-size: 1.5rem; margin-bottom: 24px; }
h3 { font-size: 1.2rem; }
h4 { font-size: 1rem; margin-bottom: 12px; }
main { padding-bottom: 48px; }
main > section { padding-block: 24px 40px; }
main > section + section {
  padding-top: 48px;
  border-top: 1px solid #d8d8d8;
}
.vessel-grid, .summary-grid {
  display: grid;
  grid-template-columns: repeat(2, minmax(0, 1fr));
  gap: 28px;
}
.vessel-column > h3 {
  padding-bottom: 10px;
  border-bottom: 2px solid #222;
}
.graph-card {
  margin: 0 0 22px;
  padding: 14px;
  border: 1px solid #dedede;
  border-radius: 8px;
  background: #fff;
}
.paired-graphs {
  display: grid;
  grid-template-columns: minmax(0, 3fr) minmax(0, 2fr);
  gap: 14px;
}
.paired-graphs .graph-card { min-width: 0; }
.paired-graphs .graph-card figcaption { font-size: 0.85rem; }
.graph-card figcaption {
  margin-bottom: 12px;
  font-weight: 600;
}
.graph-card img {
  display: block;
  width: 100%;
  min-height: 160px;
  max-height: 480px;
  object-fit: contain;
}
.graph-card.unavailable div {
  display: grid;
  min-height: 180px;
  place-items: center;
  color: #777;
  background: #f4f4f4;
}
.summary-card {
  padding: 18px;
  border: 1px solid #d4d4d4;
  border-radius: 8px;
  background: #fafafa;
}
.wide-summary { margin-top: 8px; }
.summary-heading {
  display: flex;
  align-items: baseline;
  justify-content: space-between;
  gap: 20px;
}
.summary-heading p { margin: 0; }
.summary-block table { width: 100%; border-collapse: collapse; }
.summary-block th, .summary-block td {
  padding: 7px 0;
  border-bottom: 1px solid #e2e2e2;
}
.summary-block th { color: #555; font-weight: 400; text-align: left; }
.summary-block td { font-weight: 700; text-align: right; }
@media (max-width: 850px) {
  .vessel-grid, .summary-grid { grid-template-columns: 1fr; }
}
@media (max-width: 520px) {
  .paired-graphs { grid-template-columns: 1fr; }
}
"""


__all__ = ["build_velocity_bvr_html"]
