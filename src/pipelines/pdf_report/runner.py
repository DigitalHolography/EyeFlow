"""Moments PDF and HTML report generation runner."""

from __future__ import annotations

import numpy as np

from input_output.output_manager import OutputType
from input_output.reports import write_html_report, write_pdf_report
from input_output.writers.artifact_names import prefixed_artifact_path
from input_output.writers.png import REPORT_IMAGES_STATE
from pipeline_engine import ProductSpec, StageGraph, StageSpec
from pipeline_engine.context import require_managed_runtime

from .artifacts import export_moments_report_figures
from .data import ReportData, extract_report_data
from .html import build_velocity_bvr_html
from .pages import build_velocity_bvr_pages

MOMENTS_VARIANT = "doppler_moments"


def _report_data(ctx) -> None:
    variant = getattr(getattr(ctx, "execution_variant", None), "name", None)
    if variant not in (None, MOMENTS_VARIANT):
        ctx.dag.skip(f"Moments reports are not produced for velocity workflow '{variant}'.")

    if ctx.output.h5.filename is None:
        ctx.dag.skip("No output H5 file is available for report generation.")
    if ctx.output.manager is None:
        ctx.dag.skip("No output manager is available for report generation.")

    report = extract_report_data(ctx.output.h5)
    _warn_missing_velocities(ctx, report)
    ctx.dag.handoff("report_data", report)


def _report_figures(ctx) -> None:
    export_moments_report_figures(ctx.output, ctx.dag.require("report_data"))


def _pdf_report(ctx) -> None:
    report = ctx.dag.require("report_data")
    manager = ctx.output.manager
    folder_name = manager.layout.stem
    pages = build_velocity_bvr_pages(
        report,
        folder_name,
        manager.dir_for(OutputType.PNG),
        ctx.state.get(REPORT_IMAGES_STATE, {}),
    )
    path = write_pdf_report(
        prefixed_artifact_path(
            manager.layout.ef_dir / OutputType.PDF.value / "velocityBvrReport.pdf",
            folder_name,
        ),
        pages,
    )
    ctx.log(f"Velocity/BVR PDF report generated: {path}")


def _html_report(ctx) -> None:
    report = ctx.dag.require("report_data")
    manager = ctx.output.manager
    folder_name = manager.layout.stem
    document = build_velocity_bvr_html(
        report,
        folder_name,
        manager.dir_for(OutputType.PNG),
        ctx.state.get(REPORT_IMAGES_STATE, {}),
    )
    path = write_html_report(
        prefixed_artifact_path(
            manager.layout.ef_dir / OutputType.HTML.value / "velocityBvrReport.html",
            folder_name,
        ),
        document,
    )
    ctx.log(f"Velocity/BVR HTML report generated: {path}")


def pdf_report_graph(*, options=()):
    """Build independent PDF, HTML, and reusable-figure report stages."""

    return StageGraph(
        (
            StageSpec(
                "data",
                _report_data,
                products=(ProductSpec("report_data"),),
            ),
            StageSpec(
                "figures",
                _report_figures,
                requires=("report_data",),
                artifact=True,
            ),
            StageSpec(
                "pdf",
                _pdf_report,
                requires=("report_data",),
                artifact=True,
            ),
            StageSpec(
                "html",
                _html_report,
                requires=("report_data",),
                artifact=True,
            ),
        )
    )


def run_pdf_report(ctx):
    require_managed_runtime(ctx, "pdf_report")
    return pdf_report_graph(options=ctx.options_for("pdf_report"))


def _warn_missing_velocities(ctx, report: ReportData) -> None:
    """Explain missing traces without substituting unfiltered velocity."""

    for vessel in ("artery", "vein"):
        values = report.velocity_per_beat.get(vessel)
        if values is None:
            ctx.log_warning(
                f"Report: {vessel} filtered per-beat velocity is unavailable; "
                "no unfiltered-velocity fallback was used."
            )
        elif not np.isfinite(values).any():
            ctx.log_warning(
                f"Report: {vessel} filtered per-beat velocity has no finite samples."
            )
