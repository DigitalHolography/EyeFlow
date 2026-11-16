"""Moments report generation runner."""

from __future__ import annotations

import numpy as np

from input_output.output_manager import OutputType
from input_output.reports import write_html_report, write_pdf_report
from input_output.writers.artifact_names import prefixed_artifact_path
from input_output.writers.png import REPORT_IMAGES_STATE

from .artifacts import export_moments_report_figures
from .data import ReportData, extract_report_data
from .html import build_velocity_bvr_html
from .pages import build_velocity_bvr_pages

MOMENTS_VARIANT = "doppler_moments"


def run_pdf_report(ctx) -> None:
    """Write moments velocity/BVR reports as independent PDF and HTML outputs."""
    variant = getattr(getattr(ctx, "execution_variant", None), "name", None)
    if variant not in (None, MOMENTS_VARIANT):
        ctx.log(
            f"Moments reports skipped for velocity workflow '{variant}'."
        )
        return None

    output_h5_path = ctx.output.h5.filename
    if output_h5_path is None:
        ctx.log_warning("No output H5 file available for report generation")
        return None

    output_manager = ctx.output.manager
    if output_manager is None:
        ctx.log_warning("No output manager available for report generation")
        return None

    layout = output_manager.layout
    folder_name = layout.stem
    report_dir = layout.ef_dir / OutputType.PDF.value
    html_dir = layout.ef_dir / OutputType.HTML.value

    try:
        report = extract_report_data(ctx.output.h5)
        _warn_missing_velocities(ctx, report)
        export_moments_report_figures(ctx.output, report)
        png_dir = output_manager.dir_for(OutputType.PNG)
        report_images = ctx.state.get(REPORT_IMAGES_STATE, {})
        pages = build_velocity_bvr_pages(
            report, folder_name, png_dir, report_images,
        )
        path = write_pdf_report(
            prefixed_artifact_path(
                report_dir / "velocityBvrReport.pdf", folder_name,
            ),
            pages,
        )
        ctx.log(f"Velocity/BVR PDF report generated: {path}")
        html_path = write_html_report(
            prefixed_artifact_path(
                html_dir / "velocityBvrReport.html", folder_name,
            ),
            build_velocity_bvr_html(
                report, folder_name, png_dir, report_images,
            ),
        )
        ctx.log(f"Velocity/BVR HTML report generated: {html_path}")
        return None
    except Exception as e:
        ctx.log_error(f"Report generation failed: {e}")
        raise


def _warn_missing_velocities(ctx, report: ReportData) -> None:
    """Explain missing report traces without substituting unfiltered velocity."""
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
