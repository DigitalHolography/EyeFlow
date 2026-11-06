"""PDF report generation runner."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from input_output.output_manager import OutputType
from input_output.reports import write_pdf_report
from input_output.writers.artifact_names import prefixed_artifact_path
from input_output.writers.png import REPORT_IMAGES_STATE

from .data import ReportData, extract_report_data
from .pages import build_report_pages


def run_pdf_report(ctx) -> None:
    """Assemble and save the selected workflow's PDF report."""
    output_h5_path = ctx.output.h5.filename
    if output_h5_path is None:
        ctx.log_warning("No output H5 file available for PDF generation")
        return None

    output_manager = ctx.output.manager
    if output_manager is None:
        ctx.log_warning("No output manager available for PDF generation")
        return None

    layout = output_manager.layout
    folder_name = layout.stem
    report_dir = output_manager.dir_for(OutputType.PDF)
    png_dir = output_manager.dir_for(OutputType.PNG)
    report_images = ctx.state.get(REPORT_IMAGES_STATE, {})

    try:
        report = extract_report_data(ctx.output.h5)
        _warn_missing_velocities(ctx, report)
        pages = build_report_pages(report, folder_name, png_dir, report_images)
        pdf_path = write_pdf_report(
            prefixed_artifact_path(Path(report_dir) / "report.pdf", folder_name), pages,
        )
        ctx.log(f"PDF report generated: {pdf_path}")
        return None
    except Exception as e:
        ctx.log_error(f"PDF generation failed: {e}")
        raise


def _warn_missing_velocities(ctx, report: ReportData) -> None:
    """Explain missing report traces without substituting unfiltered velocity."""
    for vessel in ("artery", "vein"):
        values = report.velocity_per_beat.get(vessel)
        if values is None:
            ctx.log_warning(
                f"PDF: {vessel} filtered per-beat velocity is unavailable; "
                "no unfiltered-velocity fallback was used."
            )
        elif not np.isfinite(values).any():
            ctx.log_warning(
                f"PDF: {vessel} filtered per-beat velocity has no finite samples."
            )
