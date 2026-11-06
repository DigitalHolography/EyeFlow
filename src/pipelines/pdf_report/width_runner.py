"""Once-per-run PDF for persisted motion and lumen-width results."""

from __future__ import annotations

from input_output.output_manager import OutputType
from input_output.reports import write_pdf_report
from input_output.writers.artifact_names import prefixed_artifact_path

from .width_data import extract_width_report_data
from .width_pages import build_width_page


def run_width_pdf_report(ctx) -> None:
    """Write the third report only when width or motion metrics exist."""
    if ctx.output.manager is None:
        ctx.log_warning("No output manager available for width PDF generation")
        return
    report = extract_width_report_data(ctx.output.h5)
    if report is None:
        ctx.log("Width PDF skipped: no displacement-profile or spatial-gradient metrics.")
        return
    stem = ctx.output.manager.layout.stem
    path = prefixed_artifact_path(
        ctx.output.manager.dir_for(OutputType.PDF) / "width_report.pdf", stem,
    )
    write_pdf_report(path, [build_width_page(report, stem)])
    ctx.log(f"Width PDF report generated: {path}")
