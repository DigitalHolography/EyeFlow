"""PDF report generation runner."""

from __future__ import annotations

import numpy as np

from input_output.output_manager import OutputType
from input_output.reports import write_pdf_report
from input_output.writers.artifact_names import prefixed_artifact_path
from input_output.writers.png import REPORT_IMAGES_STATE

from .data import ReportData, extract_report_data
from .pages import build_velocity_bvr_pages
from .spectra import build_spectra_gaps_pages, extract_spectra_gaps_report_data

MOMENTS_VARIANT = "doppler_moments"


def run_pdf_report(ctx) -> None:
    """Write the moments velocity/BVR report and one shared spectra report."""
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
    report_dir = layout.ef_dir / OutputType.PDF.value
    variant = getattr(getattr(ctx, "execution_variant", None), "name", None)

    try:
        if variant in (None, MOMENTS_VARIANT):
            report = extract_report_data(ctx.output.h5)
            _warn_missing_velocities(ctx, report)
            pages = build_velocity_bvr_pages(
                report, folder_name, output_manager.dir_for(OutputType.PNG),
                ctx.state.get(REPORT_IMAGES_STATE, {}),
            )
            path = write_pdf_report(
                prefixed_artifact_path(
                    report_dir / "velocityBvrReport.pdf", folder_name,
                ),
                pages,
            )
            ctx.log(f"Velocity/BVR PDF report generated: {path}")

        spectra_path = prefixed_artifact_path(
            report_dir / "spectraGapsReport.pdf", folder_name,
        )
        spectra = extract_spectra_gaps_report_data(ctx.output.h5)
        if spectra is not None and not spectra_path.exists():
            write_pdf_report(
                spectra_path, build_spectra_gaps_pages(spectra, folder_name),
            )
            ctx.log(f"Spectra/gaps PDF report generated: {spectra_path}")
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
