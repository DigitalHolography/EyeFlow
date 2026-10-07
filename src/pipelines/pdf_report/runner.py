"""PDF report generation runner."""

from __future__ import annotations

from pathlib import Path

from input_output.inputs import sidecar_dir_for_h5
from input_output.output_manager import OutputType
from input_output.reports import generate_a4_report


def run_pdf_report(ctx) -> None:
    """Generate an A4 PDF report from the output H5 file."""
    ctx.require_inputs("hd", "dv")

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

    hd_png_dir = (
        sidecar_dir_for_h5(ctx.inputs.hd.filename, "png")
        if ctx.inputs.hd.filename is not None
        else None
    )
    mask_dir = hd_png_dir / "mask" if hd_png_dir is not None else None

    try:
        report_options = {}
        processing_root = getattr(ctx.output.h5, "processing_root", None)
        if processing_root is not None:
            report_options["processing_root"] = processing_root
        pdf_path = generate_a4_report(
            output_h5_path=Path(output_h5_path),
            output_dir=report_dir,
            folder_name=folder_name,
            png_dir=png_dir if png_dir.exists() else None,
            hd_png_dir=hd_png_dir if hd_png_dir and hd_png_dir.exists() else None,
            mask_dir=mask_dir if mask_dir and mask_dir.exists() else None,
            **report_options,
        )
        ctx.log(f"PDF report generated: {pdf_path}")
        return None
    except Exception as e:
        ctx.log_error(f"PDF generation failed: {e}")
        raise
