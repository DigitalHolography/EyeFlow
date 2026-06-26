"""PDF report generation runner."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from input_output.output_manager import OutputType
from input_output.reports import generate_a4_report
from pipeline_engine.imports import ProcessResult


def run_pdf_report(ctx) -> ProcessResult:
    """Generate an A4 PDF report from the output H5 file."""
    ctx.require_inputs("hd", "dv")
    
    output_h5_path = ctx.output.h5.filename
    if output_h5_path is None:
        ctx.log("No output H5 file available for PDF generation")
        return ProcessResult(metrics={})

    output_manager = ctx.output.manager
    if output_manager is None:
        ctx.log("No output manager available for PDF generation")
        return ProcessResult(metrics={})

    layout = output_manager.layout
    folder_name = layout.stem
    ef_dir = layout.ef_dir
    report_dir = output_manager.dir_for(OutputType.PDF)
    png_dir = ef_dir / "png"

    hd_png_dir = _hd_png_dir(ctx.inputs.hd.filename)
    mask_dir = hd_png_dir / "mask"
    
    try:
        pdf_path = generate_a4_report(
            output_h5_paths=[Path(output_h5_path)],
            output_dir=report_dir,
            png_dir=png_dir if png_dir.exists() else None,
            hd_png_dir=hd_png_dir if hd_png_dir.exists() else None,
            mask_dir=mask_dir if mask_dir.exists() else None,
            folder_name=folder_name,
        )
        ctx.log(f"PDF report generated: {pdf_path}")
        
        return ProcessResult(
            metrics={
                "pdf_report_path": str(pdf_path),
                "pdf_report_generated": np.uint8(1),
            }
        )
    except Exception as e:
        ctx.log(f"PDF generation failed: {e}")
        return ProcessResult(
            metrics={
                "pdf_report_path": "",
                "pdf_report_generated": np.uint8(0),
            }
        )


def _hd_png_dir(hd_h5_filename: str | None) -> Path:
    if hd_h5_filename is None:
        return Path()
    return Path(hd_h5_filename).parent.parent / "png"
