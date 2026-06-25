"""PDF report generation runner."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from input_output.reports import generate_a4_report
from pipeline_engine.imports import ProcessResult


def run_pdf_report(ctx) -> ProcessResult:
    """Generate an A4 PDF report from the output H5 file."""
    ctx.require_inputs("hd", "dv")
    
    output_h5_path = ctx.output.h5.filename
    if output_h5_path is None:
        ctx.log("No output H5 file available for PDF generation")
        return ProcessResult(metrics={})
    
    holo_path = ctx.inputs.hd.filename
    if holo_path is None:
        ctx.log("No HoloDoppler input available for PNG paths")
        return ProcessResult(metrics={})
    
    holo_path_obj = Path(holo_path)
    folder_name = holo_path_obj.stem
    output_dir = Path(output_h5_path).parent.parent
    
    ef_dir = output_dir / f"{folder_name}_EF"
    png_dir = ef_dir / "png"
    hd_png_dir = Path(holo_path_obj.parent) / "HD" / "png"
    mask_dir = hd_png_dir / "mask"
    
    try:
        pdf_path = generate_a4_report(
            output_h5_paths=[Path(output_h5_path)],
            output_dir=ef_dir,
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