"""A4 PDF report generator mimicking MATLAB generateA4Report."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from PIL import Image

from input_output.writers.artifact_names import prefixed_artifact_path
from input_output.writers.png import ReportImagePaths

from .pdf_data import (
    _extract_parameters_from_h5,
    _format_parameters_for_display,
    _velocity_dataset_unit,
)


def generate_a4_report(
    output_h5_path: Path,
    output_dir: Path,
    folder_name: str,
    png_dir: Path | None = None,
    hd_png_dir: Path | None = None,
    mask_dir: Path | None = None,
    processing_root: str | None = None,
    report_images: ReportImagePaths | None = None,
) -> Path:
    """
    Generate an A4 PDF report mimicking the MATLAB generateA4Report function.
    
    This creates a comprehensive report with:
    - Title with folder name
    - Image grid (vessel map, RI plot, systole indices for artery and vein)
    - Parameters table with computed metrics

    Args:
        output_h5_path: Output H5 file from processing
        output_dir: Directory to save the PDF
        folder_name: Name of the folder (for title)
        png_dir: Directory containing PNG outputs
        hd_png_dir: Directory containing HoloDoppler M0 PNGs
        mask_dir: Directory containing mask PNGs
        processing_root: Optional workflow processing group to read
        report_images: Actual PNG paths registered by this workflow's writers

    Returns:
        Path to the generated PDF file
    """
    output_h5_path = Path(output_h5_path)
    output_dir = Path(output_dir)
    png_dir = Path(png_dir) if png_dir is not None else None
    hd_png_dir = Path(hd_png_dir) if hd_png_dir is not None else None
    mask_dir = Path(mask_dir) if mask_dir is not None else None
    images = report_images or {}

    # Create output path
    pdf_path = prefixed_artifact_path(output_dir / "report.pdf", folder_name)
    pdf_path.parent.mkdir(parents=True, exist_ok=True)

    # Extract parameters from the H5 file
    parameters = _extract_parameters_from_h5(output_h5_path, processing_root=processing_root)
    
    # Create the report
    with PdfPages(pdf_path) as pdf:
        # Create A4 figure (21cm x 29.7cm)
        fig = plt.figure(figsize=(21/2.54, 29.7/2.54), dpi=150)
        fig.patch.set_facecolor('white')
        
        # Create grid layout: 16 rows, 2 columns
        gs = fig.add_gridspec(
            16, 2,
            left=1.5/2.54/fig.get_size_inches()[0],
            right=1 - 1.5/2.54/fig.get_size_inches()[0],
            bottom=0.5/2.54/fig.get_size_inches()[1],
            top=1 - 1.5/2.54/fig.get_size_inches()[1],
            hspace=0.05,
            wspace=0.1
        )
        
        # === TITLE ===
        _add_title(fig, gs[0, :], folder_name)
        
        # === IMAGE GRID ===
        vessel_types = ['artery', 'vein']
        
        for col, vessel_type in enumerate(vessel_types):
            # Top row: vessel map (4 rows high)
            ax1 = fig.add_subplot(gs[1:5, col])
            ax1.axis('off')
            im1 = _try_load_vessel_image(
                png_dir, mask_dir, hd_png_dir, folder_name, vessel_type,
                report_images=images,
            )
            _show_report_image(ax1, im1)
            
            # Middle row: RI plot (3 rows high)
            ax2 = fig.add_subplot(gs[5:8, col])
            ax2.axis('off')
            im2 = _try_load_ri_image(png_dir, folder_name, vessel_type, report_images=images)
            if im2 is not None:
                _show_report_image(ax2, im2, zoom=1.02, right_pad_fraction=0.04)
            
            # Bottom row: systole indices (3 rows high)
            ax3 = fig.add_subplot(gs[8:11, col])
            ax3.axis('off')
            im3 = _try_load_systole_image(
                png_dir, folder_name, vessel_type, report_images=images
            )
            if im3 is not None:
                _show_report_image(ax3, im3, zoom=1.02, right_pad_fraction=0.04)
        
        # === PARAMETERS SECTION ===
        ax_params = fig.add_subplot(gs[11:16, :])
        ax_params.axis('off')
        _add_parameters_section(ax_params, parameters)
        
        # Save the figure
        pdf.savefig(fig, dpi=300)
        plt.close(fig)
    
    return pdf_path


def _add_title(fig, grid_spec, folder_name: str) -> None:
    """Add title to the report."""
    ax = fig.add_subplot(grid_spec)
    ax.axis('off')
    ax.text(
        0.5, 0.5, folder_name,
        fontsize=16, fontweight='bold',
        ha='center', va='center',
        transform=ax.transAxes,
        wrap=True
    )


def _show_report_image(
    ax,
    image: np.ndarray,
    *,
    zoom: float = 1.04,
    right_pad_fraction: float = 0.0,
) -> None:
    """Draw an image without distorting its pixel aspect ratio."""
    height, width = np.asarray(image).shape[:2]
    ax.imshow(image, aspect="equal")
    ax.set_anchor("C")
    ax.set_adjustable("box")
    ax.margins(0)

    if zoom > 1.0 and width > 0 and height > 0:
        x_center = (width - 1) / 2
        y_center = (height - 1) / 2
        half_width = width / (2 * zoom)
        half_height = height / (2 * zoom)
        right_pad = width * max(right_pad_fraction, 0.0)
        ax.set_xlim(x_center - half_width, x_center + half_width + right_pad)
        ax.set_ylim(y_center + half_height, y_center - half_height)


def _add_parameters_section(ax, parameters: dict[str, Any]) -> None:
    """Add parameters section to the report."""
    param_texts = _format_parameters_for_display(parameters)
    
    num_params = len(param_texts)
    if num_params == 0:
        ax.text(0.02, 0.5, "No parameters extracted from H5 files",
                fontsize=10, va='center')
        return
    
    # Title
    ax.text(0.02, 0.98, 'Computed Parameters:',
            fontsize=14, fontweight='bold', va='top')
    
    if num_params > 6:
        # Two columns
        split_idx = num_params // 2
        col1 = param_texts[:split_idx]
        col2 = param_texts[split_idx:]
        
        y_start = 0.86
        y_step = 0.055
        
        for i, text in enumerate(col1):
            if y_start - i * y_step > 0.05:
                ax.text(0.02, y_start - i * y_step, text,
                        fontsize=9, va='top')
        
        for i, text in enumerate(col2):
            if y_start - i * y_step > 0.05:
                ax.text(0.52, y_start - i * y_step, text,
                        fontsize=9, va='top')
    else:
        # Single column
        y_start = 0.86
        y_step = 0.08
        
        for i, text in enumerate(param_texts):
            if y_start - i * y_step > 0.05:
                ax.text(0.02, y_start - i * y_step, text,
                        fontsize=9, va='top')


def _try_load_vessel_image(
    png_dir, mask_dir, hd_png_dir, folder_name, vessel_type, *,
    report_images: ReportImagePaths | None = None,
):
    """Try multiple patterns for vessel image."""
    possible_paths = []
    registered = (report_images or {}).get(("vessel_map", vessel_type))
    if registered is not None:
        possible_paths.append(Path(registered))
    
    if png_dir:
        possible_paths.append(
            png_dir / f"{folder_name}_{vessel_type}_seg_map_bkg.png"
        )
    
    if mask_dir:
        possible_paths.append(mask_dir / f"{folder_name}_M0_{vessel_type}.png")
        possible_paths.append(mask_dir / f"M0_{vessel_type}.png")
    
    if hd_png_dir:
        possible_paths.append(hd_png_dir / f"{folder_name}_M0.png")
        possible_paths.append(hd_png_dir / "M0.png")
    
    return _load_or_placeholder(possible_paths)


def _try_load_ri_image(
    png_dir, folder_name, vessel_type, *, report_images: ReportImagePaths | None = None
):
    """Load the RI image when the PNG output directory is available."""
    registered = (report_images or {}).get(("ri", vessel_type))
    paths = [Path(registered)] if registered is not None else []
    if png_dir:
        paths.append(png_dir / f"{folder_name}_RI_v_{vessel_type}.png")
    if not paths:
        return None
    return _load_or_placeholder(paths)


def _try_load_systole_image(
    png_dir, folder_name, vessel_type, *, report_images: ReportImagePaths | None = None
):
    """Load the systole image when the PNG output directory is available."""
    registered = (report_images or {}).get(("systole", vessel_type))
    paths = [Path(registered)] if registered is not None else []
    if png_dir:
        paths.append(png_dir / f"{folder_name}_find_systoles_indices_{vessel_type}.png")
    if not paths:
        return None
    return _load_or_placeholder(paths)


def _load_or_placeholder(paths: list[Path]) -> np.ndarray:
    """Load first existing image, or return a white placeholder."""
    for path in paths:
        if path.exists():
            try:
                with Image.open(path) as image:
                    return np.asarray(image.convert("RGB"))
            except Exception:
                continue
    
    # Create white placeholder
    placeholder = np.ones((200, 200, 3), dtype=np.uint8) * 255
    return placeholder
