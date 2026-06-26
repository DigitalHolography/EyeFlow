"""A4 PDF report generator mimicking MATLAB generateA4Report."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from PIL import Image

from input_output.writers.h5 import open_h5


def generate_a4_report(
    output_h5_paths: list[Path],
    output_dir: Path,
    settings: dict[str, Any] | None = None,
    png_dir: Path | None = None,
    hd_png_dir: Path | None = None,
    mask_dir: Path | None = None,
    folder_name: str | None = None,
    include_errors: bool = False,
    include_qc_scores: bool = False,
) -> Path:
    """
    Generate an A4 PDF report mimicking the MATLAB generateA4Report function.
    
    This creates a comprehensive report with:
    - Title with folder name
    - Image grid (vessel map, RI plot, systole indices for artery and vein)
    - Parameters table with computed metrics
    - Optional error section
    - Optional quality control scores
    
    Args:
        output_h5_paths: List of output H5 files from processing
        output_dir: Directory to save the PDF
        settings: Optional settings used for the run
        png_dir: Directory containing PNG outputs
        hd_png_dir: Directory containing HoloDoppler M0 PNGs
        mask_dir: Directory containing mask PNGs
        folder_name: Name of the folder (for title)
        include_errors: Whether to include error section
        include_qc_scores: Whether to include quality control scores
    
    Returns:
        Path to the generated PDF file
    """
    if not output_h5_paths:
        raise ValueError("No output H5 files provided for report generation")
    
    # Determine folder name if not provided
    if folder_name is None:
        folder_name = output_h5_paths[0].parent.name
    
    # Create output path
    pdf_path = output_dir / f"{folder_name}_report.pdf"
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Extract parameters from H5 files
    parameters = _extract_parameters_from_h5(output_h5_paths)
    
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
            im1 = _try_load_vessel_image(png_dir, mask_dir, hd_png_dir, folder_name, vessel_type)
            if im1 is not None:
                _show_report_image(ax1, im1)
            
            # Middle row: RI plot (3 rows high)
            ax2 = fig.add_subplot(gs[5:8, col])
            ax2.axis('off')
            im2 = _try_load_ri_image(png_dir, folder_name, vessel_type)
            if im2 is not None:
                _show_report_image(ax2, im2, zoom=1.02, right_pad_fraction=0.04)
            
            # Bottom row: systole indices (3 rows high)
            ax3 = fig.add_subplot(gs[8:11, col])
            ax3.axis('off')
            im3 = _try_load_systole_image(png_dir, folder_name, vessel_type)
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


def _try_load_vessel_image(png_dir, mask_dir, hd_png_dir, folder_name, vessel_type):
    """Try multiple patterns for vessel image."""
    possible_paths = []
    
    if png_dir:
        patterns = [
            f"{folder_name}_{vessel_type}_seg_map_bkg.png",
            f"{folder_name}_vessel_map_{vessel_type}.png",
            f"{folder_name}_{vessel_type}_map.png",
            f"{vessel_type}_seg_map_bkg.png",
            f"{vessel_type}_map.png",
            f"vessel_map_{vessel_type}.png",
        ]
        for pattern in patterns:
            possible_paths.append(png_dir / pattern)
    
    if mask_dir:
        possible_paths.append(mask_dir / f"{folder_name}_M0_{vessel_type}.png")
        possible_paths.append(mask_dir / f"M0_{vessel_type}.png")
    
    if hd_png_dir:
        possible_paths.append(hd_png_dir / f"{folder_name}_M0.png")
        possible_paths.append(hd_png_dir / "M0.png")
    
    return _load_or_placeholder(possible_paths)


def _try_load_ri_image(png_dir, folder_name, vessel_type):
    """Try multiple patterns for RI image."""
    if not png_dir:
        return None
    
    patterns = [
        f"{folder_name}_RI_v_{vessel_type}.png",
        f"RI_{vessel_type}.png",
        f"{vessel_type}_RI.png",
    ]
    paths = [png_dir / pattern for pattern in patterns]
    return _load_or_placeholder(paths)


def _try_load_systole_image(png_dir, folder_name, vessel_type):
    """Try multiple patterns for systole image."""
    if not png_dir:
        return None
    
    patterns = [
        f"{folder_name}_find_systoles_indices_{vessel_type}.png",
        f"find_systoles_{vessel_type}.png",
        f"systoles_{vessel_type}.png",
    ]
    paths = [png_dir / pattern for pattern in patterns]
    return _load_or_placeholder(paths)


def _load_or_placeholder(paths: list[Path | None]) -> np.ndarray | None:
    """Load first existing image, or return a white placeholder."""
    for path in paths:
        if path is not None and path.exists():
            try:
                img = Image.open(path)
                if img.mode != 'RGB':
                    img = img.convert('RGB')
                return np.array(img)
            except Exception:
                continue
    
    # Create white placeholder
    placeholder = np.ones((200, 200, 3), dtype=np.uint8) * 255
    return placeholder


def _extract_parameters_from_h5(h5_paths: list[Path]) -> dict[str, Any]:
    """Extract parameters from H5 output files."""
    params = {}
    
    if not h5_paths:
        return params
    
    try:
        with open_h5(h5_paths[0], 'r') as f:
            # Velocity metrics
            _extract_velocity_metrics(f, params)
            
            # Heart rate
            _extract_heart_rate(f, params)
            
            # Metrics from waveform_shape_metrics
            _extract_waveform_metrics(f, params)
            
            # Attributes
            _extract_attributes(f, params)
    
    except Exception as e:
        print(f"Warning: Could not extract parameters from H5: {e}")
    
    return params


def _extract_velocity_metrics(f: h5py.File, params: dict[str, Any]) -> None:
    """Extract velocity metrics from H5 file."""
    # Try different possible paths
    artery_paths = [
        "Artery/VelocityPerBeat/VelocitySignalPerBeat/value",
        "Artery/velocity/signal/value",
        "artery/velocity/perbeat/signal/value",
    ]
    
    for path in artery_paths:
        dataset = f.get(path)
        if dataset is not None:
            arr = np.asarray(dataset)
            if arr.size > 0:
                params["Average_Arterial_Velocity"] = {"value": float(np.mean(arr)), "unit": "mm/s"}
                params["Max_Arterial_Velocity"] = {"value": float(np.max(arr)), "unit": "mm/s"}
                params["Min_Arterial_Velocity"] = {"value": float(np.min(arr)), "unit": "mm/s"}
                break
    
    vein_paths = [
        "Vein/VelocityPerBeat/VelocitySignalPerBeat/value",
        "Vein/velocity/signal/value",
        "vein/velocity/perbeat/signal/value",
    ]
    
    for path in vein_paths:
        dataset = f.get(path)
        if dataset is not None:
            arr = np.asarray(dataset)
            if arr.size > 0:
                params["Average_Venous_Velocity"] = {"value": float(np.mean(arr)), "unit": "mm/s"}
                params["Max_Venous_Velocity"] = {"value": float(np.max(arr)), "unit": "mm/s"}
                params["Min_Venous_Velocity"] = {"value": float(np.min(arr)), "unit": "mm/s"}
                break


def _extract_heart_rate(f: h5py.File, params: dict[str, Any]) -> None:
    """Extract heart rate from H5 file."""
    beat_paths = [
        "Artery/VelocityPerBeat/beatPeriodSeconds/value",
        "perbeat/beat_period_seconds/value",
    ]
    
    for path in beat_paths:
        dataset = f.get(path)
        if dataset is not None:
            arr = np.asarray(dataset)
            if arr.size > 0:
                params["heart_beat"] = {"value": 60.0 / float(np.mean(arr)), "unit": "bpm"}
                break


def _extract_waveform_metrics(f: h5py.File, params: dict[str, Any]) -> None:
    """Extract metrics from waveform_shape_metrics output."""
    # Look for Metrics group
    metrics = f.get("Metrics")
    if not metrics:
        return
    
    # Map to MATLAB parameter names
    mapping = {
        "artery_resistivity_index": "ARI",
        "vein_resistivity_index": "VRI",
        "artery_pulsatility_index": "API",
        "vein_pulsatility_index": "VPI",
        "dicrotic_notch_visibility": "DicroticNotchVisibility",
        "artery_flow_rate_mean": "Average_Arterial_Volume_Rate",
        "vein_flow_rate_mean": "Average_Venous_Volume_Rate",
        "systole_duration": "SystoleDuration",
        "diastole_duration": "DiastoleDuration",
        "artery_time_peak_to_descent": "TimePeakToDescent",
        "vein_time_to_peak_from_min": "VeinTimeToPeakFromMin",
    }
    
    # Iterate through metrics
    for metric_name, dataset in metrics.items():
        if isinstance(dataset, h5py.Dataset):
            arr = np.asarray(dataset)
            if arr.size > 0:
                matlab_name = mapping.get(metric_name)
                if matlab_name:
                    params[matlab_name] = {
                        "value": float(arr[0] if arr.size == 1 else arr),
                        "unit": ""
                    }


def _extract_attributes(f: h5py.File, params: dict[str, Any]) -> None:
    """Extract parameters from H5 attributes."""
    attr_mapping = {
        "UnixTimestampFirst": "UnixTimestampFirst",
        "UnixTimestampLast": "UnixTimestampLast",
        "SystoleDuration": "SystoleDuration",
        "DiastoleDuration": "DiastoleDuration",
    }
    
    for attr_name, param_name in attr_mapping.items():
        if attr_name in f.attrs:
            value = f.attrs[attr_name]
            if isinstance(value, (int, float, str)):
                params[param_name] = {"value": value, "unit": ""}


def _format_parameters_for_display(parameters: dict[str, Any]) -> list[str]:
    """Format parameters for display in the report."""
    formats = {
        "UnixTimestampFirst": "First Unix Timestamp = %d",
        "UnixTimestampLast": "Last Unix Timestamp = %d",
        "heart_beat": "HR = %.1f",
        "Average_Arterial_Velocity": "Avg Arterial Velocity = %.2f",
        "Max_Arterial_Velocity": "Max Arterial Velocity = %.2f",
        "Min_Arterial_Velocity": "Min Arterial Velocity = %.2f",
        "Average_Venous_Velocity": "Avg Venous Velocity = %.2f",
        "Max_Venous_Velocity": "Max Venous Velocity = %.2f",
        "Min_Venous_Velocity": "Min Venous Velocity = %.2f",
        "TimePeakToDescent": "Time Peak to Descent = %.2f",
        "VeinTimeToPeakFromMin": "Time to Peak from Min Vein = %.2f",
        "DicroticNotchVisibility": "Dicrotic Notch Visibility = %.0f",
        "Average_Arterial_Volume_Rate": "Avg Arterial Volume Rate = %.2f",
        "Average_Venous_Volume_Rate": "Avg Venous Volume Rate = %.2f",
        "SystoleDuration": "Systole Duration = %.2f",
        "DiastoleDuration": "Diastole Duration = %.2f",
        "ARI": "Arterial Resistivity Index = %.2f",
        "VRI": "Venous Resistivity Index = %.2f",
        "API": "Arterial Pulsatility Index = %.2f",
        "VPI": "Venous Pulsatility Index = %.2f",
    }
    
    formatted = []
    
    # Order parameters following MATLAB convention
    order = [
        "UnixTimestampFirst", "UnixTimestampLast",
        "heart_beat",
        "Average_Arterial_Velocity", "Max_Arterial_Velocity", "Min_Arterial_Velocity",
        "Average_Venous_Velocity", "Max_Venous_Velocity", "Min_Venous_Velocity",
        "ARI", "VRI", "API", "VPI",
        "TimePeakToDescent", "VeinTimeToPeakFromMin",
        "DicroticNotchVisibility",
        "Average_Arterial_Volume_Rate", "Average_Venous_Volume_Rate",
        "SystoleDuration", "DiastoleDuration",
    ]
    
    # Add parameters in order
    for key in order:
        if key in parameters:
            value = parameters[key]
            fmt = formats.get(key, "%s = %g")
            formatted.append(_format_value(value, fmt))
    
    # Add any remaining parameters
    for key, value in parameters.items():
        if key not in order:
            fmt = formats.get(key, "%s = %g")
            formatted.append(_format_value(value, fmt))
    
    return formatted


def _format_value(value: dict | Any, fmt: str) -> str:
    """Format a parameter value using the specified format."""
    if isinstance(value, dict):
        val = value.get("value", 0)
        unit = value.get("unit", "")
        try:
            if isinstance(val, (int, float)):
                formatted = fmt % val
            else:
                formatted = fmt.replace("%s", str(val)).replace("%g", str(val))
        except (TypeError, ValueError):
            formatted = f"{val}"
        return f"{formatted} {unit}".strip()
    return str(value)
