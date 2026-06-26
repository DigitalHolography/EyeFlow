"""Helper functions for report generation."""

from pathlib import Path
from typing import Any

import h5py
import numpy as np
from PIL import Image

from input_output.writers.h5 import open_h5


def load_or_placeholder(paths: list[Path | None]) -> np.ndarray | None:
    """
    Load first existing image, or return a white placeholder.
    
    Args:
        paths: List of potential image paths
    
    Returns:
        Image as numpy array, or placeholder if none exist
    """
    for path in paths:
        if path is not None and path.exists():
            try:
                return np.array(Image.open(path))
            except Exception:
                continue
    
    # Create white placeholder (200x200 RGB)
    placeholder = np.ones((200, 200, 3), dtype=np.uint8) * 255
    return placeholder


def extract_parameters_from_h5(h5_paths: list[Path]) -> dict[str, Any]:
    """
    Extract parameters from H5 output files.
    
    Maps MATLAB parameter names to values extracted from H5 datasets.
    
    Args:
        h5_paths: List of H5 file paths
    
    Returns:
        Dictionary of parameter names to values
    """
    params = {}
    
    if not h5_paths:
        return params
    
    # Try to load from the first H5 file
    try:
        with open_h5(h5_paths[0], 'r') as f:
            # === Velocity metrics ===
            _extract_velocity_metrics(f, params)
            
            # === Heart rate ===
            _extract_heart_rate(f, params)
            
            # === Metrics group ===
            _extract_metrics_group(f, params)
            
            # === Attributes ===
            _extract_attributes(f, params)
    
    except Exception as e:
        print(f"Warning: Could not extract parameters from H5: {e}")
    
    return params


def _extract_velocity_metrics(f: h5py.File, params: dict[str, Any]) -> None:
    """Extract velocity metrics from H5 file."""
    # Artery velocity
    artery_velocity = f.get("Artery/VelocityPerBeat/VelocitySignalPerBeat/value")
    if artery_velocity is not None:
        arr = np.asarray(artery_velocity)
        if arr.size > 0:
            params["Average_Arterial_Velocity"] = {"value": float(np.mean(arr)), "unit": "mm/s"}
            params["Max_Arterial_Velocity"] = {"value": float(np.max(arr)), "unit": "mm/s"}
            params["Min_Arterial_Velocity"] = {"value": float(np.min(arr)), "unit": "mm/s"}
    
    # Vein velocity
    vein_velocity = f.get("Vein/VelocityPerBeat/VelocitySignalPerBeat/value")
    if vein_velocity is not None:
        arr = np.asarray(vein_velocity)
        if arr.size > 0:
            params["Average_Venous_Velocity"] = {"value": float(np.mean(arr)), "unit": "mm/s"}
            params["Max_Venous_Velocity"] = {"value": float(np.max(arr)), "unit": "mm/s"}
            params["Min_Venous_Velocity"] = {"value": float(np.min(arr)), "unit": "mm/s"}


def _extract_heart_rate(f: h5py.File, params: dict[str, Any]) -> None:
    """Extract heart rate from H5 file."""
    beat_period = f.get("Artery/VelocityPerBeat/beatPeriodSeconds/value")
    if beat_period is not None:
        arr = np.asarray(beat_period)
        if arr.size > 0:
            params["heart_beat"] = {"value": 60.0 / float(np.mean(arr)), "unit": "bpm"}


def _extract_metrics_group(f: h5py.File, params: dict[str, Any]) -> None:
    """Extract metrics from the Metrics group."""
    metrics = f.get("Metrics")
    if not metrics:
        return
    
    # Map internal names to MATLAB parameter names
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


def format_parameters_for_display(parameters: dict[str, Any]) -> list[str]:
    """
    Format parameters for display in the report.
    
    Args:
        parameters: Dictionary of parameter names to values
    
    Returns:
        List of formatted parameter strings
    """
    # Define display formats for each parameter
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
    
    # Add any remaining parameters not in the order list
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