"""Read and format HDF5 measurements for the A4 report."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import h5py
import numpy as np

from input_output.schema import EyeFlowOutputPaths
from input_output.schema.provenance import read_dataset_attrs
from input_output.writers.h5 import open_h5


_VELOCITY_LEGACY_PATHS = {
    "artery": (
        "Artery/VelocityPerBeat/VelocitySignalPerBeat/value",
        "Artery/velocity/signal/value",
        "artery/velocity/perbeat/signal/value",
    ),
    "vein": (
        "Vein/VelocityPerBeat/VelocitySignalPerBeat/value",
        "Vein/velocity/signal/value",
        "vein/velocity/perbeat/signal/value",
    ),
}
_BEAT_LEGACY_PATHS = (
    "Artery/VelocityPerBeat/beatPeriodSeconds/value",
    "perbeat/beat_period_seconds/value",
)
_LEGACY_METRIC_NAMES = {
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
_FORMATS = {
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
_DISPLAY_ORDER = (
    "UnixTimestampFirst", "UnixTimestampLast", "heart_beat",
    "Average_Arterial_Velocity", "Max_Arterial_Velocity", "Min_Arterial_Velocity",
    "Average_Venous_Velocity", "Max_Venous_Velocity", "Min_Venous_Velocity",
    "ARI", "VRI", "API", "VPI", "TimePeakToDescent", "VeinTimeToPeakFromMin",
    "DicroticNotchVisibility", "Average_Arterial_Volume_Rate",
    "Average_Venous_Volume_Rate", "SystoleDuration", "DiastoleDuration",
)


def extract_parameters_from_h5(
    h5_path: Path, *, processing_root: str | None = None
) -> dict[str, Any]:
    """Read current and recognized legacy paths without changing report values."""
    params: dict[str, Any] = {}
    try:
        with open_h5(h5_path, "r") as h5file:
            if processing_root is not None:
                from input_output.writers.h5_output import PipelineH5Output

                h5file = PipelineH5Output(h5file, processing_root=processing_root)
            schema = EyeFlowOutputPaths.active()
            _extract_velocity_metrics(h5file, schema, params)
            _extract_heart_rate(h5file, schema, params)
            _extract_waveform_metrics(h5file, schema, params)
            _extract_attributes(h5file, params)
    except Exception as exc:
        print(f"Warning: Could not extract parameters from H5: {exc}")
    return params


def _first_nonempty_dataset(h5file, paths) -> h5py.Dataset | None:
    for path in paths:
        dataset = h5file.get(path)
        if isinstance(dataset, h5py.Dataset):
            if dataset.size:
                return dataset
    return None


def _extract_velocity_metrics(h5file, schema, params) -> None:
    for vessel, paths in (
        ("artery", (
            schema.artery_per_beat_filtered.velocity_signal_band_limited,
            schema.artery_per_beat.velocity_signal,
            *_VELOCITY_LEGACY_PATHS["artery"],
        )),
        ("vein", (
            schema.vein_per_beat_filtered.velocity_signal_band_limited,
            schema.vein_per_beat.velocity_signal,
            *_VELOCITY_LEGACY_PATHS["vein"],
        )),
    ):
        dataset = _first_nonempty_dataset(h5file, paths)
        if dataset is None:
            continue
        values = np.asarray(dataset)
        unit = _velocity_dataset_unit(h5file, dataset)
        label = "Arterial" if vessel == "artery" else "Venous"
        for name, operation in (("Average", np.mean), ("Max", np.max), ("Min", np.min)):
            params[f"{name}_{label}_Velocity"] = {
                "value": float(operation(values)), "unit": unit
            }


def _velocity_dataset_unit(h5file, dataset: h5py.Dataset) -> str:
    """Use dataset provenance, with the output root as fallback."""
    attrs = read_dataset_attrs(dataset)
    value = attrs.get(
        "unit", attrs.get("velocity_unit", h5file.attrs.get("velocity_unit", "mm/s"))
    )
    if isinstance(value, bytes):
        value = value.decode("utf-8")
    return str(value)


def _extract_heart_rate(h5file, schema, params) -> None:
    heart_rate = _first_nonempty_dataset(h5file, (schema.cardiac_cycle.spectral_heart_rate_bpm,))
    if heart_rate is not None:
        params["heart_beat"] = {"value": float(np.nanmean(np.asarray(heart_rate))), "unit": "bpm"}
        return
    periods = _first_nonempty_dataset(
        h5file, (schema.cardiac_cycle.systolic_cycle_duration_seconds, *_BEAT_LEGACY_PATHS)
    )
    if periods is not None:
        params["heart_beat"] = {"value": 60.0 / float(np.mean(np.asarray(periods))), "unit": "bpm"}


def _extract_waveform_metrics(h5file, schema, params) -> None:
    for vessel, prefix in (("artery", "A"), ("vein", "V")):
        root = f"{schema.waveform_shape_metrics_root}/{vessel}/global/raw"
        for metric in ("RI", "PI"):
            dataset = h5file.get(f"{root}/{metric}")
            if not isinstance(dataset, h5py.Dataset):
                continue
            values = np.asarray(dataset, dtype=float)
            finite = values[np.isfinite(values)]
            if finite.size:
                params[f"{prefix}{metric}"] = {"value": float(np.mean(finite)), "unit": ""}

    legacy = h5file.get("Metrics")
    if not isinstance(legacy, h5py.Group):
        return
    for name, dataset in legacy.items():
        target = _LEGACY_METRIC_NAMES.get(name)
        if target and isinstance(dataset, h5py.Dataset):
            values = np.asarray(dataset)
            if values.size:
                params[target] = {"value": float(np.nanmean(values)), "unit": ""}


def _extract_attributes(h5file, params) -> None:
    for name in (
        "UnixTimestampFirst", "UnixTimestampLast", "SystoleDuration", "DiastoleDuration"
    ):
        if name in h5file.attrs:
            value = h5file.attrs[name]
            if isinstance(value, (int, float, str)):
                params[name] = {"value": value, "unit": ""}


def format_parameters_for_display(parameters: dict[str, Any]) -> list[str]:
    """Format known parameters first, leaving other extracted values visible."""
    ordered = [name for name in _DISPLAY_ORDER if name in parameters]
    ordered.extend(name for name in parameters if name not in _DISPLAY_ORDER)
    return [_format_value(parameters[name], _FORMATS.get(name, "%s = %g")) for name in ordered]


def _format_value(value: dict | Any, fmt: str) -> str:
    if not isinstance(value, dict):
        return str(value)
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


_extract_parameters_from_h5 = extract_parameters_from_h5
_format_parameters_for_display = format_parameters_for_display
