"""Select measurements and quality indicators for one PDF workflow."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from input_output.readers.nodes import is_dataset
import numpy as np

from input_output.schema import EyeFlowOutputPaths
from input_output.schema.provenance import read_dataset_attrs


@dataclass(frozen=True)
class VolumeRateReportData:
    beat_curves: np.ndarray
    phase_curve: np.ndarray
    overall_value: float | None
    unit: str


@dataclass(frozen=True)
class ReportData:
    parameters: dict[str, Any] = field(default_factory=dict)
    velocity_per_beat: dict[str, np.ndarray] = field(default_factory=dict)
    velocity_units: dict[str, str] = field(default_factory=dict)
    mask_volume_rates: dict[str, VolumeRateReportData] = field(default_factory=dict)


def extract_report_data(h5file, *, include_bvr: bool = True) -> ReportData:
    """Select velocity and BVR inputs from the moments output view."""
    params: dict[str, Any] = {}
    schema = EyeFlowOutputPaths.active()
    velocities, units = _extract_velocity_metrics(h5file, schema, params)
    _extract_heart_rate(h5file, schema, params)
    _extract_waveform_metrics(h5file, schema, params, velocities)
    _extract_attributes(h5file, params)
    return ReportData(
        parameters=params,
        velocity_per_beat=velocities,
        velocity_units=units,
        mask_volume_rates=(
            _extract_mask_volume_rates(h5file, schema) if include_bvr else {}
        ),
    )


def _dataset(h5file, path):
    node = h5file.get(path)
    return node if is_dataset(node) and node.size else None


def _extract_velocity_metrics(h5file, schema, params) -> tuple[dict[str, np.ndarray], dict[str, str]]:
    velocities: dict[str, np.ndarray] = {}
    units: dict[str, str] = {}
    for vessel, path in (
        ("artery", schema.artery_per_beat_filtered.velocity_signal),
        ("vein", schema.vein_per_beat_filtered.velocity_signal),
    ):
        dataset = _dataset(h5file, path)
        if dataset is None:
            continue
        values = np.abs(np.asarray(dataset, dtype=float))
        if values.ndim != 2 or not all(values.shape):
            continue
        unit = _velocity_dataset_unit(h5file, dataset)
        finite = values[np.isfinite(values)]
        velocities[vessel] = values
        units[vessel] = unit
        if not finite.size:
            continue
        label = "Arterial" if vessel == "artery" else "Venous"
        for name, operation in (("Average", np.mean), ("Max", np.max), ("Min", np.min)):
            params[f"{name}_{label}_Velocity"] = {
                "value": float(operation(finite)), "unit": unit
            }
    return velocities, units


def _velocity_dataset_unit(h5file, dataset) -> str:
    """Use dataset provenance, with the output root as fallback."""
    attrs = read_dataset_attrs(dataset)
    value = attrs.get(
        "unit", attrs.get("velocity_unit", h5file.attrs.get("velocity_unit", "mm/s"))
    )
    if isinstance(value, bytes):
        value = value.decode("utf-8")
    return str(value)


def _extract_heart_rate(h5file, schema, params) -> None:
    heart_rate = _dataset(h5file, schema.cardiac_cycle.spectral_heart_rate_bpm)
    if heart_rate is not None:
        finite = np.asarray(heart_rate, dtype=float)
        finite = finite[np.isfinite(finite)]
        if finite.size:
            params["heart_beat"] = {"value": float(np.mean(finite)), "unit": "bpm"}


def _extract_waveform_metrics(h5file, schema, params, velocities) -> None:
    for vessel, prefix in (("artery", "A"), ("vein", "V")):
        if vessel not in velocities:
            continue
        root = f"{schema.waveform_shape_metrics_root}/{vessel}/global/raw"
        for metric in ("RI", "PI"):
            dataset = _dataset(h5file, f"{root}/{metric}")
            if dataset is None:
                continue
            values = np.asarray(dataset, dtype=float)
            finite = values[np.isfinite(values)]
            if finite.size:
                params[f"{prefix}{metric}"] = {"value": float(np.mean(finite)), "unit": ""}

def _extract_attributes(h5file, params) -> None:
    for name in (
        "UnixTimestampFirst", "UnixTimestampLast", "SystoleDuration", "DiastoleDuration"
    ):
        if name in h5file.attrs:
            value = h5file.attrs[name]
            if isinstance(value, (int, float, str)):
                params[name] = {"value": value, "unit": ""}


def _extract_mask_volume_rates(
    h5file, schema,
) -> dict[str, VolumeRateReportData]:
    rates = {}
    for vessel, paths in (
        ("artery", schema.blood_volume_rate.artery),
        ("vein", schema.blood_volume_rate.vein),
    ):
        total = _dataset(h5file, paths.total_masked_edges)
        if total is not None and total.ndim == 2 and all(total.shape):
            rates[vessel] = _volume_rate_data(
                np.asarray(total, dtype=float), total,
            )
    return rates


def _volume_rate_data(beat_curves: np.ndarray, dataset) -> VolumeRateReportData:
    phase_curve = np.ma.median(
        np.ma.masked_invalid(beat_curves), axis=1,
    ).filled(np.nan)
    finite = phase_curve[np.isfinite(phase_curve)]
    return VolumeRateReportData(
        beat_curves=np.asarray(beat_curves, dtype=float),
        phase_curve=np.asarray(phase_curve, dtype=float),
        overall_value=float(np.median(finite)) if finite.size else None,
        unit=str(read_dataset_attrs(dataset).get("unit", "mm^3/s")),
    )


