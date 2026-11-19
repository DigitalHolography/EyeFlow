"""Select measurements and quality indicators for the moments reports."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from calculations.blood_flow_velocity.signal_analysis.waveform.cycles import (
    cycle_extrema,
)
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
class SystoleReportData:
    velocity: np.ndarray
    systolic_boundaries: np.ndarray
    maxima: np.ndarray
    minima: np.ndarray
    sample_interval_seconds: float | None
    unit: str


@dataclass(frozen=True)
class VesselGeometryReportData:
    """Compact geometry statistics used by the second report page."""

    sampled_branches: int
    valid_sampling_regions: int
    mean_lumen_diameter_microns: float | None
    median_lumen_diameter_microns: float | None
    lumen_diameter_standard_deviation_microns: float | None = None


@dataclass(frozen=True)
class ReportData:
    parameters: dict[str, Any] = field(default_factory=dict)
    velocity_per_beat: dict[str, np.ndarray] = field(default_factory=dict)
    velocity_per_beat_unfiltered: dict[str, np.ndarray] = field(default_factory=dict)
    velocity_units: dict[str, str] = field(default_factory=dict)
    systole_signals: dict[str, SystoleReportData] = field(default_factory=dict)
    mask_volume_rates: dict[str, VolumeRateReportData] = field(default_factory=dict)
    vessel_geometry: dict[str, VesselGeometryReportData] = field(default_factory=dict)


def extract_report_data(h5file, *, include_bvr: bool = True) -> ReportData:
    """Select velocity and BVR inputs from the moments output view."""
    params: dict[str, Any] = {}
    schema = EyeFlowOutputPaths.active()
    velocities, units = _extract_velocity_metrics(h5file, schema, params)
    unfiltered_velocities, _ = _extract_velocity_arrays(
        h5file,
        (
            ("artery", schema.artery_per_beat.velocity_signal),
            ("vein", schema.vein_per_beat.velocity_signal),
        ),
    )
    _extract_heart_rate(h5file, schema, params)
    _extract_waveform_metrics(h5file, schema, params, velocities)
    _extract_attributes(h5file, params)
    return ReportData(
        parameters=params,
        velocity_per_beat=velocities,
        velocity_per_beat_unfiltered=unfiltered_velocities,
        velocity_units=units,
        systole_signals=_extract_systole_signals(h5file, schema),
        mask_volume_rates=(
            _extract_mask_volume_rates(h5file, schema) if include_bvr else {}
        ),
        vessel_geometry=_extract_vessel_geometry(h5file, schema),
    )
def _dataset(h5file, path):
    node = h5file.get(path)
    return node if is_dataset(node) and node.size else None


def _extract_velocity_metrics(
    h5file, schema, params,
) -> tuple[dict[str, np.ndarray], dict[str, str]]:
    velocities, units = _extract_velocity_arrays(
        h5file,
        (
            ("artery", schema.artery_per_beat_filtered.velocity_signal),
            ("vein", schema.vein_per_beat_filtered.velocity_signal),
        ),
    )
    for vessel, values in velocities.items():
        unit = units[vessel]
        finite = values[np.isfinite(values)]
        if not finite.size:
            continue
        label = "Arterial" if vessel == "artery" else "Venous"
        operations = (("Average", np.mean), ("Max", np.max), ("Min", np.min))
        for name, operation in operations:
            params[f"{name}_{label}_Velocity"] = {
                "value": float(operation(finite)), "unit": unit
            }
    return velocities, units


def _extract_velocity_arrays(
    h5file, vessel_paths,
) -> tuple[dict[str, np.ndarray], dict[str, str]]:
    velocities: dict[str, np.ndarray] = {}
    units: dict[str, str] = {}
    for vessel, path in vessel_paths:
        dataset = _dataset(h5file, path)
        if dataset is None:
            continue
        values = np.abs(np.asarray(dataset, dtype=float))
        if values.ndim != 2 or not all(values.shape):
            continue
        unit = _velocity_dataset_unit(h5file, dataset)
        velocities[vessel] = values
        units[vessel] = unit
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


def _extract_systole_signals(
    h5file, schema,
) -> dict[str, SystoleReportData]:
    boundaries_dataset = _dataset(
        h5file, schema.cardiac_cycle.systolic_peak_frame_indices
    )
    if boundaries_dataset is None:
        return {}

    raw_boundaries = np.asarray(boundaries_dataset, dtype=float).reshape(-1)
    raw_boundaries = raw_boundaries[np.isfinite(raw_boundaries)]
    if not raw_boundaries.size:
        return {}

    datasets = {}
    for vessel, path in (
        ("artery", schema.analysis.retinal_artery_velocity_signal_filtered),
        ("vein", schema.analysis.retinal_vein_velocity_signal_filtered),
    ):
        dataset = _dataset(h5file, path)
        if dataset is None:
            continue
        values = np.asarray(dataset, dtype=float)
        if values.ndim == 1 and values.size:
            datasets[vessel] = (values, dataset)
    if not datasets:
        return {}

    reference = (
        datasets["artery"][0]
        if "artery" in datasets
        else next(iter(datasets.values()))[0]
    )
    boundaries = raw_boundaries.astype(np.int64)
    boundaries = np.unique(boundaries[(boundaries >= 0) & (boundaries < reference.size)])
    if not boundaries.size:
        return {}
    maxima, minima = cycle_extrema(reference, boundaries)
    sample_interval_seconds = _positive_scalar(
        _dataset(h5file, f"{schema.meta_root}/FrameIntervalSeconds/value")
    )

    result = {}
    for vessel, (values, dataset) in datasets.items():
        result[vessel] = SystoleReportData(
            velocity=values,
            systolic_boundaries=boundaries,
            maxima=maxima,
            minima=minima,
            sample_interval_seconds=sample_interval_seconds,
            unit=_velocity_dataset_unit(h5file, dataset),
        )
    return result


def _positive_scalar(dataset) -> float | None:
    if dataset is None:
        return None
    values = np.asarray(dataset, dtype=float)
    if values.size != 1:
        return None
    value = float(values.reshape(()))
    return value if np.isfinite(value) and value > 0.0 else None


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


def _extract_vessel_geometry(
    h5file, schema,
) -> dict[str, VesselGeometryReportData]:
    pitch_dataset = _dataset(h5file, schema.segmentation.pixel_pitch_m)
    pitch_m = None
    if pitch_dataset is not None:
        pitch = np.asarray(pitch_dataset, dtype=float)
        if pitch.size == 1:
            candidate = float(pitch.reshape(()))
            if np.isfinite(candidate) and candidate > 0.0:
                pitch_m = candidate

    geometry: dict[str, VesselGeometryReportData] = {}
    for vessel, paths in (
        ("artery", schema.segmentation.artery),
        ("vein", schema.segmentation.vein),
    ):
        dataset = _dataset(h5file, paths.lumen_diameter)
        if dataset is None:
            continue
        values = np.asarray(dataset, dtype=float)
        if values.ndim != 2 or not all(values.shape):
            continue
        finite_mask = np.isfinite(values)
        finite_microns = (
            np.asarray([], dtype=float)
            if pitch_m is None
            else values[finite_mask] * pitch_m * 1e6
        )
        geometry[vessel] = VesselGeometryReportData(
            sampled_branches=int(np.count_nonzero(np.any(finite_mask, axis=1))),
            valid_sampling_regions=int(np.count_nonzero(finite_mask)),
            mean_lumen_diameter_microns=(
                float(np.mean(finite_microns)) if finite_microns.size else None
            ),
            median_lumen_diameter_microns=(
                float(np.median(finite_microns)) if finite_microns.size else None
            ),
            lumen_diameter_standard_deviation_microns=(
                float(np.std(finite_microns)) if finite_microns.size else None
            ),
        )
    return geometry


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
