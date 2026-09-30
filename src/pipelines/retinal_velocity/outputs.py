"""HDF5 output packing for base retinal-velocity products."""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np

from input_output.schema import EyeFlowOutputPaths

from .models import RetinalVelocity
from .semantics import velocity_dataset_attrs


def pack_retinal_velocity_outputs(
    velocity: RetinalVelocity,
    output_paths: EyeFlowOutputPaths | str | None = None,
) -> dict[str, object]:
    """Pack average maps and cardiac-cycle timing from the core result."""

    schema = _resolve_output_paths(output_paths)
    analysis_paths = schema.analysis
    cycle_paths = schema.cardiac_cycle
    spectral = velocity.cardiac_cycle.spectral
    velocity_attrs = velocity_dataset_attrs(velocity)
    frequency_attrs = {
        key: value
        for key, value in velocity_attrs.items()
        if key not in {"unit", "velocity_quantity"}
    }
    frequency_attrs.update({"unit": "Hz", "quantity": "rms_frequency"})
    outputs = {
        analysis_paths.velocity_map_avg: metric_value(
            velocity.velocity_average,
            dim_desc=("y", "x"),
            attrs=velocity_attrs,
        ),
        analysis_paths.fRMS_avg: metric_value(
            velocity.frms_average,
            dim_desc=("y", "x"),
            attrs=frequency_attrs,
        ),
        analysis_paths.fRMS_bkg_avg: metric_value(
            velocity.frms_background_average,
            dim_desc=("y", "x"),
            attrs=frequency_attrs,
        ),
        analysis_paths.delta_fRMS_avg: metric_value(
            velocity.delta_frms_average,
            dim_desc=("y", "x"),
            attrs=frequency_attrs,
        ),
        cycle_paths.systolic_peak_frame_indices: metric_value(
            velocity.cycle_boundary_indexes,
        ),
        cycle_paths.systolic_cycle_duration_seconds: metric_value(
            velocity.cycle_durations_seconds,
            unit="s",
        ),
        cycle_paths.spectral_fundamental_frequency_hz: metric_value(
            spectral.fundamental_hz,
            unit="Hz",
        ),
        cycle_paths.spectral_heart_rate_bpm: metric_value(
            spectral.heart_rate_bpm,
            unit="bpm",
        ),
        cycle_paths.spectral_heart_rate_standard_error_bpm: metric_value(
            spectral.heart_rate_ste_bpm,
            unit="bpm",
        ),
        cycle_paths.spectral_period_seconds: metric_value(
            spectral.period_seconds,
            unit="s",
        ),
    }
    return outputs


def metric_value(
    data,
    *,
    unit: str | None = None,
    dim_desc: Iterable[str] | None = None,
    attrs: dict[str, object] | None = None,
):
    output_attrs: dict[str, object] = dict(attrs or {})
    if unit:
        output_attrs["unit"] = unit
    if dim_desc:
        output_attrs["dimDesc"] = list(dim_desc)
    value = metric_data(data)
    return (value, output_attrs) if output_attrs else value


def metric_data(data):
    if isinstance(data, bool):
        return data
    if isinstance(data, float):
        return np.float32(data)
    if isinstance(data, int):
        return np.int32(data)
    if isinstance(data, complex):
        return np.complex64(data)
    value = np.asarray(data)
    if value.dtype.kind == "f":
        return value.astype(np.float32, copy=False)
    if value.dtype.kind == "c":
        return value.astype(np.complex64, copy=False)
    if value.dtype.kind == "i":
        return value.astype(np.int32, copy=False)
    if value.dtype.kind == "u":
        return value.astype(np.uint32, copy=False)
    return value


def _resolve_output_paths(
    output_paths: EyeFlowOutputPaths | str | None,
) -> EyeFlowOutputPaths:
    if isinstance(output_paths, EyeFlowOutputPaths):
        return output_paths
    return EyeFlowOutputPaths.active(output_paths)


__all__ = ["metric_data", "metric_value", "pack_retinal_velocity_outputs"]
