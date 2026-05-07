"""Domain-specific metric dictionaries for EyeFlow HDF5 output paths."""

from __future__ import annotations

from collections.abc import Iterable, Mapping

import numpy as np

from calculations.blood_flow_velocity import PerBeatAnalysisResult

from .schema import (
    DOPPLER_VIEW_ANALYSIS_SCHEMA,
    EyeFlowOutputSchemaVariant,
    VelocityPerBeatOutputPaths,
    output_schema_variant,
)

ZERO_BASED_INDEX_PATHS = frozenset(
    {
        DOPPLER_VIEW_ANALYSIS_SCHEMA.dataset_path("beat_indices"),
    }
)


def systolic_index_base_for_path(path: str) -> int | None:
    from .hdf5 import normalize_h5_path

    normalized = normalize_h5_path(path)
    return 0 if normalized in ZERO_BASED_INDEX_PATHS else None


def pack_velocity_per_beat_outputs(
    result: PerBeatAnalysisResult,
    schema_variant: EyeFlowOutputSchemaVariant | str | None = None,
) -> dict[str, object]:
    schema = _resolve_output_schema(schema_variant)
    metrics = {
        schema.beat_period_idx: _metric_value(
            result.beat_period_idx,
            unit="frame",
            dim_desc=("beat",),
        ),
        schema.beat_period_seconds: _metric_value(
            result.beat_period_seconds,
            unit="s",
            dim_desc=("beat",),
        ),
    }
    metrics.update(_pack_vessel_outputs(schema.artery_per_beat, result.artery))
    metrics.update(_pack_vessel_outputs(schema.vein_per_beat, result.vein))
    return metrics


def pack_dopplerview_analysis_outputs(
    analysis: Mapping[str, object],
    schema_variant: EyeFlowOutputSchemaVariant | str | None = None,
) -> dict[str, object]:
    paths = _resolve_output_schema(schema_variant).analysis
    metrics = {
        paths.retinal_artery_velocity_signal: _metric_value(
            analysis["retinal_artery_velocity_signal"],
            unit="mm/s",
        ),
        paths.retinal_vein_velocity_signal: _metric_value(
            analysis["retinal_vein_velocity_signal"],
            unit="mm/s",
        ),
        paths.velocity_map_avg: _metric_value(analysis["velocity_map_avg"]),
        paths.fRMS_avg: _metric_value(analysis["fRMS_avg"]),
        paths.fRMS_bkg_avg: _metric_value(analysis["fRMS_bkg_avg"]),
        paths.velocitysignal_per_beat: _metric_value(
            analysis["retinal_artery_velocity_signal_filtered_perbeat"],
            unit="mm/s",
        ),
        paths.velocitysignal_filtered: _metric_value(
            analysis["retinal_artery_velocity_signal_filtered"],
            unit="mm/s",
        ),
        paths.beat_indices: _metric_value(analysis["beat_indices"]),
        paths.time_per_beat: _metric_value(
            analysis["time_per_beat"],
            unit="s",
        ),
    }
    if paths.retinal_velocity_array is not None:
        metrics[paths.retinal_velocity_array] = _metric_value(
            analysis["retinal_vessel_velocity"],
            unit="mm/s",
        )
    return metrics


def _resolve_output_schema(
    schema_variant: EyeFlowOutputSchemaVariant | str | None,
) -> EyeFlowOutputSchemaVariant:
    if isinstance(schema_variant, EyeFlowOutputSchemaVariant):
        return schema_variant
    return output_schema_variant(schema_variant)


def _pack_vessel_outputs(
    paths: VelocityPerBeatOutputPaths,
    vessel,
) -> dict[str, object]:
    signal = vessel.signal
    metrics = {
        paths.velocity_signal: _metric_value(
            signal.velocity_signal_per_beat,
            unit="mm/s",
            dim_desc=("beat", "sample"),
        ),
        paths.velocity_signal_fft_abs: _metric_value(
            np.abs(signal.velocity_signal_per_beat_fft),
            unit="a.u.",
            dim_desc=("beat", "frequency_bin"),
        ),
        paths.velocity_signal_fft_arg: _metric_value(
            np.angle(signal.velocity_signal_per_beat_fft),
            unit="rad",
            dim_desc=("beat", "frequency_bin"),
        ),
        paths.velocity_signal_band_limited: _metric_value(
            signal.velocity_signal_per_beat_band_limited,
            unit="mm/s",
            dim_desc=("beat", "sample"),
        ),
    }
    if vessel.segments is not None:
        metrics.update(_pack_vessel_segment_outputs(paths, vessel.segments))
    return metrics


def _pack_vessel_segment_outputs(
    paths: VelocityPerBeatOutputPaths,
    segments,
) -> dict[str, object]:
    if paths.segment_velocity_signal is None:
        return {}
    if paths.segment_velocity_signal_band_limited is None:
        return {}
    return {
        paths.segment_velocity_signal: _segment_metric_value(
            segments.velocity_signal_per_beat_per_segment,
            unit="mm/s",
        ),
        paths.segment_velocity_signal_band_limited: _segment_metric_value(
            segments.velocity_signal_per_beat_per_segment_band_limited,
            unit="mm/s",
        ),
    }


def _segment_metric_value(data, *, unit: str):
    value, dim_desc = _drop_singleton_radius_axis(data)
    return _metric_value(value, unit=unit, dim_desc=dim_desc)


def _drop_singleton_radius_axis(data) -> tuple[np.ndarray, tuple[str, ...]]:
    value = _metric_data(data)
    if value.ndim == 4 and value.shape[-1] == 1:
        return value[..., 0], ("sample", "beat", "branch")
    return value, ("sample", "beat", "branch", "radius")


def _metric_value(
    data,
    *,
    unit: str | None = None,
    dim_desc: Iterable[str] | None = None,
):
    attrs: dict[str, object] = {}
    if unit:
        attrs["unit"] = unit
    if dim_desc:
        attrs["dimDesc"] = list(dim_desc)
    data = _metric_data(data)
    return (data, attrs) if attrs else data


def _metric_data(data):
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
