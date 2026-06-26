"""Output packing for per-beat velocity calculations."""

from collections.abc import Iterable

import numpy as np

from calculations.blood_flow_velocity import PerBeatAnalysisResult
from input_output.schema import EyeFlowOutputPaths, VelocityPerBeatOutputPaths


def pack_velocity_per_beat_outputs(
    result: PerBeatAnalysisResult,
    output_paths: EyeFlowOutputPaths | str | None = None,
) -> dict[str, object]:
    schema = _resolve_output_paths(output_paths)
    metrics = {
        schema.beat_period_idx: metric_value(
            _matlab_row_vector(result.beat_period_idx),
            unit="frame",
            dim_desc=("row", "beat"),
        ),
        schema.beat_period_seconds: metric_value(
            _matlab_row_vector(result.beat_period_seconds),
            unit="s",
            dim_desc=("row", "beat"),
        ),
    }
    metrics.update(_pack_vessel_outputs(schema.artery_per_beat, result.artery))
    metrics.update(_pack_vessel_outputs(schema.vein_per_beat, result.vein))
    return metrics


def _pack_vessel_outputs(
    paths: VelocityPerBeatOutputPaths,
    vessel,
) -> dict[str, object]:
    signal = vessel.signal
    metrics = {
        paths.velocity_signal: metric_value(
            signal.velocity_signal_per_beat,
            unit="mm/s",
            dim_desc=("beat", "sample"),
        ),
        paths.velocity_signal_fft_abs: metric_value(
            np.abs(signal.velocity_signal_per_beat_fft),
            unit="a.u.",
            dim_desc=("beat", "frequency_bin"),
        ),
        paths.velocity_signal_fft_arg: metric_value(
            np.angle(signal.velocity_signal_per_beat_fft),
            unit="rad",
            dim_desc=("beat", "frequency_bin"),
        ),
        paths.velocity_signal_band_limited: metric_value(
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
    value = metric_data(data)
    if value.ndim != 4:
        raise ValueError(
            "segment per-beat outputs must have shape "
            "(sample, beat, branch, radius)."
        )
    return metric_value(
        value,
        unit=unit,
        dim_desc=("sample", "beat", "branch", "radius"),
    )


def _matlab_row_vector(data) -> np.ndarray:
    return np.asarray(data).reshape(1, -1)


def _resolve_output_paths(
    output_paths: EyeFlowOutputPaths | str | None,
) -> EyeFlowOutputPaths:
    if isinstance(output_paths, EyeFlowOutputPaths):
        return output_paths
    return EyeFlowOutputPaths.active(output_paths)


def metric_value(
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
    data = metric_data(data)
    return (data, attrs) if attrs else data


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
