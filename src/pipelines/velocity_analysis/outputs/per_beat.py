"""Output packing for per-beat velocity calculations."""

import numpy as np

from calculations.blood_flow_velocity import PerBeatAnalysisResult
from calculations.blood_flow_velocity.signal_analysis.per_beat.segments import (
    SEGMENT_PER_BEAT_DIM_DESC,
)
from input_output.schema import EyeFlowOutputPaths, VelocityPerBeatOutputPaths
from pipelines.velocity.models import RetinalVelocity
from pipelines.velocity.outputs import metric_data, metric_value
from pipelines.velocity.semantics import velocity_dataset_attrs

from .paths import resolve_output_paths


def pack_velocity_per_beat_outputs(
    result: PerBeatAnalysisResult,
    output_paths: EyeFlowOutputPaths | str | None = None,
    *,
    velocity: RetinalVelocity | None = None,
) -> dict[str, object]:
    schema = resolve_output_paths(output_paths)
    velocity_attrs = velocity_dataset_attrs(velocity)
    metrics = {
        schema.beat_period_seconds: metric_value(
            _matlab_row_vector(result.beat_period_seconds),
            unit="s",
            dim_desc=("row", "beat"),
        ),
    }
    metrics.update(
        _pack_vessel_outputs(
            schema.artery_per_beat,
            result.artery,
            velocity_attrs,
        )
    )
    metrics.update(
        _pack_vessel_outputs(
            schema.vein_per_beat,
            result.vein,
            velocity_attrs,
        )
    )
    metrics.update(
        _pack_safe_segment_outputs(
            schema.artery_per_beat_safe,
            result.artery.safe_segments,
            velocity_attrs,
        )
    )
    metrics.update(
        _pack_safe_segment_outputs(
            schema.vein_per_beat_safe,
            result.vein.safe_segments,
            velocity_attrs,
        )
    )
    return metrics


def _pack_vessel_outputs(
    paths: VelocityPerBeatOutputPaths,
    vessel,
    velocity_attrs: dict[str, object],
) -> dict[str, object]:
    signal = vessel.signal
    metrics = {
        paths.velocity_signal: metric_value(
            signal.velocity_signal_per_beat,
            attrs=velocity_attrs,
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
            attrs=velocity_attrs,
            dim_desc=("beat", "sample"),
        ),
    }
    if vessel.segments is not None:
        metrics.update(
            _pack_vessel_segment_outputs(paths, vessel.segments, velocity_attrs)
        )
    return metrics


def _pack_vessel_segment_outputs(
    paths: VelocityPerBeatOutputPaths,
    segments,
    velocity_attrs: dict[str, object],
) -> dict[str, object]:
    if paths.segment_velocity_signal is None:
        return {}
    if paths.segment_velocity_signal_band_limited is None:
        return {}
    return {
        paths.segment_velocity_signal: _segment_metric_value(
            segments.velocity_signal_per_beat_per_segment,
            attrs=velocity_attrs,
        ),
        paths.segment_velocity_signal_band_limited: _segment_metric_value(
            segments.velocity_signal_per_beat_per_segment_band_limited,
            attrs=velocity_attrs,
        ),
    }


def _pack_safe_segment_outputs(
    paths,
    segments,
    velocity_attrs: dict[str, object],
) -> dict[str, object]:
    if segments is None or paths.velocity_signal is None:
        return {}
    if paths.velocity_signal_band_limited is None:
        return {}
    return {
        paths.velocity_signal: _segment_metric_value(
            segments.velocity_signal_per_beat_per_segment,
            attrs=velocity_attrs,
        ),
        paths.velocity_signal_band_limited: _segment_metric_value(
            segments.velocity_signal_per_beat_per_segment_band_limited,
            attrs=velocity_attrs,
        ),
    }


def _segment_metric_value(data, *, attrs: dict[str, object]):
    value = metric_data(data)
    if value.ndim != 4:
        raise ValueError(
            "segment per-beat outputs must have shape "
            "(sample, beat, branch, radius)."
        )
    return metric_value(
        value,
        attrs=attrs,
        dim_desc=SEGMENT_PER_BEAT_DIM_DESC,
    )


def _matlab_row_vector(data) -> np.ndarray:
    return np.asarray(data).reshape(1, -1)


__all__ = ["pack_velocity_per_beat_outputs"]
