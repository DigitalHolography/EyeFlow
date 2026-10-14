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
    metrics = {}
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
    for paths, vessel in (
        (schema.artery_per_beat_filtered, result.artery_filtered),
        (schema.vein_per_beat_filtered, result.vein_filtered),
    ):
        if vessel is not None:
            metrics.update(_pack_vessel_outputs(paths, vessel, velocity_attrs))
    for path, payload in metrics.items():
        attrs = payload[1]
        source_path = path
        for representation in ("BandLimited", "FFTAbs", "FFTPhase"):
            source_path = source_path.replace(f"/{representation}/", "/")
        if source_path == path:
            source_path = path.replace("/PerBeat/", "/Continuous/")
        attrs["source_path"] = "/" + source_path
        attrs["beat_interpolation"] = "interpft_Nfft_plus_one_drop_endpoint"
        attrs["phase_endpoint"] = "excluded"
        attrs["signal_source"] = "filtered" if "/Filtered/" in path else "raw"
        if "/Filtered/" in path:
            attrs.update(
                filter_method="butterworth_filtfilt",
                filter_order=4,
                filter_cutoff_hz=15.0,
                filter_application="continuous_before_beat_subdivision",
            )
        if "/BandLimited/" in path:
            sample_count = np.asarray(payload[0]).shape[0 if "/Segments/" in path else 1]
            attrs.update(
                band_limit_harmonics=result.band_limited_harmonics,
                band_limit_method="fourier_harmonic_truncation",
                band_limit_highest_retained_harmonic=min(
                    result.band_limited_harmonics, sample_count // 2
                ),
                band_limit_includes_DC=True,
                band_limit_reconstruction="matlab_magnitude",
                band_limit_boundary="periodic_within_each_beat",
            )
    return metrics


def pack_velocity_per_beat_inputs(
    result: PerBeatAnalysisResult,
    output_paths: EyeFlowOutputPaths | str | None = None,
    *,
    velocity: RetinalVelocity,
) -> dict[str, object]:
    """Pack per-beat waveforms plus canonical cardiac-cycle durations."""

    schema = resolve_output_paths(output_paths)
    inputs = pack_velocity_per_beat_outputs(
        result,
        schema,
        velocity=velocity,
    )
    inputs[schema.cardiac_cycle.systolic_cycle_duration_seconds] = metric_value(
        velocity.cycle_durations_seconds,
        unit="s",
    )
    return inputs


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
        metrics.update(_pack_vessel_segment_outputs(paths, vessel.segments, velocity_attrs))
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


def _segment_metric_value(data, *, attrs: dict[str, object]):
    value = metric_data(data)
    if value.ndim != 4:
        raise ValueError("segment per-beat outputs must have shape (sample, beat, branch, radius).")
    return metric_value(
        value,
        attrs=attrs,
        dim_desc=SEGMENT_PER_BEAT_DIM_DESC,
    )


__all__ = [
    "pack_velocity_per_beat_inputs",
    "pack_velocity_per_beat_outputs",
]
