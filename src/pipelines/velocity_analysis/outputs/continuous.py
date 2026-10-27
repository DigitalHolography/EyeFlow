"""Whole-vessel continuous velocity output packing."""

import numpy as np

from input_output.schema import EyeFlowOutputPaths
from input_output.schema.serialization import serialize_segment_velocity
from pipelines.velocity.models import RetinalVelocity
from pipelines.velocity.outputs import metric_value
from pipelines.velocity.semantics import velocity_dataset_attrs

from ..analysis.filtering import lowpass_velocity_signals
from .paths import resolve_output_paths


def pack_continuous_velocity_outputs(
    velocity: RetinalVelocity,
    output_paths: EyeFlowOutputPaths | str | None = None,
) -> dict[str, object]:
    """Pack raw and filtered artery and vein velocity signals."""
    schema = resolve_output_paths(output_paths)
    paths = schema.analysis
    attrs = velocity_dataset_attrs(velocity)
    filtered_attrs = {**attrs, **continuous_filter_attrs(velocity.dt_seconds)}
    artery_raw = velocity.continuous("artery", raw=True)
    vein_raw = velocity.continuous("vein", raw=True)
    artery_filtered = velocity.continuous("artery")
    vein_filtered = velocity.continuous("vein")
    return {
        paths.retinal_artery_velocity_signal: metric_value(
            artery_raw,
            attrs=attrs,
        ),
        paths.retinal_vein_velocity_signal: metric_value(
            vein_raw,
            attrs=attrs,
        ),
        paths.retinal_artery_velocity_signal_filtered: metric_value(
            artery_filtered,
            attrs=filtered_attrs,
        ),
        paths.retinal_vein_velocity_signal_filtered: metric_value(
            vein_filtered,
            attrs=filtered_attrs,
        ),
    }


def pack_segment_velocity_outputs(
    artery_segments,
    vein_segments,
    output_paths: EyeFlowOutputPaths | str | None = None,
    *,
    source_data=None,
    velocity: RetinalVelocity | None = None,
) -> dict[str, object]:
    """Pack continuous segment velocity signals without beat decomposition."""

    schema = resolve_output_paths(output_paths)
    attrs = velocity_dataset_attrs(velocity)
    return {
        **_pack_segment_velocity_output(
            artery_segments,
            schema.artery_segments,
            source_data,
            attrs,
        ),
        **_pack_segment_velocity_output(
            vein_segments,
            schema.vein_segments,
            source_data,
            attrs,
        ),
    }


def _pack_segment_velocity_output(
    segments,
    paths,
    source_data,
    attrs: dict[str, object],
) -> dict[str, object]:
    if segments is None or paths.velocity_signal is None:
        return {}
    profile = segments.profile
    if np.asarray(profile.topology.native.branch_ids).size == 0:
        return {}

    values = np.asarray(profile.segment_signal, dtype=np.float32)
    if values.ndim != 3:
        raise ValueError(
            f"segment velocity must have shape (radius, branch, frame), got {values.shape}."
        )

    outputs = {
        paths.velocity_signal: metric_value(
            serialize_segment_velocity(values),
            attrs=attrs,
            dim_desc=("frame", "branch", "radius"),
        )
    }
    if paths.velocity_signal_filtered is not None:
        filtered = getattr(segments, "filtered_signal", None)
        if filtered is None:
            filtered = _lowpass_segment_velocity(values, source_data)
        outputs[paths.velocity_signal_filtered] = metric_value(
            serialize_segment_velocity(filtered),
            attrs={
                **attrs,
                **continuous_filter_attrs(float(source_data.source.holodoppler.timing.dt_seconds)),
            },
            dim_desc=("frame", "branch", "radius"),
        )
    return outputs


def _lowpass_segment_velocity(values: np.ndarray, source_data) -> np.ndarray:
    timing = source_data.source.holodoppler.timing
    if timing is None:
        raise ValueError("Segment filtered velocity requires source timing.")

    return lowpass_velocity_signals(
        values,
        dt_seconds=float(timing.dt_seconds),
    )


def continuous_filter_attrs(dt_seconds: float) -> dict[str, object]:
    """Describe the continuous Butterworth stage independently of bandlimiting."""
    return {
        "signal_source": "filtered",
        "filter_method": "butterworth_filtfilt",
        "filter_order": 4,
        "filter_cutoff_hz": 15.0,
        "filter_padding": "odd_reflection",
        "filter_padding_samples": 12,
        "filter_application": "continuous_before_beat_subdivision",
        "filter_short_series_policy": "leave_unchanged_when_samples_le_12",
        "dt_seconds": float(dt_seconds),
    }
