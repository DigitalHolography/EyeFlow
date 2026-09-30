"""Whole-vessel continuous velocity output packing."""

import numpy as np

from calculations.math import butter_lowpass_filtfilt
from input_output.schema import EyeFlowOutputPaths
from pipelines.retinal_velocity.models import RetinalVelocity
from pipelines.retinal_velocity.outputs import metric_value
from pipelines.retinal_velocity.semantics import (
    resolve_velocity_semantics,
    velocity_dataset_attrs,
)
from pipelines.retinal_velocity.signal_processing import (
    DEFAULT_VELOCITY_SIGNAL_LOWPASS_HZ,
)


def pack_continuous_velocity_outputs(
    velocity: RetinalVelocity,
    output_paths: EyeFlowOutputPaths | str | None = None,
) -> dict[str, object]:
    """Pack raw and band-limited artery and vein velocity signals."""
    schema = _resolve_output_paths(output_paths)
    paths = schema.analysis
    attrs = velocity_dataset_attrs(velocity)
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
        paths.retinal_artery_velocity_signal_band_limited: metric_value(
            artery_filtered,
            attrs=attrs,
        ),
        paths.retinal_vein_velocity_signal_band_limited: metric_value(
            vein_filtered,
            attrs=attrs,
        ),
    }
def pack_segment_velocity_outputs(
    artery_segments,
    vein_segments,
    output_paths: EyeFlowOutputPaths | str | None = None,
    *,
    source_data=None,
    velocity_analysis: RetinalVelocity | None = None,
) -> dict[str, object]:
    """Pack continuous segment velocity signals without beat decomposition."""

    schema = _resolve_output_paths(output_paths)
    unit = resolve_velocity_semantics(velocity_analysis).unit
    return {
        **_pack_segment_velocity_output(
            artery_segments,
            schema.artery_segments,
            source_data,
            unit,
        ),
        **_pack_segment_velocity_output(
            vein_segments,
            schema.vein_segments,
            source_data,
            unit,
        ),
    }


def _pack_segment_velocity_output(
    segments,
    paths,
    source_data,
    unit: str,
) -> dict[str, object]:
    if segments is None or paths.velocity_signal is None:
        return {}
    if np.asarray(segments.branch_ids).size == 0:
        return {}

    values = np.asarray(segments.projected_signal, dtype=np.float32)
    if values.ndim != 3:
        raise ValueError(
            "segment velocity must have shape (radius, branch, frame), "
            f"got {values.shape}."
        )

    outputs = {
        paths.velocity_signal: metric_value(
            values.transpose(2, 1, 0),
            unit=unit,
            dim_desc=("frame", "branch", "radius"),
        )
    }
    if paths.velocity_signal_band_limited is not None:
        band_limited = _lowpass_segment_velocity(values, source_data)
        outputs[paths.velocity_signal_band_limited] = metric_value(
            band_limited.transpose(2, 1, 0),
            unit=unit,
            dim_desc=("frame", "branch", "radius"),
        )
    return outputs


def _lowpass_segment_velocity(values: np.ndarray, source_data) -> np.ndarray:
    timing = source_data.source.holodoppler.timing
    if timing is None:
        raise ValueError("Segment band-limited velocity requires source timing.")

    filtered = np.full(values.shape, np.nan, dtype=np.float32)
    for radius_index in range(values.shape[0]):
        for branch_index in range(values.shape[1]):
            signal = values[radius_index, branch_index]
            if not np.any(np.isfinite(signal)):
                continue
            filtered[radius_index, branch_index] = butter_lowpass_filtfilt(
                signal,
                dt_seconds=np.float32(timing.dt_seconds),
                lowpass_freq_hz=np.float32(DEFAULT_VELOCITY_SIGNAL_LOWPASS_HZ),
                order=4,
            )
    return filtered


def _resolve_output_paths(
    output_paths: EyeFlowOutputPaths | str | None,
) -> EyeFlowOutputPaths:
    if isinstance(output_paths, EyeFlowOutputPaths):
        return output_paths
    return EyeFlowOutputPaths.active(output_paths)
