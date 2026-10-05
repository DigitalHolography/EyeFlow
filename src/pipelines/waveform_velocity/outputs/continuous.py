"""Whole-vessel continuous velocity output packing."""

import numpy as np

from input_output.schema import EyeFlowOutputPaths
from pipelines.retinal_velocity.models import RetinalVelocity
from pipelines.retinal_velocity.outputs import metric_value
from pipelines.retinal_velocity.semantics import velocity_dataset_attrs

from ..analysis.filtering import lowpass_velocity_signals
from .paths import resolve_output_paths


def pack_continuous_velocity_outputs(
    velocity: RetinalVelocity,
    output_paths: EyeFlowOutputPaths | str | None = None,
) -> dict[str, object]:
    """Pack raw and band-limited artery and vein velocity signals."""
    schema = resolve_output_paths(output_paths)
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

    schema = resolve_output_paths(output_paths)
    attrs = velocity_dataset_attrs(velocity_analysis)
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
            "segment velocity must have shape (radius, branch, frame), "
            f"got {values.shape}."
        )

    outputs = {
        paths.velocity_signal: metric_value(
            values.transpose(2, 1, 0),
            attrs=attrs,
            dim_desc=("frame", "branch", "radius"),
        )
    }
    if paths.velocity_signal_band_limited is not None:
        band_limited = _lowpass_segment_velocity(values, source_data)
        outputs[paths.velocity_signal_band_limited] = metric_value(
            band_limited.transpose(2, 1, 0),
            attrs=attrs,
            dim_desc=("frame", "branch", "radius"),
        )
    return outputs


def _lowpass_segment_velocity(values: np.ndarray, source_data) -> np.ndarray:
    timing = source_data.source.holodoppler.timing
    if timing is None:
        raise ValueError("Segment band-limited velocity requires source timing.")

    return lowpass_velocity_signals(
        values,
        dt_seconds=float(timing.dt_seconds),
    )
