"""Orchestrate selectable waveform-shape metric products."""

from pipelines.waveform_velocity_core.runner import (
    VELOCITY_PER_BEAT_OUTPUTS_STATE,
    WAVEFORM_CONTEXT_STATE,
)

from .outputs import pack_waveform_shape_outputs


def run_waveform_shape_metrics(ctx) -> dict[str, object]:
    """Calculate only the selected waveform-shape metric products."""
    selected = ctx.options_for("waveform_shape_metrics")
    if not selected:
        return {}

    context = _required_state(ctx, WAVEFORM_CONTEXT_STATE)
    velocity_outputs = _required_state(ctx, VELOCITY_PER_BEAT_OUTPUTS_STATE)
    outputs = pack_waveform_shape_outputs(
        velocity_outputs,
        context.source_data,
        context.artery_segment_result,
        context.vein_segment_result,
        include_per_beat="per_beat" in selected,
        include_hemifield="hemifield" in selected,
    )
    ctx.state.set("waveform_shape_metric_outputs", outputs)
    return outputs


def _required_state(ctx, key: str):
    value = ctx.state.get(key)
    if value is None:
        raise RuntimeError(
            f"Required pipeline state '{key}' is unavailable; "
            "check the pipeline DAG dependencies."
        )
    return value
