"""Orchestrate absolute waveform metric products."""

from pipelines.waveform_velocity_core.runner import (
    VELOCITY_PER_BEAT_OUTPUTS_STATE,
    WAVEFORM_CONTEXT_STATE,
)

from .outputs import pack_absolute_waveform_outputs


def run_absolute_waveform_metrics(ctx) -> dict[str, object]:
    """Calculate absolute metrics from shared per-beat velocity outputs."""
    selected = ctx.options_for("absolute_waveform_metrics")
    if not selected:
        return {}

    velocity_outputs = _required_state(ctx, VELOCITY_PER_BEAT_OUTPUTS_STATE)
    context = (
        _required_state(ctx, WAVEFORM_CONTEXT_STATE)
        if "hemifield" in selected
        else None
    )
    outputs = pack_absolute_waveform_outputs(
        velocity_outputs,
        source_data=context.source_data if context is not None else None,
        artery_segments=(
            context.artery_segment_result if context is not None else None
        ),
        vein_segments=(
            context.vein_segment_result if context is not None else None
        ),
        include_per_beat="per_beat" in selected,
        include_segments="segments" in selected,
        include_hemifield="hemifield" in selected,
    )
    ctx.state.set("absolute_waveform_metric_outputs", outputs)
    return outputs


def _required_state(ctx, key: str):
    value = ctx.state.get(key)
    if value is None:
        raise RuntimeError(
            f"Required pipeline state '{key}' is unavailable; "
            "check the pipeline DAG dependencies."
        )
    return value
