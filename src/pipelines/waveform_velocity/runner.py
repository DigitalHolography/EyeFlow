"""Orchestrate selectable velocity output products."""

from pipelines.waveform_velocity_core.runner import (
    VELOCITY_PER_BEAT_OUTPUTS_STATE,
    VELOCITY_PER_BEAT_RESULT_STATE,
    WAVEFORM_CONTEXT_STATE,
)
from pipelines.waveform_velocity_core.per_beat import run_velocity_per_beat_metrics

from .continuous import pack_continuous_velocity_outputs
from .hemifield import pack_hemifield_velocity_outputs
from .profiles import pack_cross_section_profile_outputs


def run_waveform_velocity(ctx) -> dict[str, object]:
    """Publish base velocity plus the selected derived velocity products."""
    context = _required_state(ctx, WAVEFORM_CONTEXT_STATE)
    selected = ctx.options_for("waveform_velocity")
    metrics = pack_continuous_velocity_outputs(context.dopplerview_analysis)

    per_beat_result = ctx.state.get(VELOCITY_PER_BEAT_RESULT_STATE)
    velocity_outputs = ctx.state.get(VELOCITY_PER_BEAT_OUTPUTS_STATE, {})
    if "per_beat" in selected:
        if per_beat_result is None:
            per_beat_result, velocity_outputs = run_velocity_per_beat_metrics(context)
            ctx.state.set(VELOCITY_PER_BEAT_RESULT_STATE, per_beat_result)
            ctx.state.set(VELOCITY_PER_BEAT_OUTPUTS_STATE, velocity_outputs)
        metrics.update(velocity_outputs)

    if "velocity_profiles" in selected:
        cycle_boundaries = (
            per_beat_result.cycle_boundary_indexes
            if per_beat_result is not None
            else context.per_beat_analysis.cycle_boundary_indexes
        )
        index_base = (
            0
            if per_beat_result is not None
            else int(context.source_data.provenance["beat_index_base"])
        )
        metrics.update(
            pack_cross_section_profile_outputs(
                context.artery_segment_result,
                context.vein_segment_result,
                cycle_boundaries,
                index_base=index_base,
            )
        )

    if "hemifield" in selected:
        metrics.update(
            pack_hemifield_velocity_outputs(
                velocity_outputs,
                context.source_data,
                context.artery_segment_result,
                context.vein_segment_result,
            )
        )

    return metrics


def _required_state(ctx, key: str):
    value = ctx.state.get(key)
    if value is None:
        raise RuntimeError(
            f"Required pipeline state '{key}' is unavailable; "
            "check the pipeline DAG dependencies."
        )
    return value
