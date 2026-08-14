"""Orchestrate low-rank waveform decomposition products."""

from pipelines.waveform_velocity_core.runner import (
    VELOCITY_PER_BEAT_OUTPUTS_STATE,
    WAVEFORM_CONTEXT_STATE,
)

from .outputs import pack_lowrank_waveform_decomposition_outputs

LOWRANK_WAVEFORM_OUTPUTS_STATE = "lowrank_waveform_decomposition_outputs"


def run_lowrank_waveform_decomposition(ctx) -> dict[str, object]:
    """Calculate joint and per-beat low-rank products from segment waveforms."""
    velocity_outputs = ctx.state.get(VELOCITY_PER_BEAT_OUTPUTS_STATE)
    if velocity_outputs is None:
        raise RuntimeError(
            f"Required pipeline state '{VELOCITY_PER_BEAT_OUTPUTS_STATE}' is "
            "unavailable; check the pipeline DAG dependencies."
        )
    selected = ctx.options_for("lowrank_waveform_decomposition")
    include_quadrants = "quadrants" in selected
    context = ctx.state.get(WAVEFORM_CONTEXT_STATE) if include_quadrants else None
    if include_quadrants and context is None:
        raise RuntimeError(
            f"Required pipeline state '{WAVEFORM_CONTEXT_STATE}' is unavailable; "
            "quadrant low-rank outputs require segment geometry."
        )
    outputs = pack_lowrank_waveform_decomposition_outputs(
        velocity_outputs,
        vein_flag=True,
        include_quadrants=include_quadrants,
        source_data=context.source_data if context is not None else None,
        artery_segments=(
            context.artery_segment_result if context is not None else None
        ),
        vein_segments=context.vein_segment_result if context is not None else None,
    )
    ctx.state.set(LOWRANK_WAVEFORM_OUTPUTS_STATE, outputs)
    return outputs


__all__ = [
    "LOWRANK_WAVEFORM_OUTPUTS_STATE",
    "run_lowrank_waveform_decomposition",
]
