"""Orchestrate low-rank waveform decomposition products."""

from pipelines.velocity_analysis import velocity_analysis
from pipelines.velocity_analysis.outputs import (
    pack_velocity_per_beat_inputs,
)

from .outputs import pack_lowrank_waveform_decomposition_outputs

LOWRANK_WAVEFORM_OUTPUTS_STATE = "lowrank_waveform_decomposition_outputs"


def run_lowrank_waveform_decomposition(ctx) -> dict[str, object]:
    """Calculate joint and per-beat low-rank products from segment waveforms."""
    analysis = velocity_analysis(ctx)
    velocity_outputs = pack_velocity_per_beat_inputs(
        analysis.per_beat_result,
        velocity=analysis.velocity,
    )
    selected = ctx.options_for("lowrank_waveform_decomposition")
    include_quadrants = "quadrants" in selected
    outputs = pack_lowrank_waveform_decomposition_outputs(
        velocity_outputs,
        vein_flag=True,
        include_quadrants=include_quadrants,
        artery_segments=(
            analysis.artery_segments if include_quadrants else None
        ),
        vein_segments=analysis.vein_segments if include_quadrants else None,
    )
    ctx.state.set(LOWRANK_WAVEFORM_OUTPUTS_STATE, outputs)
    return outputs


__all__ = [
    "LOWRANK_WAVEFORM_OUTPUTS_STATE",
    "run_lowrank_waveform_decomposition",
]
