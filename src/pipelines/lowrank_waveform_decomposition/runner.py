"""Orchestrate low-rank waveform decomposition products."""

from pipelines.waveform_velocity import waveform_velocity
from pipelines.waveform_velocity.outputs import (
    pack_velocity_per_beat_outputs,
)

from .outputs import pack_lowrank_waveform_decomposition_outputs

LOWRANK_WAVEFORM_OUTPUTS_STATE = "lowrank_waveform_decomposition_outputs"


def run_lowrank_waveform_decomposition(ctx) -> dict[str, object]:
    """Calculate joint and per-beat low-rank products from segment waveforms."""
    waveform = waveform_velocity(ctx)
    velocity_outputs = pack_velocity_per_beat_outputs(
        waveform.per_beat_result,
        velocity_analysis=waveform.retinal_velocity,
    )
    selected = ctx.options_for("lowrank_waveform_decomposition")
    include_quadrants = "quadrants" in selected
    outputs = pack_lowrank_waveform_decomposition_outputs(
        velocity_outputs,
        vein_flag=True,
        include_quadrants=include_quadrants,
        artery_segments=(
            waveform.artery_segments if include_quadrants else None
        ),
        vein_segments=waveform.vein_segments if include_quadrants else None,
    )
    ctx.state.set(LOWRANK_WAVEFORM_OUTPUTS_STATE, outputs)
    return outputs


__all__ = [
    "LOWRANK_WAVEFORM_OUTPUTS_STATE",
    "run_lowrank_waveform_decomposition",
]
