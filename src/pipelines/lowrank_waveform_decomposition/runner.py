"""Orchestrate low-rank waveform decomposition products."""

from pipelines.waveform_velocity_core.runner import (
    VELOCITY_PER_BEAT_OUTPUTS_STATE,
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
    outputs = pack_lowrank_waveform_decomposition_outputs(
        velocity_outputs,
        vein_flag="veins" in selected,
    )
    ctx.state.set(LOWRANK_WAVEFORM_OUTPUTS_STATE, outputs)
    return outputs


__all__ = [
    "LOWRANK_WAVEFORM_OUTPUTS_STATE",
    "run_lowrank_waveform_decomposition",
]
