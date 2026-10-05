"""Orchestrate absolute waveform metric products."""

from pipelines.waveform_velocity import waveform_velocity
from pipelines.waveform_velocity.outputs import (
    pack_velocity_per_beat_outputs,
)

from .outputs import pack_absolute_waveform_outputs


def run_absolute_waveform_metrics(ctx) -> dict[str, object]:
    """Calculate default global metrics and selected regional products."""
    selected = ctx.options_for("absolute_waveform_metrics")
    waveform = waveform_velocity(ctx)
    velocity_outputs = pack_velocity_per_beat_outputs(
        waveform.per_beat_result,
        velocity_analysis=waveform.retinal_velocity,
    )
    outputs = pack_absolute_waveform_outputs(
        velocity_outputs,
        source_data=waveform.source_data if "quadrants" in selected else None,
        artery_segments=(
            waveform.artery_segments if "quadrants" in selected else None
        ),
        vein_segments=(
            waveform.vein_segments if "quadrants" in selected else None
        ),
        include_per_beat=True,
        include_segments="segments" in selected,
        include_quadrants="quadrants" in selected,
    )
    ctx.state.set("absolute_waveform_metric_outputs", outputs)
    return outputs
