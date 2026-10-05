"""Orchestrate selectable waveform-shape metric products."""

from pipelines.waveform_velocity import waveform_velocity
from pipelines.waveform_velocity.outputs import (
    pack_velocity_per_beat_outputs,
)

from .outputs import pack_waveform_shape_outputs


def run_waveform_shape_metrics(ctx) -> dict[str, object]:
    """Calculate default global metrics and selected regional products."""
    selected = ctx.options_for("waveform_shape_metrics")
    waveform = waveform_velocity(ctx)
    velocity_outputs = pack_velocity_per_beat_outputs(
        waveform.per_beat_result,
        velocity_analysis=waveform.retinal_velocity,
    )
    outputs = pack_waveform_shape_outputs(
        velocity_outputs,
        waveform.source_data,
        waveform.artery_segments,
        waveform.vein_segments,
        include_per_beat=True,
        include_segments="segments" in selected,
        include_quadrants="quadrants" in selected,
    )
    ctx.state.set("waveform_shape_metric_outputs", outputs)
    return outputs
