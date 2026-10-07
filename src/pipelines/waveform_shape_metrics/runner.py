"""Orchestrate selectable waveform-shape metric products."""

from pipelines.velocity_analysis import velocity_analysis
from pipelines.velocity_analysis.outputs import (
    pack_velocity_per_beat_inputs,
)

from .outputs import pack_waveform_shape_outputs


def run_waveform_shape_metrics(ctx) -> dict[str, object]:
    """Calculate default global metrics and selected regional products."""
    selected = ctx.options_for("waveform_shape_metrics")
    analysis = velocity_analysis(ctx)
    velocity_outputs = pack_velocity_per_beat_inputs(
        analysis.per_beat_result,
        velocity=analysis.velocity,
    )
    outputs = pack_waveform_shape_outputs(
        velocity_outputs,
        analysis.source_data,
        analysis.artery_segments,
        analysis.vein_segments,
        include_per_beat=True,
        include_segments="segments" in selected,
        include_quadrants="quadrants" in selected,
    )
    ctx.state.set("waveform_shape_metric_outputs", outputs)
    return outputs
