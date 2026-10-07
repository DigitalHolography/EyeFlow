"""Orchestrate absolute waveform metric products."""

from pipelines.velocity_analysis import velocity_analysis
from pipelines.velocity_analysis.outputs import (
    pack_velocity_per_beat_inputs,
)

from .outputs import pack_absolute_waveform_outputs


def run_absolute_waveform_metrics(ctx) -> dict[str, object]:
    """Calculate default global metrics and selected regional products."""
    selected = ctx.options_for("absolute_waveform_metrics")
    analysis = velocity_analysis(ctx)
    velocity_outputs = pack_velocity_per_beat_inputs(
        analysis.per_beat_result,
        velocity=analysis.velocity,
    )
    outputs = pack_absolute_waveform_outputs(
        velocity_outputs,
        source_data=analysis.source_data if "quadrants" in selected else None,
        artery_segments=(
            analysis.artery_segments if "quadrants" in selected else None
        ),
        vein_segments=(
            analysis.vein_segments if "quadrants" in selected else None
        ),
        include_per_beat=True,
        include_segments="segments" in selected,
        include_quadrants="quadrants" in selected,
    )
    ctx.state.set("absolute_waveform_metric_outputs", outputs)
    return outputs
