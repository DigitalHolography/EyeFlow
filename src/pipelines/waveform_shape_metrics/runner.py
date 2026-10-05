"""Orchestrate selectable waveform-shape metric products."""

from pipelines.waveform_velocity import waveform_velocity
from pipelines.waveform_velocity.outputs import (
    pack_velocity_per_beat_outputs,
)

from .outputs import pack_waveform_shape_outputs


def run_waveform_shape_metrics(ctx) -> dict[str, object]:
    """Calculate only the selected waveform-shape metric products."""
    selected = ctx.options_for("waveform_shape_metrics")
    report_required = ctx.pipeline_scheduled("pdf_report")
    if not selected and not report_required:
        return {}

    waveform = waveform_velocity(ctx)
    velocity_outputs = pack_velocity_per_beat_outputs(
        waveform.require_per_beat()
    )
    outputs = pack_waveform_shape_outputs(
        velocity_outputs,
        waveform.source_data,
        waveform.artery_segments,
        waveform.vein_segments,
        include_per_beat="per_beat" in selected or report_required,
        include_segments="segments" in selected,
        include_quadrants="quadrants" in selected,
    )
    ctx.state.set("waveform_shape_metric_outputs", outputs)
    return outputs
