"""Compose the public waveform-shape metric output groups."""

from __future__ import annotations

from collections.abc import Mapping

from input_output import EyeFlowOutputPaths

from .metrics.runner import run_waveform_shape_metric_calculations
from .metrics.quadrants import pack_quadrant_metrics


def pack_waveform_shape_outputs(
    metrics: Mapping[str, object],
    source_data,
    artery_segments,
    vein_segments,
    output_paths: EyeFlowOutputPaths | str | None = None,
    *,
    include_per_beat: bool = True,
    include_segments: bool = True,
    include_quadrants: bool = True,
) -> dict[str, object]:
    """Build only the selected waveform-shape metric output groups.

    The calculator-relative global and by-segment groups are first added below
    the configured waveform-shape root.  The same calculated segment metrics
    are then reused for the quadrant-indexed groups.
    """
    if not include_per_beat and not include_quadrants:
        return {}

    global_and_by_segment_outputs = run_waveform_shape_metric_calculations(
        metrics,
        output_paths,
        include_segments=include_segments or include_quadrants,
    )
    metrics_with_waveform_outputs = dict(metrics)
    metrics_with_waveform_outputs.update(global_and_by_segment_outputs)
    quadrant_outputs = pack_quadrant_metrics(
        metrics_with_waveform_outputs,
        source_data,
        artery_segments,
        vein_segments,
        output_paths,
    ) if include_quadrants else {}
    return {
        **{
            key: value
            for key, value in global_and_by_segment_outputs.items()
            if include_per_beat
            and (include_segments or "/by_segment/" not in key)
        },
        **quadrant_outputs,
    }
