"""Compose the public waveform-shape metric output groups."""

from __future__ import annotations

from collections.abc import Mapping

from input_output import EyeFlowOutputPaths

from .metrics.runner import run_waveform_shape_metric_calculations
from .metrics.hemifield import pack_hemifield_metrics


def pack_waveform_shape_outputs(
    metrics: Mapping[str, object],
    source_data,
    artery_segments,
    vein_segments,
    output_paths: EyeFlowOutputPaths | str | None = None,
    *,
    include_per_beat: bool = True,
    include_hemifield: bool = True,
) -> dict[str, object]:
    """Build only the selected waveform-shape metric output groups.

    The calculator-relative global and by-segment groups are first added below
    the configured waveform-shape root.  The same calculated segment metrics
    are then reused for the region-indexed hemifield groups.
    """
    if not include_per_beat and not include_hemifield:
        return {}

    global_and_by_segment_outputs = run_waveform_shape_metric_calculations(
        metrics,
        output_paths,
    )
    metrics_with_waveform_outputs = dict(metrics)
    metrics_with_waveform_outputs.update(global_and_by_segment_outputs)
    hemifield_outputs = pack_hemifield_metrics(
        metrics_with_waveform_outputs,
        source_data,
        artery_segments,
        vein_segments,
        output_paths,
    ) if include_hemifield else {}
    return {
        **(global_and_by_segment_outputs if include_per_beat else {}),
        **hemifield_outputs,
    }
