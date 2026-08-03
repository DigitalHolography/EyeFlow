"""Compose the public waveform-shape metric output groups."""

from __future__ import annotations

from collections.abc import Mapping

from input_output import EyeFlowOutputPaths

from .metrics.runner import run_waveform_shape_metric_calculations
from .velocity.hemifield import (
    pack_hemifield_metrics,
    pack_hemifield_velocity_outputs,
)
from .velocity.segmentation import pack_segmentation_outputs


def pack_waveform_shape_outputs(
    metrics: Mapping[str, object],
    source_data,
    artery_segments,
    vein_segments,
    output_paths: EyeFlowOutputPaths | str | None = None,
) -> dict[str, object]:
    """Build global, by-segment, and hemifield waveform-shape outputs.

    The calculator-relative global and by-segment groups are first added below
    the configured waveform-shape root.  The same calculated segment metrics
    are then reused for the region-indexed hemifield groups.  The spatial
    segmentation outputs are composed here as well so the pipeline has one
    public waveform-output boundary.
    """
    segmentation_outputs = pack_segmentation_outputs(
        source_data,
        artery_segments,
        vein_segments,
        output_paths,
    )
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
    )
    hemifield_velocity_outputs = pack_hemifield_velocity_outputs(
        metrics,
        source_data,
        artery_segments,
        vein_segments,
        output_paths,
    )
    return {
        **segmentation_outputs,
        **global_and_by_segment_outputs,
        **hemifield_outputs,
        **hemifield_velocity_outputs,
    }
