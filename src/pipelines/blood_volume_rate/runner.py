"""Orchestrate selected blood-volume-rate output families."""

from __future__ import annotations

from pipelines.spatial_gradient_moment0.runner import (
    SPATIAL_GRADIENT_PRODUCTS_STATE,
    SpatialGradientProducts,
)
from pipelines.topology_core.runner import prepared_topologies
from pipelines.waveform_velocity_core.runner import (
    VELOCITY_PER_BEAT_OUTPUTS_STATE,
    WAVEFORM_CONTEXT_STATE,
)

from .outputs import pack_gradient_edge_outputs, pack_mask_derived_outputs


def run_blood_volume_rate(ctx) -> dict[str, object]:
    selected = ctx.options_for("blood_volume_rate")
    if not selected:
        return {}

    context = ctx.state.get(WAVEFORM_CONTEXT_STATE)
    if context is None:
        raise RuntimeError("Blood-volume rate requires waveform velocity state.")

    outputs: dict[str, object] = {}
    if "gradient_edges" in selected:
        gradients = ctx.state.get(SPATIAL_GRADIENT_PRODUCTS_STATE)
        if not isinstance(gradients, SpatialGradientProducts):
            raise RuntimeError(
                "Gradient-derived blood-volume rate requires spatial-gradient state."
            )
        outputs.update(
            pack_gradient_edge_outputs(
                context.artery_segment_result,
                context.vein_segment_result,
                gradients,
                context.per_beat_analysis.cycle_boundary_indexes,
                index_base=int(context.source_data.provenance["beat_index_base"]),
            )
        )

    if "masked_edges" in selected:
        velocity_outputs = ctx.state.get(VELOCITY_PER_BEAT_OUTPUTS_STATE)
        if not isinstance(velocity_outputs, dict):
            raise RuntimeError(
                "Mask-derived blood-volume rate requires per-beat velocity state."
            )
        outputs.update(
            pack_mask_derived_outputs(
                prepared_topologies(ctx),
                velocity_outputs,
                pixel_size_mm=float(
                    context.source_data.cross_section_settings.pixel_size_mm
                ),
            )
        )
    return outputs


__all__ = ["run_blood_volume_rate"]
