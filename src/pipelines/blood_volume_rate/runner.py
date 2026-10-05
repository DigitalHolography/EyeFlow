"""Orchestrate selected blood-volume-rate output families."""

from __future__ import annotations

from input_output.schema import EyeFlowOutputPaths
from pipelines.spatial_gradient_moment0.runner import (
    SPATIAL_GRADIENT_PRODUCTS_STATE,
    SpatialGradientProducts,
)
from pipelines.topology_core.runner import prepared_topologies
from pipelines.waveform_velocity import waveform_velocity
from pipelines.waveform_velocity.outputs import (
    pack_velocity_per_beat_outputs,
)

from .outputs import (
    export_blood_volume_rate_signals,
    export_lumen_diameter_distributions,
    pack_gradient_edge_outputs,
    pack_mask_derived_outputs,
)


def run_blood_volume_rate(ctx) -> dict[str, object]:
    selected = ctx.options_for("blood_volume_rate")
    if not selected:
        return {}

    waveform = waveform_velocity(ctx)

    outputs: dict[str, object] = {}
    if "gradient_edges" in selected:
        gradients = ctx.state.get(SPATIAL_GRADIENT_PRODUCTS_STATE)
        if not isinstance(gradients, SpatialGradientProducts):
            raise RuntimeError(
                "Gradient-derived blood-volume rate requires spatial-gradient state."
            )
        outputs.update(
            pack_gradient_edge_outputs(
                waveform.artery_segments,
                waveform.vein_segments,
                gradients,
                waveform.cycle_boundary_indexes,
                index_base=0,
                gradient_sources_persisted=ctx.pipeline_targeted(
                    "spatial_gradient_moment0"
                ),
            )
        )

    if "masked_edges" in selected:
        velocity_outputs = pack_velocity_per_beat_outputs(
            waveform.per_beat_result,
            velocity_analysis=waveform.retinal_velocity,
        )
        mask_outputs = pack_mask_derived_outputs(
            prepared_topologies(ctx),
            velocity_outputs,
            pixel_size_mm=float(
                waveform.source_data.profile_settings.pixel_size_mm
            ),
        )
        outputs.update(mask_outputs)
        if ctx.output.available:
            schema = EyeFlowOutputPaths.active()
            export_blood_volume_rate_signals(
                ctx.output,
                mask_outputs[schema.blood_volume_rate.artery.total_masked_edges],
                mask_outputs[schema.blood_volume_rate.vein.total_masked_edges],
            )
    if ctx.output.available:
        schema = EyeFlowOutputPaths.active()
        export_lumen_diameter_distributions(
            ctx.output,
            ctx.output.h5.array(schema.segmentation.artery.lumen_diameter),
            ctx.output.h5.array(schema.segmentation.vein.lumen_diameter),
            pixel_pitch_m=ctx.output.h5.read(schema.segmentation.pixel_pitch_m),
        )
    return outputs


__all__ = ["run_blood_volume_rate"]
