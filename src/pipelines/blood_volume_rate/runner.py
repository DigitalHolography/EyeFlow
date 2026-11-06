"""Executable gradient- and mask-derived blood-volume-rate branches."""

from dataclasses import replace

from data_manager import InvocationResource
from input_output.schema import EyeFlowOutputPaths
from pipeline_engine import ProductSpec, StageGraph, StageSpec
from pipeline_engine.context import require_managed_runtime
from pipelines.spatial_gradient_moment0.runner import (
    SPATIAL_GRADIENT_PRODUCTS_STATE,
    SpatialGradientProducts,
)
from pipelines.topology_core.runner import prepared_topologies
from pipelines.velocity_analysis import velocity_analysis

from .outputs import (export_blood_volume_rate_signals, export_lumen_diameter_distributions,
                      pack_gradient_edge_outputs, pack_mask_derived_outputs)


def _inputs(ctx):
    analysis = velocity_analysis(ctx)
    ctx.dag.handoff("analysis", replace(analysis, source_data=InvocationResource(analysis.source_data)))


def _gradient_edges(ctx):
    analysis = ctx.dag.require("analysis")
    gradients = ctx.state.get(SPATIAL_GRADIENT_PRODUCTS_STATE)
    if not isinstance(gradients, SpatialGradientProducts):
        raise RuntimeError("Gradient-derived blood-volume rate requires spatial-gradient state.")
    outputs = pack_gradient_edge_outputs(analysis.artery_segments, analysis.vein_segments,
        gradients, analysis.cycle_boundary_indexes, index_base=0,
        gradient_sources_persisted=ctx.pipeline_targeted("spatial_gradient_moment0"))
    ctx.dag.handoff("gradient_edges", outputs)


def _masked_edges(ctx):
    analysis = ctx.dag.require("analysis")
    values = analysis.waveform_inputs
    ctx.dag.handoff("masked_edges", pack_mask_derived_outputs(prepared_topologies(ctx), values,
        pixel_size_mm=float(analysis.source_data.profile_settings.pixel_size_mm)))


def _signals(ctx):
    if not ctx.output.available:
        return
    paths = EyeFlowOutputPaths.active().blood_volume_rate
    values = ctx.dag.require("masked_edges")
    export_blood_volume_rate_signals(ctx.output,
        values[paths.artery.total_masked_edges], values[paths.vein.total_masked_edges])


def _lumen_distributions(ctx):
    if not ctx.output.available:
        return
    schema = EyeFlowOutputPaths.active()
    export_lumen_diameter_distributions(ctx.output,
        ctx.output.h5.array(schema.segmentation.artery.lumen_diameter),
        ctx.output.h5.array(schema.segmentation.vein.lumen_diameter),
        pixel_pitch_m=ctx.output.h5.read(schema.segmentation.pixel_pitch_m))


def blood_volume_rate_graph(*, options=()):
    branches = [name for name in ("gradient_edges", "masked_edges") if name in options]
    if not branches:
        return StageGraph(())
    stages = [StageSpec("inputs", _inputs, products=(ProductSpec("analysis"),))]
    callbacks = {"gradient_edges": _gradient_edges, "masked_edges": _masked_edges}
    for name in branches:
        stages.append(StageSpec(name, callbacks[name], requires=("analysis",),
            products=(ProductSpec(name, pack=lambda value, ctx: value),)))
    if "masked_edges" in branches:
        stages.append(StageSpec("signals", _signals, requires=("masked_edges",), artifact=True))
    stages.append(StageSpec("lumen_distributions", _lumen_distributions,
                            after=tuple(branches), artifact=True))
    return StageGraph(tuple(stages))


def run_blood_volume_rate(ctx):
    require_managed_runtime(ctx, "blood_volume_rate")
    return blood_volume_rate_graph(options=ctx.options_for("blood_volume_rate"))
