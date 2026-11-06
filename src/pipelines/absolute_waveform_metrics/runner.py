"""Executable global, segment and quadrant waveform metrics."""

from functools import partial
from dataclasses import replace

from data_manager import InvocationResource
from pipeline_engine import ProductSpec, StageGraph, StageSpec
from pipeline_engine.context import require_managed_runtime
from pipelines.velocity_analysis import velocity_analysis

from .outputs import pack_absolute_waveform_outputs
from .quadrants import pack_quadrant_metrics


def _inputs(ctx):
    analysis = velocity_analysis(ctx)
    values = analysis.waveform_inputs
    ctx.dag.handoff("inputs", (replace(analysis, source_data=InvocationResource(analysis.source_data)), values))


def _selected_metrics(values, ctx):
    return {key: value for key, value in values.items()
            if "segments" in ctx.options_for("absolute_waveform_metrics") or "/by_segment/" not in key}


def _metrics(branch, ctx):
    _, values = ctx.dag.require("inputs")
    computed = pack_absolute_waveform_outputs(
        values, include_segments=branch == "segment_metrics",
        include_global=branch == "global_metrics", include_quadrants=False,
    )
    computed = ctx.dag.handoff(branch, computed)
    ctx.state.set("absolute_waveform_metric_outputs", {
        **ctx.state.get("absolute_waveform_metric_outputs", {}),
        **_selected_metrics(computed, ctx),
    })


def _quadrants(ctx):
    analysis, inputs = ctx.dag.require("inputs")
    values = {**inputs, **ctx.dag.require("global_metrics"), **ctx.dag.require("segment_metrics")}
    from input_output.schema import EyeFlowOutputPaths
    outputs = pack_quadrant_metrics(values, analysis.artery_segments, analysis.vein_segments, EyeFlowOutputPaths.active())
    outputs = ctx.dag.handoff("quadrants", outputs)
    ctx.state.set("absolute_waveform_metric_outputs", {**ctx.state["absolute_waveform_metric_outputs"], **outputs})


def absolute_waveform_graph(*, options=()):
    segments = bool({"segments", "quadrants"} & set(options))
    stages = [
        StageSpec("inputs", _inputs, products=(ProductSpec("inputs"),)),
        StageSpec("global_metrics", partial(_metrics, "global_metrics"), requires=("inputs",),
                  products=(ProductSpec("global_metrics", pack=_selected_metrics),)),
    ]
    if segments:
        stages.append(StageSpec("segment_metrics", partial(_metrics, "segment_metrics"), requires=("inputs",),
            products=(ProductSpec("segment_metrics", pack=_selected_metrics),)))
    if "quadrants" in options:
        stages.append(StageSpec("quadrants", _quadrants, requires=("inputs", "global_metrics", "segment_metrics"),
            products=(ProductSpec("quadrants", pack=lambda value, ctx: value),)))
    return StageGraph(tuple(stages))


def run_absolute_waveform_metrics(ctx):
    require_managed_runtime(ctx, "absolute_waveform_metrics")
    return absolute_waveform_graph(options=ctx.options_for("absolute_waveform_metrics"))
