"""Executable global, segment and quadrant waveform metrics."""

from functools import partial

from data_manager import InvocationResource
from pipeline_engine import ProductSpec, StageGraph, StageSpec
from pipeline_engine.context import require_managed_runtime
from pipelines.velocity_analysis import velocity_analysis

from .outputs import pack_waveform_shape_outputs
from .metrics.quadrants import pack_quadrant_metrics


def _inputs(ctx):
    analysis = velocity_analysis(ctx)
    values = analysis.waveform_inputs
    ctx.dag.handoff("inputs", InvocationResource((analysis, values)))


def _selected_metrics(values, ctx):
    return {key: value for key, value in values.items()
            if "segments" in ctx.options_for("waveform_shape_metrics") or "/by_segment/" not in key}


def _metrics(branch, ctx):
    analysis, values = ctx.dag.require("inputs")
    computed = pack_waveform_shape_outputs(
        values, analysis.source_data, analysis.artery_segments, analysis.vein_segments,
        include_segments=branch == "segment_metrics",
        include_global=branch == "global_metrics", include_quadrants=False,
    )
    ctx.dag.handoff(branch, InvocationResource(computed))
    ctx.state.set("waveform_shape_metric_outputs", {
        **ctx.state.get("waveform_shape_metric_outputs", {}),
        **_selected_metrics(computed, ctx),
    })


def _quadrants(ctx):
    analysis, inputs = ctx.dag.require("inputs")
    values = {**inputs, **ctx.dag.require("global_metrics"), **ctx.dag.require("segment_metrics")}
    outputs = pack_quadrant_metrics(values, analysis.artery_segments, analysis.vein_segments)
    ctx.dag.handoff("quadrants", outputs)
    ctx.state.set("waveform_shape_metric_outputs", {**ctx.state["waveform_shape_metric_outputs"], **outputs})


def waveform_shape_graph(*, options=()):
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


def run_waveform_shape_metrics(ctx):
    require_managed_runtime(ctx, "waveform_shape_metrics")
    return waveform_shape_graph(options=ctx.options_for("waveform_shape_metrics"))
