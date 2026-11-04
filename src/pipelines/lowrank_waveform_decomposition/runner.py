"""Executable low-rank global and optional quadrant decomposition."""

from functools import partial

from data_manager import InvocationResource
from pipeline_engine import ProductSpec, StageGraph, StageSpec
from pipeline_engine.context import require_managed_runtime
from pipelines.velocity_analysis import velocity_analysis

from .outputs import pack_lowrank_waveform_decomposition_outputs

LOWRANK_WAVEFORM_OUTPUTS_STATE = "lowrank_waveform_decomposition_outputs"


def _inputs(ctx):
    analysis = velocity_analysis(ctx)
    values = analysis.waveform_inputs
    ctx.dag.handoff("inputs", InvocationResource((analysis, values)))


def _decompose(quadrants, ctx):
    analysis, values = ctx.dag.require("inputs")
    outputs = pack_lowrank_waveform_decomposition_outputs(values,
        include_global=not quadrants, include_quadrants=quadrants,
        artery_segments=analysis.artery_segments if quadrants else None,
        vein_segments=analysis.vein_segments if quadrants else None)
    ctx.dag.handoff("quadrants" if quadrants else "metrics", outputs)
    ctx.state.set(LOWRANK_WAVEFORM_OUTPUTS_STATE,
        {**ctx.state.get(LOWRANK_WAVEFORM_OUTPUTS_STATE, {}), **outputs})


def lowrank_waveform_graph(*, options=()):
    packed = lambda value, ctx: value
    stages = [StageSpec("inputs", _inputs, products=(ProductSpec("inputs"),)),
              StageSpec("metrics", partial(_decompose, False), requires=("inputs",),
                        products=(ProductSpec("metrics", pack=packed),))]
    if "quadrants" in options:
        stages.append(StageSpec("quadrants", partial(_decompose, True), requires=("inputs",),
                                products=(ProductSpec("quadrants", pack=packed),)))
    return StageGraph(tuple(stages))


def run_lowrank_waveform_decomposition(ctx):
    require_managed_runtime(ctx, "lowrank_waveform_decomposition")
    return lowrank_waveform_graph(options=ctx.options_for("lowrank_waveform_decomposition"))
