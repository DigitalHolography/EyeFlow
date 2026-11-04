"""Importable workers used to exercise actual Windows spawn and queue transport."""

import os
from pathlib import Path
from time import monotonic, sleep

import numpy as np

from input_output.output_manager import OutputType
from pipeline_engine.planning.base import PIPELINE_REGISTRY, pipeline
from pipeline_engine.stage_graph import StageGraph, StageSpec, ProductSpec


@pipeline("process_fixture", execution_mode="process", run_once=True)
def process_fixture(ctx):
    def source(ctx):
        values = ctx.inputs.hd.dataset("moment0")[:]
        assert not values.flags.writeable
        ctx.log("Process fixture source is ready.")
        marker = ctx.inputs.hd.h5.h5file.attrs.get("marker")
        if marker:
            Path(marker + "-started").write_text("started")
            deadline = monotonic() + 20
            while not Path(marker).exists():
                if monotonic() > deadline:
                    raise TimeoutError("Coordinator did not execute independent work")
                sleep(0.01)
        ctx.state.set("from_worker", 7)
        ctx.dag.handoff("video", values + 2)

    def independent(ctx):
        ctx.dag.handoff("independent", np.array([17]))
        artifact = ctx.output.path_for(OutputType.PNG, "worker.png")
        artifact.parent.mkdir(parents=True, exist_ok=True)
        artifact.write_bytes(b"fixture artifact")

    def optional_failure(ctx):
        if ctx.inputs.hd.h5.h5file.attrs.get("crash", False):
            os._exit(23)
        raise ValueError("optional stage failed")

    return StageGraph(
        (
            StageSpec(
                "source", source, products=(ProductSpec("video", path="Worker/video"),)
            ),
            StageSpec("optional", optional_failure, requires=("video",)),
            StageSpec(
                "independent",
                independent,
                products=(ProductSpec("independent", path="Worker/independent"),),
            ),
        )
    )


PROCESS_DESCRIPTOR = PIPELINE_REGISTRY["process_fixture"]


@pipeline("stage_completion_fixture", execution_mode="process", run_once=True)
def stage_completion_fixture(ctx):
    def publish(ctx):
        ctx.dag.handoff("value", 3)
        ctx.state.set("completed_stage", True)
        artifact = ctx.output.path_for(OutputType.PNG, "completed.png")
        artifact.parent.mkdir(parents=True, exist_ok=True)
        artifact.write_bytes(b"completed stage")

    def crash(ctx):
        os._exit(23)

    return StageGraph((
        StageSpec("stage_completion_fixture", publish,
                     products=(ProductSpec("value", path="Stage/value"),)),
        StageSpec("crash", crash, after=("stage_completion_fixture",)),
    ))


STAGE_COMPLETION_DESCRIPTOR = PIPELINE_REGISTRY["stage_completion_fixture"]
