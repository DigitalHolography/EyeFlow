"""Importable workers used to exercise actual Windows spawn and queue transport."""

import os
from pathlib import Path
from time import monotonic, sleep

import numpy as np

from input_output.output_manager import OutputType
from pipeline_engine.base import PIPELINE_REGISTRY, pipeline
from pipeline_engine.graphs import PipelineGraph, PipelineNode, PipelineProduct


@pipeline("process_fixture", execution_mode="process", run_once=True)
def process_fixture(ctx):
    def source(ctx):
        values = ctx.inputs.hd.dataset("moment0")[:]
        assert not values.flags.writeable
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

    return PipelineGraph(
        (
            PipelineNode(
                "source", source, outputs=(PipelineProduct("video", path="Worker/video"),)
            ),
            PipelineNode("optional", optional_failure, requires=("video",)),
            PipelineNode(
                "independent",
                independent,
                outputs=(PipelineProduct("independent", path="Worker/independent"),),
            ),
        )
    )


PROCESS_DESCRIPTOR = PIPELINE_REGISTRY["process_fixture"]
