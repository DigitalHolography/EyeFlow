"""Committed handoffs, dependency failures and workflow-scoped stage graphs."""

from concurrent.futures import ThreadPoolExecutor
import json
from types import SimpleNamespace

import h5py
import numpy as np
import pytest

from input_output.runs.layout import HoloRunLayout
from input_output.output_manager import OutputManager, OutputType
from pipeline_engine import (
    ExecutionVariant, PipelineDescriptor, PipelineGraph, PipelineNode, PipelineProduct,
    current_pipeline_dag, handoff,
)
from pipeline_engine.runtime import PipelineRunError, run_pipelines_to_output


def _descriptor(name, run, **kwargs):
    return PipelineDescriptor(
        name, name, True, pipeline_factory=lambda: SimpleNamespace(name=name, run=run, source_path=__file__),
        **kwargs,
    )


def _manager(root):
    return OutputManager(HoloRunLayout.from_holo(root / "scan.holo"))


def _run(manager, descriptors):
    return run_pipelines_to_output(
        output_manager=manager, pipelines=descriptors, holodoppler_h5=None, doppler_vision_h5=None,
    )


def test_committed_handoff_survives_later_failure_and_only_missing_consumers_skip(tmp_path):
    calls = []

    def produce(ctx):
        assert current_pipeline_dag() is ctx.dag
        values = np.array([3., 4.])
        handoff("metric", values)
        values[:] = 99
        raise ValueError("later optional calculation failed")

    def consume(ctx):
        calls.append("consumer")
        values = ctx.dag.require("metric")
        np.testing.assert_array_equal(values, [3., 4.])
        values[:] = 200
        np.testing.assert_array_equal(ctx.dag.require("metric"), [3., 4.])
        ctx.dag.handoff("derived", values + 1)

    def missing(ctx):
        pytest.fail("The missing product consumer should not run")

    graph = PipelineGraph((
        PipelineNode("producer", produce, outputs=(
            PipelineProduct("metric", "Scientific/Metric", {"unit": "Hz", "dimDesc": ["time"]}),
            PipelineProduct("unavailable", "Scientific/Missing"),
        )),
        PipelineNode("consumer", consume, requires=("metric",), outputs=(PipelineProduct("derived", "Scientific/Derived"),)),
        PipelineNode("missing", missing, requires=("unavailable",)),
        PipelineNode("independent", lambda ctx: calls.append("independent")),
    ))
    with pytest.raises(PipelineRunError, match="later optional calculation failed") as error:
        _run(_manager(tmp_path), [_descriptor("science", lambda ctx: graph)])
    with h5py.File(error.value.output_h5_path) as h5:
        np.testing.assert_array_equal(h5["Scientific/Metric"][:], [3., 4.])
        np.testing.assert_array_equal(h5["Scientific/Derived"][:], [201., 201.])
        assert h5["Scientific/Metric"].attrs["unit"] == "Hz"
        assert "Scientific/Missing" not in h5
        assert json.loads(h5.attrs["stage_status"])["science"] == {
            "producer": "failed", "consumer": "succeeded", "missing": "skipped", "independent": "succeeded",
        }
    assert set(calls) == {"consumer", "independent"}
    with pytest.raises(RuntimeError, match="no pipeline DAG"):
        current_pipeline_dag()


def test_failed_packing_is_atomic_even_when_caught_and_direct_writes_are_rejected(tmp_path):
    def produce(ctx):
        with pytest.raises(TypeError):
            ctx.dag.handoff("bad", 1)
        with pytest.raises(RuntimeError, match="handoff"):
            ctx.output.h5.write("Bypass/value", 1)
        ctx.dag.handoff("good", 2)
        with pytest.raises(ValueError, match="duplicate"):
            ctx.dag.handoff("good", 3)
        with pytest.raises(ValueError, match="Undeclared"):
            ctx.dag.handoff("undeclared", 3)

    graph = PipelineGraph((PipelineNode("producer", produce, outputs=(
        PipelineProduct("bad", pack=lambda value, ctx: {"Partial/value": 1, "Invalid/value": object()}, optional=True),
        PipelineProduct("good", "Good/value"),
    )),))
    path = _run(_manager(tmp_path), [_descriptor("science", lambda ctx: graph)])
    with h5py.File(path) as h5:
        assert h5["Good/value"][()] == 2
        assert "Partial" not in h5 and "Invalid" not in h5 and "Bypass" not in h5


def test_moments_and_bandratio_graphs_keep_handoffs_and_failures_isolated(tmp_path):
    def source(ctx):
        ctx.publish_execution_variants([
            ExecutionVariant(name, {"number": value}, root, folder, {})
            for name, value, root, folder in (
                ("moments", 7, "ProcessingMoments", "moments"),
                ("bands", 9, "ProcessingBandRatio", "bandratio"),
            )
        ])

    def produce(ctx):
        ctx.dag.handoff("metric", ctx.state["number"])
        if ctx.execution_variant.name == "moments":
            raise ValueError("moments optional failure")

    def consume(ctx):
        ctx.dag.handoff("derived", ctx.dag.require("metric") + 1)

    graph = PipelineGraph((
        PipelineNode("producer", produce, outputs=(PipelineProduct("metric", "ProcessingMoments/Metric/value", {"unit": "1"}),)),
        PipelineNode("consumer", consume, requires=("metric",), outputs=(PipelineProduct("derived", "ProcessingMoments/Derived/value"),)),
    ))
    path = _run(_manager(tmp_path), [
        _descriptor("source", source, produces_execution_variants=True),
        _descriptor("science", lambda ctx: graph, dag_requires=("source",)),
    ])
    with h5py.File(path) as h5:
        assert h5["ProcessingMoments/Metric/value"][()] == 7
        assert h5["ProcessingBandRatio/Metric/value"][()] == 9
        assert h5["ProcessingMoments/Derived/value"][()] == 8
        assert h5["ProcessingBandRatio/Derived/value"][()] == 10
        statuses = json.loads(h5.attrs["stage_status"])
        assert statuses["science [moments]"]["producer"] == "failed"
        assert statuses["science [bands]"]["producer"] == "succeeded"
        products = json.loads(h5.attrs["stage_products"])
        assert products["science [moments]"]["metric"]["datasets"] == ["ProcessingMoments/Metric/value"]
        assert products["science [bands]"]["metric"]["datasets"] == ["ProcessingBandRatio/Metric/value"]


def test_graph_definition_rejects_invalid_dependencies():
    noop = lambda ctx: None
    graph = PipelineGraph((
        PipelineNode("a", noop, after=("b",)),
        PipelineNode("b", noop, after=("a",)),
    ))
    with pytest.raises((ValueError, RuntimeError)):
        graph.plan()


def test_graph_report_waits_for_declared_artifacts_and_writes_intermediates_immediately(tmp_path, monkeypatch):
    monkeypatch.setattr("pipeline_engine.background_artifacts.ProcessPoolExecutor", lambda **kwargs: ThreadPoolExecutor(max_workers=1))
    manager = _manager(tmp_path)

    def produce(ctx):
        ctx.dag.handoff("metric", 7)

    def artifact(ctx):
        assert ctx.dag.require("metric") == 7
        ctx.output.write_png(np.ones((2, 2), dtype=np.uint8), "upstream.png")

    def report(ctx):
        assert manager.path_for(OutputType.PNG, "upstream.png").is_file()
        intermediate = ctx.output.write_png(np.zeros((2, 2), dtype=np.uint8), "report.png")
        assert intermediate.is_file()

    graph = PipelineGraph((
        PipelineNode("producer", produce, outputs=(PipelineProduct("metric"),)),
        PipelineNode("figures", artifact, requires=("metric",), artifact=True),
        PipelineNode("report", report, after=("figures",), wait_for_artifacts=True),
    ))
    _run(manager, [_descriptor("science", lambda ctx: graph)])
