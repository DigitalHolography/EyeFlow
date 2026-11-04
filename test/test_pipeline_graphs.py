"""Committed handoffs, dependency failures and workflow-scoped stage graphs."""

from concurrent.futures import ThreadPoolExecutor
from copy import copy
import json
from threading import Event
from types import SimpleNamespace

import h5py
import numpy as np
import pytest

from input_output.runs.layout import HoloRunLayout
from input_output.output_manager import OutputManager, OutputType
from pipeline_engine import (
    ExecutionVariant, PipelineDescriptor, StageGraph, StageSpec, ProductSpec,
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
        values = ctx.dag.require_mutable("metric")
        np.testing.assert_array_equal(values, [3., 4.])
        values[:] = 200
        np.testing.assert_array_equal(ctx.dag.require("metric"), [3., 4.])
        ctx.dag.handoff("derived", values + 1)

    def missing(ctx):
        pytest.fail("The missing product consumer should not run")

    graph = StageGraph((
        StageSpec("producer", produce, products=(
            ProductSpec("metric", "Scientific/Metric", {"unit": "Hz", "dimDesc": ["time"]}),
            ProductSpec("unavailable", "Scientific/Missing"),
        )),
        StageSpec("consumer", consume, requires=("metric",), products=(ProductSpec("derived", "Scientific/Derived"),)),
        StageSpec("missing", missing, requires=("unavailable",)),
        StageSpec("independent", lambda ctx: calls.append("independent")),
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


def test_handoff_launches_ready_artifacts_before_its_producer_returns(tmp_path, monkeypatch):
    monkeypatch.setattr("pipeline_engine.background_artifacts.ProcessPoolExecutor", lambda **kwargs: ThreadPoolExecutor(max_workers=1))
    manager = _manager(tmp_path)
    started = Event()

    def export(snapshot, output):
        started.set()
        return [output.write_png(np.zeros((1, 1), dtype=np.uint8), "snapshot.png")]

    def produce(ctx):
        ctx.dag.handoff("metric", 7)
        assert started.wait(10), "Export must start before the producer finishes"
        ctx._background_artifacts.collect(wait=True)
        assert manager.path_for(OutputType.PNG, "snapshot.png").is_file()
        raise ValueError("failure after committed data")

    def artifacts(ctx):
        assert ctx.dag.require("metric") == 7
        ctx.defer_artifacts("png", export, lambda path: path.write_bytes(b"data"))

    graph = StageGraph((
        StageSpec("producer", produce, products=(ProductSpec("metric", "Scientific/Metric"),)),
        StageSpec("figures", artifacts, requires=("metric",), artifact=True),
    ))
    with pytest.raises(PipelineRunError, match="failure after committed"):
        _run(manager, [_descriptor("science", lambda ctx: graph)])
    assert manager.path_for(OutputType.PNG, "snapshot.png").is_file()


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

    graph = StageGraph((StageSpec("producer", produce, products=(
        ProductSpec("bad", pack=lambda value, ctx: {"Partial/value": 1, "Invalid/value": object()}, optional=True),
        ProductSpec("good", "Good/value"),
    )),))
    path = _run(_manager(tmp_path), [_descriptor("science", lambda ctx: graph)])
    with h5py.File(path) as h5:
        assert h5["Good/value"][()] == 2
        assert "Partial" not in h5 and "Invalid" not in h5 and "Bypass" not in h5


def test_copied_or_worker_context_cannot_handoff_into_another_scope(tmp_path):
    def produce(ctx):
        clone = copy(ctx)
        with pytest.raises(RuntimeError, match="No active pipeline DAG"):
            clone.dag.handoff("metric", 9)
        with ThreadPoolExecutor(max_workers=1) as pool:
            with pytest.raises(RuntimeError, match="owning pipeline thread"):
                pool.submit(ctx.dag.handoff, "metric", 9).result()
            with pytest.raises(RuntimeError, match="no pipeline DAG"):
                pool.submit(handoff, "metric", 9).result()
        ctx.dag.handoff("metric", 4)

    graph = StageGraph((StageSpec("producer", produce, products=(ProductSpec("metric", "Metric/value"),)),))
    path = _run(_manager(tmp_path), [_descriptor("science", lambda ctx: graph)])
    with h5py.File(path) as h5:
        assert h5["Metric/value"][()] == 4


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

    graph = StageGraph((
        StageSpec("producer", produce, products=(ProductSpec("metric", "ProcessingMoments/Metric/value", {"unit": "1"}),)),
        StageSpec("consumer", consume, requires=("metric",), products=(ProductSpec("derived", "ProcessingMoments/Derived/value"),)),
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

    graph = StageGraph((
        StageSpec("producer", produce, products=(ProductSpec("metric"),)),
        StageSpec("figures", artifact, requires=("metric",), artifact=True),
        StageSpec("report", report, after=("figures",), wait_for_artifacts=True),
    ))
    _run(manager, [_descriptor("science", lambda ctx: graph)])


def test_handoff_cannot_replace_a_parent_of_an_already_committed_dataset(tmp_path):
    def produce(ctx):
        ctx.dag.handoff("leaf", 7)
        with pytest.raises(ValueError, match="non-overlapping"):
            ctx.dag.handoff("parent", 9)

    graph = StageGraph((StageSpec("producer", produce, products=(
        ProductSpec("leaf", "Group/Metric"), ProductSpec("parent", pack=lambda value, ctx: {"Group": value}, optional=True),
    )),))
    path = _run(_manager(tmp_path), [_descriptor("science", lambda ctx: graph)])
    with h5py.File(path) as h5:
        assert h5["Group/Metric"][()] == 7


def test_missing_output_optional_resolution_after_success_and_cleanup_once(tmp_path):
    calls, captured = [], []

    def producer(ctx):
        captured.append(ctx.dag)
        with pytest.raises(ValueError, match="did not declare"):
            ctx.dag.get("undeclared")
        # Missing required output fails this stage after it returns normally.

    def optional(ctx):
        assert ctx.dag.get("x") is None
        calls.append("optional")

    graph = StageGraph((
        StageSpec("producer", producer, products=(ProductSpec("x"),)),
        StageSpec("required", lambda ctx: pytest.fail("Missing required product ran"), requires=("x",)),
        StageSpec("optional", optional, optional_requires=("x",)),
        StageSpec("after", lambda ctx: pytest.fail("Failed stage completion was ignored"), after=("producer",)),
        StageSpec("independent", lambda ctx: calls.append("independent")),
    ), cleanup=lambda: calls.append("cleanup"))
    with pytest.raises(PipelineRunError, match="Missing handoffs"):
        _run(_manager(tmp_path), [_descriptor("science", lambda ctx: graph)])
    assert calls == ["optional", "independent", "cleanup"]
    with pytest.raises(RuntimeError, match="completed"):
        captured[0].handoff("x", 3)


def test_optional_artifact_waits_for_producer_resolution(tmp_path):
    calls = []
    def producer(ctx):
        ctx.dag.handoff("x", 3)
        assert not calls, "An optional pending product was treated as absent"
        calls.append("producer")
    def artifact(ctx):
        assert ctx.dag.require("x") == 3
        assert ctx.dag.get("y") is None
        calls.append("artifact")
    graph = StageGraph((
        StageSpec("producer", producer, products=(ProductSpec("x"), ProductSpec("y", optional=True))),
        StageSpec("artifact", artifact, requires=("x",), optional_requires=("y",), artifact=True),
    ))
    _run(_manager(tmp_path), [_descriptor("science", lambda ctx: graph)])
    assert calls == ["producer", "artifact"]
