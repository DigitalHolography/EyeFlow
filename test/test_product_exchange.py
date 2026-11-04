"""Logical identity, terminal readiness and coordinator ownership boundaries."""

import json
from dataclasses import replace
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pytest

from pipeline_engine.base import PipelineDescriptor
from pipeline_engine.dag import PipelineDAG
from pipeline_engine.output_model import OutputModel, Record
from pipeline_engine.product_exchange import PipelineInvocationId, ProductId, ProductRequirement, ProductExchangeRegistry, ProductStatus
from pipeline_engine.product_exchange.acceptance import ProductAcceptance, accept_result
from pipeline_engine.protocol import InvocationResult
from pipeline_engine.stage_graph import ProductSpec, StageGraph, StageSpec
from pipeline_engine.stage_graph.journal import StageJournal
from pipeline_engine.storage import ProductStorage


def test_identity_variants_and_finalization_do_not_overwrite_accepted_products():
    registry = ProductExchangeRegistry()
    moments = PipelineInvocationId("run", "velocity", "moments")
    bands = PipelineInvocationId("run", "velocity", "bands")
    x, y = ProductId(moments, "velocity_map"), ProductId(bands, "velocity_map")
    missing = ProductId(moments, "missing")
    for key in (x, y, missing):
        registry.declare(key)
    registry.add_dependent(x, bands)
    assert registry.dependents(x) == {bands}
    assert registry.status(x) == ProductStatus.PENDING
    registry.accept(x, "handle", "commit")
    registry.accept(x, "handle", "commit")
    with pytest.raises(ValueError, match="Conflicting"):
        registry.accept(x, "different", "commit")
    registry.finalize_producer(moments, "crash")
    assert registry.status(x) == ProductStatus.AVAILABLE
    assert registry.status(missing) == ProductStatus.FAILED
    assert registry.status(y) == ProductStatus.PENDING
    same = ProductRequirement(alias="input", producer_pipeline="velocity", product="velocity_map", variant="same")
    assert same.resolve(bands) == y
    assert replace(same, variant="shared").resolve(bands).producer.variant is None


def test_global_product_requirements_select_producer_without_whole_invocation_dependency():
    def graph(*, options=()):
        return StageGraph((StageSpec("stage", lambda ctx: None, products=(ProductSpec("x", exposed=True),)),))
    producer = PipelineDescriptor("a", "a", True, graph_factory=graph)
    consumer = PipelineDescriptor("b", "b", True, product_requires=(
        ProductRequirement(alias="input", producer_pipeline="a", product="x"),))
    dag = PipelineDAG([consumer, producer])
    assert dag.resolve_targets(["b"]).names == ("a", "b")
    assert dag.dependencies_of("b") == ("a",)
    assert dag.dependencies_of("b", include_products=False) == ()
    with pytest.raises(ValueError, match="exposed product"):
        PipelineDAG([producer, replace(consumer, product_requires=(replace(consumer.product_requires[0], product="wrong"),))])


def test_journal_merges_stages_and_variant_scopes_independently(tmp_path):
    with ProductStorage() as storage:
        model = OutputModel(storage, tmp_path / "final.h5")
        StageJournal.observe(model, "science [moments]", "a", "succeeded")
        StageJournal.observe(model, "science [bands]", "a", "running")
        StageJournal.observe(model, "other", "b", "succeeded")
        StageJournal.observe(model, "science [moments]", "c", "pending")
        StageJournal.finalize(model, "science [bands]", "crash")
        values = json.loads(model.attrs["stage_status"])
        assert values == {"science [moments]": {"a": "succeeded", "c": "pending"},
                          "science [bands]": {"a": "failed"}, "other": {"b": "succeeded"}}
        with pytest.raises(ValueError, match="another invocation"):
            StageJournal.reconcile(model, {"stage_status": json.dumps({"other": {"b": "failed"}})},
                                   allowed_scope="science [moments]")


def test_authoritative_acceptance_rejects_another_thread(tmp_path):
    with ProductStorage() as storage:
        model = OutputModel(storage, tmp_path / "final.h5")
        acceptance = ProductAcceptance(model=model, state={}, destination=tmp_path,
                                       inputs={}, registry=ProductExchangeRegistry())
        with ThreadPoolExecutor(max_workers=1) as pool:
            with pytest.raises(RuntimeError, match="coordinator thread"):
                pool.submit(acceptance.finalize, PipelineInvocationId("run", "a")).result()


def test_log_callback_failure_does_not_invalidate_committed_records(tmp_path, monkeypatch):
    def fail_log(message):
        raise RuntimeError("presentation failure")

    monkeypatch.setattr("pipeline_engine.product_exchange.acceptance.Logger.log", fail_log)
    result = InvocationResult(
        records={"": Record(), "value": Record("dataset", np.array([7]))},
        state={"updated": 10}, artifacts=(), artifact_root=tmp_path,
        messages=("presentation callback",),
    )
    with ProductStorage() as storage:
        model = OutputModel(storage, tmp_path / "final.h5")
        state = {}
        accept_result(result, model=model, state=state, destination=tmp_path, name="worker")
        np.testing.assert_array_equal(model["value"][:], [7])
        assert state["updated"] == 10


@pytest.mark.parametrize("scope", ["same", "moments", "unproduced"])
def test_runtime_same_variant_products_are_qualified_and_remain_isolated(tmp_path, scope):
    from types import SimpleNamespace
    import h5py
    from pipeline_engine.base import ExecutionVariant
    from pipeline_engine.runtime import run_pipelines_to_output
    from input_output.output_manager import OutputManager
    from input_output.runs.layout import HoloRunLayout

    def workflows(ctx):
        ctx.publish_execution_variants([
            ExecutionVariant("moments", {"number": 3}, "ProcessingMoments", "moments", {}),
            ExecutionVariant("bands", {"number": 7}, "ProcessingBandRatio", "bands", {}),
        ])
    def produce(ctx):
        ctx.dag.handoff("x", ctx.state["number"])
    def consume(ctx):
        ctx.dag.handoff("value", ctx.dag.require("input") + 1)
    def producer_graph(*, options=()):
        return StageGraph((StageSpec("produce", produce, products=(ProductSpec("x", exposed=True),)),))
    def consumer_graph(*, options=()):
        return StageGraph((StageSpec("consume", consume, requires=("input",),
            products=(ProductSpec("value", "ProcessingMoments/Value"),)),))
    def descriptor(name, callback, **metadata):
        return PipelineDescriptor(name, name, True, pipeline_factory=lambda: SimpleNamespace(name=name, run=callback), **metadata)
    final = run_pipelines_to_output(output_manager=OutputManager(HoloRunLayout.from_holo(tmp_path / "scan.holo")),
        pipelines=[descriptor("source", workflows, produces_execution_variants=True),
                   descriptor("producer", lambda ctx: producer_graph(), graph_factory=producer_graph, dag_requires=("source",)),
                   descriptor("consumer", lambda ctx: consumer_graph(), graph_factory=consumer_graph,
                       run_once=scope != "same",
                       product_requires=(ProductRequirement(alias="input", producer_pipeline="producer", product="x", variant=scope),))],
        holodoppler_h5=None, doppler_vision_h5=None)
    with h5py.File(final) as result:
        if scope == "unproduced":
            assert "ProcessingMoments/Value" not in result
            assert json.loads(result.attrs["pipeline_status"])["shared"]["consumer"] == "skipped"
        else:
            assert result["ProcessingMoments/Value"][()] == 4
        if scope == "same":
            assert result["ProcessingBandRatio/Value"][()] == 8
