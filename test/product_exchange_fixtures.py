"""Importable real-process fixtures, synchronized by manager event proxies."""

import os
from pathlib import Path

import numpy as np
import psutil

from pipeline_engine.base import PIPELINE_REGISTRY, pipeline
from pipeline_engine.product_exchange import ProductRequirement
from pipeline_engine.stage_graph import ProductSpec, StageGraph, StageSpec

ARRAY_SIZE = 2 * 1024 * 1024


def producer_graph(*, options=()):
    return StageGraph((StageSpec("produce", produce, products=(
        ProductSpec("x", "Accepted/x", exposed=True), ProductSpec("never", "Accepted/tiny", optional=True, exposed=True),
    )),))


def produce(ctx):
    sync = ctx.state["sync"]
    if "initial_array" in ctx.state:
        assert not ctx.state["initial_array"].flags.writeable
        np.testing.assert_array_equal(ctx.state["initial_array"], 3)
    sync["pid"].value = os.getpid()
    if ctx.state.get("mode") == "before_emit":
        os._exit(23)
    ctx.dag.handoff("x", np.arange(ARRAY_SIZE, dtype=np.float32))
    if ctx.state.get("mode") == "small_product":
        ctx.dag.handoff("never", np.asarray([3, 7], dtype=np.int32))
    # A transport receipt has arrived; authoritative acceptance may still be queued.
    sync["received"].set()
    if ctx.state.get("mode") == "after_receipt":
        assert sync["crash"].wait(20)
    else:
        assert sync["b_started"].wait(20), "B waited for producer invocation completion"
    ctx.runtime.work_h5.create_dataset("Unaccepted/value", data=999)
    ctx.state.set("unaccepted", 999)
    os._exit(23)


@pipeline("exchange_producer", execution_mode="process", run_once=True, graph_factory=producer_graph)
def producer(ctx):
    return producer_graph()


def consumer_graph(*, options=()):
    return StageGraph((
        StageSpec("missing", lambda ctx: None, products=(ProductSpec("missing", optional=True),)),
        StageSpec("consume", consume, requires=("input",), products=(ProductSpec("sum", "Consumer/sum"),)),
        StageSpec("optional", optional, optional_requires=("missing",), products=(ProductSpec("independent", "Consumer/independent"),)),
    ))


def consume(ctx):
    sync = ctx.state["sync"]
    values = ctx.dag.require("input")
    if "initial_array" in ctx.state:
        assert not ctx.state["initial_array"].flags.writeable
        np.testing.assert_array_equal(ctx.state["initial_array"], 3)
    assert not values.flags.writeable
    assert not Path(ctx.state["final_path"]).exists()
    np.testing.assert_array_equal(values, np.arange(ARRAY_SIZE, dtype=np.float32))
    sync["backing"].value = str(values.filename)
    sync["b_started"].set()
    # Wait for an actual producer exit, not an arbitrary delay.
    try:
        psutil.Process(sync["pid"].value).wait(timeout=20)
    except psutil.NoSuchProcess:
        pass
    assert Path(values.filename).is_file()
    np.testing.assert_array_equal(values[-1024:], np.arange(ARRAY_SIZE - 1024, ARRAY_SIZE, dtype=np.float32))
    ctx.dag.handoff("sum", np.sum(values, dtype=np.float64))
    sync["b_finished"].set()


def optional(ctx):
    assert ctx.dag.get("missing") is None
    ctx.dag.handoff("independent", 17)


@pipeline("exchange_consumer", execution_mode="process", run_once=True, graph_factory=consumer_graph,
          product_requires=(
              ProductRequirement(alias="input", producer_pipeline="exchange_producer", product="x"),
              # A different pending product would intentionally block B until A
              # finishes. The optional-failure branch is exercised separately.
          ))
def consumer(ctx):
    return consumer_graph()


PRODUCER = PIPELINE_REGISTRY["exchange_producer"]
CONSUMER = PIPELINE_REGISTRY["exchange_consumer"]


def _raise_on_coordinator():
    raise RuntimeError("coordinator acceptance fixture exception")


class RaisingState:
    def __reduce__(self):
        return _raise_on_coordinator, ()


def transport_fixture(spec, events):
    """Exercise invalid/retried real manifests through the production receiver."""
    from dataclasses import replace
    from uuid import uuid4
    from data_manager import DataManager
    from pipeline_engine.output_model import OutputModel, Record
    from pipeline_engine.protocol import InvocationResult, ProductAcceptedEvent, WorkerEvent, TransportReceipt
    from pipeline_engine.storage import ProductStorage
    from pipeline_engine.transfer import dump_result

    storage = ProductStorage(directory=spec.directory / "products")
    data = DataManager(storage, owner_id=spec.owner_id)
    model = OutputModel(storage, spec.directory / "unpublished.h5")
    artifacts = spec.directory / "artifacts"
    artifacts.mkdir()

    def send(event):
        events.put(event)
        receipt = spec.acknowledgements.get()
        assert isinstance(receipt, TransportReceipt) and receipt.commit_id == event.commit_id

    def commit(records, *, products=(), state=None, files=(), commit_id=None, kind="ProductCommit"):
        identity = commit_id or uuid4().hex
        result = InvocationResult(records, state or {}, files, artifacts,
                                  invocation_id=spec.invocation_id, commit_id=identity, products=products)
        manifest = dump_result(result, storage)
        event = ProductAcceptedEvent(kind, spec.name, "calculate", manifest,
            commit_id=identity, invocation_id=spec.invocation_id)
        send(event)
        return event

    try:
        send(WorkerEvent("WorkerStarted", spec.name, invocation_id=spec.invocation_id))
        send(WorkerEvent("StageStarted", spec.name, "calculate", invocation_id=spec.invocation_id))
        x = data.accept(np.arange(32768, dtype=np.float32))
        model.create_dataset("Accepted/x", data=data.read(x))
        first = commit(model.committed_records(), products=(("x", x),))
        mode = spec.mode
        if mode == "duplicate":
            send(first)
            commit({"": Record()}, kind="InvocationFinished")
            return
        if mode == "cancel":
            assert spec.hold.wait(30), "Coordinator did not cancel the worker"
            return
        y = data.accept(np.ones(32768, dtype=np.float32))
        records = {"": Record(), "Unaccepted": Record(), "Unaccepted/value": Record("dataset", y.storage)}
        if mode == "invalid_reference":
            y = replace(y, storage=replace(y.storage, path=str(spec.directory / "missing.npy")))
        elif mode == "merge_failure":
            records = {"": Record(), "Accepted": Record(), "Accepted/x": Record("dataset", y.storage)}
        elif mode == "artifact_failure":
            first_file, second_file = artifacts / "first.png", artifacts / "blocked/second.png"
            first_file.write_bytes(b"first")
            second_file.parent.mkdir()
            second_file.write_bytes(b"second")
            commit(records, products=(("y", y),), files=(first_file, second_file))
            return
        elif mode == "acceptance_exception":
            commit(records, products=(("y", y),), state={"raising": RaisingState()})
            return
        elif mode == "identity_mismatch":
            bad = replace(first, invocation_id=replace(spec.invocation_id, variant="other"))
            send(bad)
            return
        commit(records, products=(("y", y),), commit_id=first.commit_id if mode == "conflicting_duplicate" else None)
    finally:
        storage.close()
