"""Small importable fixtures for the real product-transport boundary."""

import os
from pathlib import Path

import numpy as np
import psutil

from pipeline_engine.planning.base import PIPELINE_REGISTRY, pipeline
from pipeline_engine.product_exchange import ProductRequirement
from pipeline_engine.stage_graph import ProductSpec, StageGraph, StageSpec

ARRAY_SIZE = 32768


def producer_graph(*, options=()):
    return StageGraph((StageSpec(
        "produce",
        produce,
        products=(ProductSpec("x", "Accepted/x", exposed=True),),
    ),))


def produce(ctx):
    sync = ctx.state["sync"]
    sync["pid"].value = os.getpid()
    ctx.dag.handoff("x", np.arange(ARRAY_SIZE, dtype=np.float32))
    sync["received"].set()
    assert sync["consumer_started"].wait(20)
    ctx.runtime.work_h5.create_dataset("Unaccepted/value", data=999)
    os._exit(23)


@pipeline(
    "exchange_producer",
    execution_mode="process",
    run_once=True,
    graph_factory=producer_graph,
)
def producer(ctx):
    return producer_graph()


def consumer_graph(*, options=()):
    return StageGraph((StageSpec(
        "consume",
        consume,
        requires=("input",),
        products=(ProductSpec("sum", "Consumer/sum"),),
    ),))


def consume(ctx):
    sync = ctx.state["sync"]
    values = ctx.dag.require("input")
    assert not values.flags.writeable
    assert not Path(ctx.state["final_path"]).exists()
    np.testing.assert_array_equal(values, np.arange(ARRAY_SIZE, dtype=np.float32))
    sync["consumer_started"].set()
    try:
        psutil.Process(sync["pid"].value).wait(timeout=20)
    except psutil.NoSuchProcess:
        pass
    np.testing.assert_array_equal(
        values[-1024:],
        np.arange(ARRAY_SIZE - 1024, ARRAY_SIZE, dtype=np.float32),
    )
    ctx.dag.handoff("sum", np.sum(values, dtype=np.float64))


@pipeline(
    "exchange_consumer",
    execution_mode="process",
    run_once=True,
    graph_factory=consumer_graph,
    product_requires=(ProductRequirement(
        alias="input",
        producer_pipeline="exchange_producer",
        product="x",
    ),),
)
def consumer(ctx):
    return consumer_graph()


PRODUCER = PIPELINE_REGISTRY["exchange_producer"]
CONSUMER = PIPELINE_REGISTRY["exchange_consumer"]


def artifact_failure_worker(spec, events):
    """Publish one valid product, then an atomically rejected artifact update."""

    from uuid import uuid4

    from data_manager import DataManager
    from pipeline_engine.data.storage import ProductStorage
    from pipeline_engine.data.transfer import dump_result
    from pipeline_engine.execution.protocol import (
        InvocationResult,
        ProductAcceptedEvent,
        TransportReceipt,
        WorkerEvent,
    )
    from pipeline_engine.output.output_model import OutputModel, Record

    storage = ProductStorage(directory=spec.directory / "products")
    data = DataManager(storage, owner_id=spec.owner_id)
    model = OutputModel(storage, spec.directory / "unpublished.h5")
    artifacts = spec.directory / "artifacts"
    artifacts.mkdir()

    def send(event):
        events.put(event)
        receipt = spec.acknowledgements.get()
        assert isinstance(receipt, TransportReceipt)
        assert receipt.commit_id == event.commit_id

    def commit(records, *, products=(), state=None, files=()):
        commit_id = uuid4().hex
        result = InvocationResult(
            records,
            state or {},
            files,
            artifacts,
            invocation_id=spec.invocation_id,
            commit_id=commit_id,
            products=products,
        )
        send(ProductAcceptedEvent(
            "ProductCommit",
            spec.name,
            "calculate",
            dump_result(result, storage),
            commit_id=commit_id,
            invocation_id=spec.invocation_id,
        ))

    try:
        send(WorkerEvent(
            "WorkerStarted", spec.name, invocation_id=spec.invocation_id
        ))
        send(WorkerEvent(
            "StageStarted",
            spec.name,
            "calculate",
            invocation_id=spec.invocation_id,
        ))
        x = data.accept(np.arange(ARRAY_SIZE, dtype=np.float32))
        with model.checkpoint(owner=spec.name):
            model.create_dataset("Accepted/x", data=data.read(x))
        commit(model.committed_records(), products=(("x", x),))

        y = data.accept(np.ones(ARRAY_SIZE, dtype=np.float32))
        first = artifacts / "first.png"
        second = artifacts / "blocked/second.png"
        first.write_bytes(b"first")
        second.parent.mkdir()
        second.write_bytes(b"second")
        commit(
            {
                "": Record(attrs={"worker_attribute": 1}),
                "Unaccepted": Record(),
                "Unaccepted/value": Record("dataset", y.storage),
            },
            products=(("y", y),),
            state={"stable": 999},
            files=(first, second),
        )
    finally:
        storage.close()
