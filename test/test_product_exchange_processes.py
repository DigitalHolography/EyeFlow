"""Minimal real-process coverage for product transport and acceptance."""

from app_settings import UserSettings

import json
import multiprocessing
from contextlib import ExitStack
from dataclasses import replace
from pathlib import Path
from time import monotonic
from types import SimpleNamespace

import h5py
import numpy as np
import pytest

from input_output.output_manager import OutputManager, OutputType
from input_output.runs.layout import HoloRunLayout
from pipeline_engine.data.storage import ProductStorage
from pipeline_engine.execution.process_executor import ActiveInvocation, ProcessExecutor
from pipeline_engine.execution.protocol import InvocationSettings
from pipeline_engine.execution.runtime import PipelineRunError, run_pipelines_to_output
from pipeline_engine.io.input.service import InputService
from pipeline_engine.output.output_model import OutputModel
from pipeline_engine.planning.base import PipelineDescriptor
from pipeline_engine.product_exchange import (
    PipelineInvocationId,
    ProductExchangeRegistry,
    ProductId,
    ProductStatus,
)
from pipeline_engine.product_exchange.acceptance import ProductAcceptance
from pipeline_engine.resources import ResourceBudget
from product_exchange_fixtures import (
    ARRAY_SIZE,
    CONSUMER,
    PRODUCER,
    artifact_failure_worker,
)


pytestmark = pytest.mark.usefixtures("inline_output_service")


def test_accepted_product_survives_producer_crash_and_is_published(tmp_path):
    context = multiprocessing.get_context("spawn")
    output = OutputManager(HoloRunLayout.from_holo(tmp_path / "scan.holo"))
    final = output.path_for(OutputType.H5)
    retained = []

    with context.Manager() as manager, ExitStack() as cleanup:
        sync = {
            "received": manager.Event(),
            "consumer_started": manager.Event(),
            "pid": manager.Value("i", 0),
        }

        def initial(ctx):
            ctx.state.set("sync", sync)
            ctx.state.set("final_path", str(final))

        descriptors = [PipelineDescriptor(
            "initial",
            "initial",
            True,
            pipeline_factory=lambda: SimpleNamespace(name="initial", run=initial),
        )]
        descriptors += [
            replace(PRODUCER, dag_requires=("initial",)),
            replace(CONSUMER, dag_requires=("initial",)),
        ]
        budget = ResourceBudget(cpu_threads=8, memory_bytes=2 * 1024**3)

        def retain(pending):
            retained.append(pending)
            cleanup.callback(pending.close)

        with budget.activate(), pytest.raises(PipelineRunError, match="23"):
            run_pipelines_to_output(
                output_manager=output,
                pipelines=descriptors,
                holodoppler_h5=None,
                dopplerview_h5=None,
                settings=UserSettings(values={}),
                on_finalization_ready=retain,
            )
        assert sync["received"].is_set() and sync["consumer_started"].is_set()
        pending = retained[0]
        root = pending.model.storage.directory
        assert root.is_dir() and not final.exists()
        pending.finish()
        assert not root.exists()
        assert budget.cpu_used == 0 and not any(budget.used.values())

    with h5py.File(final) as result:
        statuses = json.loads(result.attrs["pipeline_status"])["shared"]
        assert statuses["exchange_producer"] == "failed"
        assert statuses["exchange_consumer"] == "succeeded"
        assert "Unaccepted" not in result
        np.testing.assert_array_equal(
            result["Accepted/x"][:], np.arange(ARRAY_SIZE, dtype=np.float32)
        )
        assert result["Consumer/sum"][()] == np.float32(
            np.arange(ARRAY_SIZE, dtype=np.float64).sum()
        )


def test_failed_artifact_acceptance_keeps_prior_commit_and_releases_resources(
    tmp_path,
):
    context = multiprocessing.get_context("spawn")
    manager = OutputManager(HoloRunLayout.from_holo(tmp_path / "scan.holo"))
    destination = manager.layout.ef_dir
    destination.mkdir(parents=True)
    (destination / "blocked").write_bytes(b"existing file")
    budget = ResourceBudget(cpu_threads=4, memory_bytes=1024**3)
    registry = ProductExchangeRegistry()
    identity = PipelineInvocationId("transport-run", "transport")
    x_id, y_id = ProductId(identity, "x"), ProductId(identity, "y")
    registry.declare(x_id)
    registry.declare(y_id)

    with budget.activate(), ProductStorage(budget=budget) as storage:
        model = OutputModel(storage, manager.path_for(OutputType.H5))
        service = InputService({}, storage, context)
        state = {"stable": 7}
        acceptance = ProductAcceptance(
            model=model,
            state=state,
            destination=destination,
            inputs={},
            registry=registry,
        )
        executor = ProcessExecutor(
            context=context,
            input_service=service,
            inputs={},
            model=model,
            manager=manager,
            budget=budget,
            state={},
            settings=InvocationSettings({}, {}, (), ()),
            acceptance=acceptance,
            registry=registry,
            run_id=identity.run_id,
        )
        acknowledgements = context.Queue(maxsize=1)
        executor._acks[identity.pipeline] = acknowledgements
        request = budget.default_request(threads=2)
        assert budget.try_acquire(request.threads, request.memory)
        account = budget.memory_account()
        directory = storage.directory / "transport"
        directory.mkdir()
        spec = SimpleNamespace(
            name=identity.pipeline,
            invocation_id=identity,
            owner_id=str(storage.directory.resolve()),
            directory=directory,
            acknowledgements=acknowledgements,
        )
        process = context.Process(
            target=artifact_failure_worker,
            args=(spec, executor.events),
        )
        process.start()
        executor.active[identity.pipeline] = ActiveInvocation(
            process,
            request.threads,
            request.memory,
            account,
        )
        outcomes = []
        try:
            deadline = monotonic() + 25
            while not outcomes:
                assert monotonic() < deadline, "Worker acceptance deadlocked"
                outcomes.extend(executor.collect(wait=True))
        finally:
            executor.close()
            service.close()

        assert registry.status(x_id) == ProductStatus.AVAILABLE
        assert registry.status(y_id) == ProductStatus.FAILED
        assert "Coordinator acceptance failed" in outcomes[0].error
        assert budget.cpu_used == 0 and budget.used["transient"] == 0
        assert not (destination / "first.png").exists()
        assert "Unaccepted" not in model
        assert state == {"stable": 7}
        assert "worker_attribute" not in model.attrs
        assert not list(destination.glob(".ef-*"))
        assert Path(registry.get(x_id).storage.path).is_file()
        model.publish()
        with h5py.File(model.filename) as published:
            np.testing.assert_array_equal(
                published["Accepted/x"][:],
                np.arange(ARRAY_SIZE, dtype=np.float32),
            )
            assert "Unaccepted" not in published

    assert not any(budget.used.values())
