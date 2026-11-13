"""Cross-process recovery through actual transport, ownership and publication."""

import json
import multiprocessing
from pathlib import Path
from types import SimpleNamespace

import h5py
import numpy as np
import pytest

from input_output.output_manager import OutputManager, OutputType
from input_output.runs.layout import HoloRunLayout
from pipeline_engine.planning.base import PipelineDescriptor
from pipeline_engine.output.output_model import OutputModel
from pipeline_engine.resources import ResourceBudget
from pipeline_engine.execution.runtime import PipelineRunError, run_pipelines_to_output
from product_exchange_fixtures import ARRAY_SIZE, CONSUMER, PRODUCER, TRANSPORT_FAULTS, fault_worker


@pytest.mark.parametrize("mode", ["overlap", "after_receipt", "before_emit", "before_emit_required_only", "small_product",
                                  "snapshot_failure", "delivery_failure", "missing_receipt"])
def test_real_spawn_product_survives_producer_crash_and_final_publication(tmp_path, monkeypatch, mode):
    mp = multiprocessing.get_context("spawn")
    output = OutputManager(HoloRunLayout.from_holo(tmp_path / "scan.holo"))
    final = output.path_for(OutputType.H5)
    before_emit = mode in {"before_emit", "before_emit_required_only"}
    required_only = mode == "before_emit_required_only"
    from utils.logger import Logger
    logs = []
    monkeypatch.setattr(Logger, "_current", Logger(on_log=logs.append))
    bindings = []
    if mode in TRANSPORT_FAULTS:
        monkeypatch.setattr("pipeline_engine.execution.process_executor.run_worker", fault_worker)
    with mp.Manager() as manager:
        sync = {name: manager.Event() for name in ("received", "crash", "b_started", "b_finished")}
        sync.update(pid=manager.Value("i", 0), backing=manager.Value(str, ""))

        def initial(ctx):
            bindings.append(ctx.state.raw)
            from data_manager import DataManager
            data = DataManager(ctx.runtime.work_h5.storage)
            ctx.state.set("initial_array", data.read(data.accept(np.full(32768, 3, dtype=np.float32))))
            ctx.state.set("sync", sync)
            ctx.state.set("mode", mode)
            ctx.state.set("accepted_state", 7)
            ctx.state.set("final_path", str(final))

        def block_acceptance(ctx):
            assert sync["received"].wait(20)
            assert not sync["b_started"].is_set(), "Receipt incorrectly enabled readiness"
            sync["crash"].set()
            import psutil
            try:
                psutil.Process(sync["pid"].value).wait(timeout=20)
            except psutil.NoSuchProcess:
                pass
            assert not sync["b_started"].is_set()

        from dataclasses import replace
        descriptors = [PipelineDescriptor("initial", "initial", True,
            pipeline_factory=lambda: SimpleNamespace(name="initial", run=initial))]
        descriptors += [replace(PRODUCER, dag_requires=("initial",)), replace(CONSUMER, dag_requires=("initial",))]
        if mode == "after_receipt":
            descriptors.append(PipelineDescriptor("receipt_barrier", "receipt_barrier", True,
                dag_requires=("initial",), pipeline_factory=lambda: SimpleNamespace(name="receipt_barrier", run=block_acceptance)))
        budget = ResourceBudget(cpu_threads=8, memory_bytes=2 * 1024**3)
        retained = []
        failure = "Worker transport failed" if mode in TRANSPORT_FAULTS else "23"
        with budget.activate(), pytest.raises(PipelineRunError, match=failure):
            run_pipelines_to_output(output_manager=output, pipelines=descriptors,
                pipeline_options={"exchange_producer": (mode,),
                                  "exchange_consumer": ("required_only",) if required_only else ()},
                holodoppler_h5=None, doppler_vision_h5=None, copy_hd_passthrough=False,
                on_finalization_ready=retained.append)
        assert budget.cpu_used == 0
        assert budget.used["transient"] == 0
        pending = retained[0]
        root = pending.model.storage.directory
        assert root.is_dir()
        assert not final.exists()
        if mode in TRANSPORT_FAULTS:
            assert bindings[0]["accepted_state"] == (7 if mode == "missing_receipt" else 17)
            assert ("late_state" in bindings[0]) == (mode == "missing_receipt")
            assert output.path_for(OutputType.PNG, "accepted.png").read_bytes() == b"accepted"
            assert output.path_for(OutputType.PNG, "late.png").exists() == (mode == "missing_receipt")
        if not before_emit:
            assert sync["b_started"].is_set() and sync["b_finished"].is_set()
            assert Path(sync["backing"].value).is_file()
        pending.finish()
        assert not root.exists()
        assert sum(budget.used.values()) == 0
    with h5py.File(final) as result:
        statuses = json.loads(result.attrs["pipeline_status"])["shared"]
        assert statuses["exchange_producer"] == "failed"
        assert statuses["exchange_consumer"] == ("skipped" if required_only else "partial" if before_emit else "succeeded")
        assert "Unaccepted" not in result
        if not required_only:
            assert result["Consumer/independent"][()] == 17
        stages = json.loads(result.attrs["stage_status"])
        if before_emit:
            assert "Accepted/x" not in result
            assert stages["exchange_consumer"]["consume"] == "skipped"
        else:
            np.testing.assert_array_equal(result["Accepted/x"][:], np.arange(ARRAY_SIZE, dtype=np.float32))
            assert result["Consumer/sum"][()] == np.float32(np.arange(ARRAY_SIZE, dtype=np.float64).sum())
            assert stages["exchange_producer"]["produce"] == "failed"
            assert stages["exchange_consumer"]["consume"] == "succeeded"
        if mode == "small_product":
            np.testing.assert_array_equal(result["Accepted/tiny"][:], [3, 7])
        if mode in TRANSPORT_FAULTS:
            assert ("Late/value" in result) == (mode == "missing_receipt")
    label = "[SKIPPED]" if required_only else "[PARTIAL]" if before_emit else "[OK]"
    assert any(label + " exchange_consumer" in message for message in logs)
    assert not any("[OK] exchange_producer" in message for message in logs)
    assert any("Failed exchange_producer" in message and failure in message for message in logs)


@pytest.mark.parametrize("mode", ["duplicate", "conflicting_duplicate", "invalid_reference", "merge_failure",
                                  "artifact_failure", "acceptance_exception", "identity_mismatch", "cancel",
                                  "finish_before_exit", "cleanup_failure"])
def test_real_transport_acceptance_failures_keep_prior_commits_and_release_resources(tmp_path, mode):
    from time import monotonic
    from pipeline_engine.execution.process_executor import ProcessExecutor, ActiveInvocation
    from pipeline_engine.product_exchange import PipelineInvocationId, ProductId, ProductExchangeRegistry, ProductStatus
    from pipeline_engine.product_exchange.acceptance import ProductAcceptance
    from pipeline_engine.io.input.service import InputService
    from pipeline_engine.data.storage import ProductStorage
    from pipeline_engine.execution.protocol import InvocationSettings
    from product_exchange_fixtures import transport_fixture

    mp = multiprocessing.get_context("spawn")
    manager = OutputManager(HoloRunLayout.from_holo(tmp_path / "scan.holo"))
    destination = manager.layout.ef_dir
    destination.mkdir(parents=True)
    if mode == "artifact_failure":
        (destination / "blocked").write_bytes(b"existing file prevents second artifact publication")
    budget = ResourceBudget(cpu_threads=4, memory_bytes=1024**3)
    registry = ProductExchangeRegistry()
    identity = PipelineInvocationId("transport-run", "transport")
    x_id, y_id = ProductId(identity, "x"), ProductId(identity, "y")
    for key in (x_id, y_id):
        registry.declare(key)
    with budget.activate(), ProductStorage(budget=budget) as storage:
        model = OutputModel(storage, manager.path_for(OutputType.H5))
        service = InputService({}, storage, mp)
        state = {"stable": 7}
        acceptance = ProductAcceptance(model=model, state=state, destination=destination, inputs={}, registry=registry)
        executor = ProcessExecutor(context=mp, input_service=service, inputs={}, model=model, manager=manager,
            budget=budget, state={}, settings=InvocationSettings({}, {}, (), (), "both", 1),
            acceptance=acceptance, registry=registry, run_id=identity.run_id)
        acknowledgements = mp.Queue(maxsize=1)
        executor._acks[identity.pipeline] = acknowledgements
        request = budget.default_request(threads=2)
        assert budget.try_acquire(request.threads, request.memory)
        account = budget.memory_account()
        directory = storage.directory / "transport"
        directory.mkdir()
        spec = SimpleNamespace(name=identity.pipeline, invocation_id=identity, owner_id=str(storage.directory.resolve()),
            directory=directory, mode=mode, hold=mp.Event(), final_sent=mp.Event(), acknowledgements=acknowledgements)
        process = mp.Process(target=transport_fixture, args=(spec, executor.events))
        process.start()
        executor.active[identity.pipeline] = ActiveInvocation(process, request.threads, request.memory, account)
        outcomes = []
        try:
            deadline = monotonic() + 25
            while not outcomes:
                assert monotonic() < deadline, "Worker protocol or acceptance deadlocked"
                outcomes.extend(executor.collect(wait=True))
                if mode == "finish_before_exit" and spec.final_sent.is_set() and not spec.hold.is_set():
                    assert process.is_alive() and not outcomes
                    assert budget.cpu_used == request.threads
                    assert budget.used["transient"] == request.memory
                    spec.hold.set()
                if mode == "cancel" and registry.status(x_id) == ProductStatus.AVAILABLE:
                    break
        finally:
            executor.close()
            executor.close()
            service.close()
        assert registry.status(x_id) == ProductStatus.AVAILABLE
        assert registry.status(y_id) == ProductStatus.FAILED
        assert budget.cpu_used == 0
        assert budget.used["transient"] == 0
        if mode in {"duplicate", "finish_before_exit"}:
            assert outcomes[0].error is None
        elif mode == "cleanup_failure":
            assert "during cleanup" in outcomes[0].error
        elif mode != "cancel":
            assert "Coordinator acceptance failed" in outcomes[0].error
        assert not (destination / "first.png").exists()
        assert "Unaccepted" not in model
        assert state == {"stable": 7}
        assert "worker_attribute" not in model.attrs
        assert not list(destination.glob(".ef-*"))
        assert Path(registry.get(x_id).storage.path).is_file()
        model.publish()
        with h5py.File(model.filename) as published:
            np.testing.assert_array_equal(published["Accepted/x"][:], np.arange(32768, dtype=np.float32))
            assert "Unaccepted" not in published
    assert sum(budget.used.values()) == 0
