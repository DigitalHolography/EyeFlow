"""Local input replay, real processes, partial results and final-only publication."""

import json
from dataclasses import replace
from pathlib import Path
from time import monotonic, sleep
from types import SimpleNamespace

import h5py
import numpy as np
import pytest
from runtime_process_fixtures import PROCESS_DESCRIPTOR

from input_output.output_manager import OutputManager, OutputType
from input_output.runs.layout import HoloRunLayout
from pipeline_engine.base import PipelineDescriptor
from pipeline_engine.output_model import OutputModel
from pipeline_engine.resources import ResourceBudget
from pipeline_engine.runtime import PipelineRunError, run_pipelines_to_output
from pipeline_engine.storage import ProductStorage
from utils.logger import Logger


def descriptor(name, run, requires=(), **kwargs):
    return PipelineDescriptor(
        name,
        name,
        True,
        dag_requires=requires,
        pipeline_factory=lambda: SimpleNamespace(name=name, run=run, source_path=__file__),
        **kwargs,
    )


def test_output_journal_preserves_commits_without_reverting_other_paths(tmp_path):
    with ProductStorage(directory=tmp_path / "storage") as storage:
        model = OutputModel(storage, tmp_path / "result.h5")
        with pytest.raises(ValueError), model.checkpoint() as checkpoint:
            model.create_dataset("Committed/value", data=np.array([1, 2]))
            checkpoint.commit(("Committed/value",))
            model.create_dataset("Unfinished/value", data=99)
            raise ValueError("late failure")
        assert "Committed/value" in model
        assert "Unfinished" not in model
        assert not Path(model.filename).exists()
        model.publish()
        with h5py.File(model.filename) as output:
            np.testing.assert_array_equal(output["Committed/value"][:], [1, 2])


def test_spawn_overlaps_independent_work_and_keeps_partial_results(tmp_path):
    source_path = tmp_path / "HD.h5"
    marker = str(tmp_path / "continue")
    values = np.ones((8, 64, 64), dtype=np.float32)
    with h5py.File(source_path, "w") as source:
        source["moment0"] = values
        source.attrs.update(marker=marker, crash=True)
    manager = OutputManager(HoloRunLayout.from_holo(tmp_path / "scan.holo"))
    log_messages = []
    Logger.configure(on_log=log_messages.append)

    def initial_state(ctx):
        ctx.state.set("from_main", 5)

    def independent(ctx):
        deadline = monotonic() + 20
        while not Path(marker + "-started").exists():
            assert monotonic() < deadline, "Process did not start concurrently"
            sleep(0.01)
        assert any("Process fixture source is ready." in text for text in log_messages)
        assert not manager.path_for(OutputType.H5).exists()
        np.testing.assert_array_equal(ctx.inputs.hd.dataset("moment0")[0], values[0])
        ctx.state.set("from_main", 11)
        Path(marker).write_text("continue")
        return {"Main/value": 11}

    def downstream(ctx):
        assert ctx.state.get("from_main") == 11
        return {"Main/downstream": 12}

    def final_consumer(ctx):
        assert ctx.state.get("from_main") == 11
        assert ctx.state.get("from_worker") == 7

    budget = ResourceBudget(cpu_threads=4, memory_bytes=512 * 2**20)
    try:
        with budget.activate(), pytest.raises(PipelineRunError):
            run_pipelines_to_output(
                output_manager=manager,
                pipelines=[
                    descriptor("initial", initial_state),
                    replace(PROCESS_DESCRIPTOR, dag_requires=("initial",)),
                    descriptor("independent", independent, ("initial",)),
                    descriptor("downstream", downstream, ("independent",)),
                    descriptor(
                        "final", final_consumer, ("independent",), wait_for_artifacts=True
                    ),
                ],
                holodoppler_h5=source_path,
                doppler_vision_h5=None,
                copy_hd_passthrough=False,
            )
    finally:
        Logger.reset_current()
    with h5py.File(manager.path_for(OutputType.H5)) as output:
        np.testing.assert_array_equal(output["Worker/video"][:], values + 2)
        assert output["Main/value"][()] == 11
        assert output["Main/downstream"][()] == 12
        statuses = json.loads(output.attrs["pipeline_status"])
        assert statuses["shared"]["process_fixture"] == "failed"
    assert budget.cpu_used == sum(budget.used.values()) == 0

def test_resource_admission_shares_worker_memory_account():
    mib = 2**20
    budget = ResourceBudget(cpu_threads=2, memory_bytes=100 * mib)
    worker = budget.memory_account()
    assert worker.reserve_memory(60 * mib)
    assert not budget.reserve_memory(50 * mib, "transient")
    assert budget.used["retained"] == 60 * mib
    worker.release_memory(60 * mib)
    assert budget.try_acquire(2, 90 * mib)
    assert not budget.try_acquire(1)
    budget.release(2, 90 * mib)
    assert sum(budget.used.values()) == 0


def test_failed_final_serialization_does_not_report_a_missing_output(tmp_path, monkeypatch):
    from test_run_service import _write_input

    from pipeline_engine.run_service import execute_run, resolve_run_spec

    monkeypatch.setenv("EYEFLOW_SCRATCH_DIR", str(tmp_path / "scratch"))
    holo = _write_input(tmp_path)
    spec = resolve_run_spec(
        input_paths=[holo],
        target_names=["producer"],
        pipelines=[descriptor("producer", lambda ctx: {"value": 1})],
    )

    def fail_publish(self):
        raise OSError("destination unavailable")

    monkeypatch.setattr(OutputModel, "publish", fail_publish)
    result = execute_run(spec)
    assert not result.outputs
    assert len(result.failures) == 1
    assert "destination unavailable" in result.failures[0].message


def test_spawn_preparation_failure_keeps_independent_products(tmp_path, monkeypatch):
    from pipeline_engine.scheduler import CoordinatorScheduler

    def fail_launch(self, pipeline, threads, memory):
        raise OSError("cannot prepare worker")

    monkeypatch.setattr(CoordinatorScheduler, "_launch", fail_launch)
    monkeypatch.setenv("EYEFLOW_SCRATCH_DIR", str(tmp_path / "scratch"))
    manager = OutputManager(HoloRunLayout.from_holo(tmp_path / "scan.holo"))
    budget = ResourceBudget(cpu_threads=4, memory_bytes=512 * 2**20)
    with budget.activate(), pytest.raises(PipelineRunError, match="cannot prepare worker"):
        run_pipelines_to_output(
            output_manager=manager,
            pipelines=[PROCESS_DESCRIPTOR, descriptor("safe", lambda ctx: {"safe": 17})],
            holodoppler_h5=None,
            doppler_vision_h5=None,
        )
    with h5py.File(manager.path_for(OutputType.H5)) as output:
        assert output["safe"][()] == 17
    assert budget.cpu_used == sum(budget.used.values()) == 0
