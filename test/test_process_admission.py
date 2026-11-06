"""Coordinator admission must keep worker completion collection responsive."""

import json
import multiprocessing
from pathlib import Path
from types import SimpleNamespace

import h5py
import psutil
import pytest
from runtime_process_fixtures import STAGE_COMPLETION_DESCRIPTOR

from input_output.output_manager import OutputManager, OutputType
from input_output.runs.layout import HoloRunLayout
from pipeline_engine.planning.base import PipelineDescriptor
from pipeline_engine.resources import ResourceBudget
from pipeline_engine.execution.runtime import PipelineRunError, run_pipelines_to_output


def _run_single_cpu_case(destination):
    manager = OutputManager(HoloRunLayout.from_holo(Path(destination) / "scan.holo"))

    def observe_completed_stage(ctx):
        assert ctx.state.get("completed_stage") is True
        return {"Main/value": 11}

    safe = PipelineDescriptor(
        name="safe", description="Independent coordinator product", available=True,
        pipeline_factory=lambda: SimpleNamespace(
            name="safe", source_path=__file__, run=observe_completed_stage,
        ),
    )
    budget = ResourceBudget(cpu_threads=1, memory_bytes=512 * 2**20)
    with budget.activate(), pytest.raises(PipelineRunError, match="23"):
        run_pipelines_to_output(
            output_manager=manager,
            pipelines=[STAGE_COMPLETION_DESCRIPTOR, safe],
            holodoppler_h5=None,
            doppler_vision_h5=None,
            copy_hd_passthrough=False,
        )
    with h5py.File(manager.path_for(OutputType.H5)) as output:
        assert output["Main/value"][()] == 11
        assert output["Stage/value"][()] == 3
        statuses = json.loads(output.attrs["pipeline_status"])["shared"]
        assert statuses == {"stage_completion_fixture": "failed", "safe": "succeeded"}
        stages = json.loads(output.attrs["stage_status"])["stage_completion_fixture"]
        assert stages == {"stage_completion_fixture": "succeeded", "crash": "failed"}
    assert manager.path_for(OutputType.PNG, "completed.png").read_bytes() == b"completed stage"
    assert budget.cpu_used == 0
    assert not any(budget.used.values())


def test_single_cpu_worker_crash_does_not_block_independent_coordinator(tmp_path):
    # Bound this regression in its own process so a blocked admission fails the
    # test instead of hanging pytest. That process spawns the actual worker.
    process = multiprocessing.get_context("spawn").Process(
        target=_run_single_cpu_case, args=(str(tmp_path),),
    )
    process.start()
    try:
        process.join(timeout=25)
        assert not process.is_alive(), "Coordinator blocked instead of collecting worker completion"
        assert process.exitcode == 0
    finally:
        if process.is_alive():
            try:
                children = psutil.Process(process.pid).children(recursive=True)
            except psutil.Error:
                children = []
            for child in children:
                try:
                    child.terminate()
                except psutil.Error:
                    pass
            process.terminate()
            process.join(timeout=5)
            if process.is_alive():
                process.kill()
                process.join(timeout=5)
        process.close()
