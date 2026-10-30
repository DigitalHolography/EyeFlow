"""Local input replay, real processes, partial results and final-only publication."""

import json
import multiprocessing
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
from pipeline_engine.input_service import InputService
from pipeline_engine.output_model import OutputModel
from pipeline_engine.resources import ResourceBudget
from pipeline_engine.runtime import PipelineRunError, run_pipelines_to_output
from pipeline_engine.storage import ArrayReference, ProductStorage, is_local_path


def descriptor(name, run, requires=(), **kwargs):
    return PipelineDescriptor(
        name,
        name,
        True,
        dag_requires=requires,
        pipeline_factory=lambda: SimpleNamespace(name=name, run=run, source_path=__file__),
        **kwargs,
    )


def test_source_is_staged_once_and_replayed_read_only(tmp_path):
    source_path = tmp_path / "source.h5"
    original = np.arange(128 * 256, dtype=np.float32).reshape(128, 256)
    with h5py.File(source_path, "w") as source:
        source["video"] = original
        source["video_alias"] = source["video"]
    with h5py.File(source_path) as source, ProductStorage() as storage:
        service = InputService({"hd": source}, storage, multiprocessing.get_context("spawn"))
        first, second = service.client(), service.client()
        try:
            a, b = service.files(first)["hd"], service.files(second)["hd"]
            np.testing.assert_array_equal(a["video"][1:3], original[1:3])
            np.testing.assert_array_equal(b["video"][:], original)
            np.testing.assert_array_equal(b["video_alias"][:], original)
            assert b["video"].id == b["video_alias"].id
            assert service.source_reads == {("hd", "video"): 1}
            assert isinstance(service.cache[("hd", "video")], ArrayReference)
            assert not b["video"][:].flags.writeable
            assert is_local_path(storage.directory)
            assert all(path.suffix != ".tmp" for path in storage.directory.iterdir())
        finally:
            first.close()
            second.close()
            service.close()


def test_output_journal_preserves_commits_without_reverting_other_paths(tmp_path):
    with ProductStorage() as storage:
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


def test_reports_read_the_assembled_view_before_publication(tmp_path):
    from input_output.reports.pdf_data import extract_parameters_from_h5
    from input_output.schema import EyeFlowOutputPaths
    from input_output.writers.h5_output import PipelineH5Output

    with ProductStorage() as storage:
        model = OutputModel(storage, tmp_path / "future.h5")
        output = PipelineH5Output(model, processing_root="ProcessingMoments")
        output.write(
            EyeFlowOutputPaths.active().artery_per_beat.velocity_signal,
            np.array([1.0, 3.0], dtype=np.float32),
            unit="mm/s",
        )
        parameters = extract_parameters_from_h5(
            Path(model.filename), processing_root="ProcessingMoments", output_source=model
        )
        assert parameters["Average_Arterial_Velocity"] == {"value": 2.0, "unit": "mm/s"}
        assert not Path(model.filename).exists()


def test_final_serialization_preserves_sparse_and_same_file_virtual_storage(tmp_path):
    from input_output.schema.output_payload import DatasetValue
    from input_output.writers.h5_output import PipelineH5Output

    unique = np.arange(9 * 2048, dtype=np.float32).reshape(9, 2048)
    unique[2] = np.nan
    values = np.concatenate([unique, unique[-2:0:-1]])
    with ProductStorage() as storage:
        model = OutputModel(storage, tmp_path / "final.h5")
        PipelineH5Output(model).write_many(
            {
                "FFT": DatasetValue(
                    values,
                    {"unit": "mm/s"},
                    {
                        "mirror_fft_axis": 0,
                        "sparse_nan_chunks": True,
                        "chunks": (1, 2048),
                        "compression": "gzip",
                    },
                )
            }
        )
        model.publish()
        with h5py.File(model.filename) as output:
            assert output["FFT"].is_virtual
            assert all(source.file_name == "." for source in output["FFT"].virtual_sources())
            np.testing.assert_array_equal(output["FFT"][:], values)
            assert output["FFT"].attrs["unit"] == "mm/s"


@pytest.mark.parametrize("crash", [False, True])
def test_spawn_overlaps_independent_work_and_keeps_partial_results(tmp_path, crash):
    source_path = tmp_path / "HD.h5"
    marker = str(tmp_path / "continue")
    values = np.ones((8, 64, 64), dtype=np.float32)
    with h5py.File(source_path, "w") as source:
        source["moment0"] = values
        source.attrs["marker"] = marker
        source.attrs["crash"] = crash
    manager = OutputManager(HoloRunLayout.from_holo(tmp_path / "scan.holo"))

    def initial_state(ctx):
        ctx.state.set("from_main", 5)

    def independent(ctx):
        deadline = monotonic() + 20
        while not Path(marker + "-started").exists():
            assert monotonic() < deadline, "Process did not start concurrently"
            sleep(0.01)
        assert not manager.path_for(OutputType.H5).exists()
        # Same source, different slice: the input service must reuse local data.
        np.testing.assert_array_equal(ctx.inputs.hd.dataset("moment0")[0], values[0])
        ctx.state.set("from_main", 11)
        Path(marker).write_text("continue")
        return {"Main/value": 11}

    def downstream(ctx):
        assert ctx.state.get("from_main") == 11
        return {"Main/downstream": 12}

    def final_consumer(ctx):
        assert ctx.state.get("from_main") == 11, "Worker overwrote newer coordinator state"
        assert ctx.state.get("from_worker") == 7

    budget = ResourceBudget(cpu_threads=4, memory_bytes=512 * 2**20)
    with budget.activate(), pytest.raises(PipelineRunError):
        run_pipelines_to_output(
            output_manager=manager,
            pipelines=[
                descriptor("initial", initial_state),
                replace(PROCESS_DESCRIPTOR, dag_requires=("initial",)),
                descriptor("independent", independent, ("initial",)),
                descriptor("downstream", downstream, ("independent",)),
                descriptor("final", final_consumer, ("independent",), wait_for_artifacts=True),
            ],
            holodoppler_h5=source_path,
            doppler_vision_h5=None,
            copy_hd_passthrough=False,
        )
    with h5py.File(manager.path_for(OutputType.H5)) as output:
        np.testing.assert_array_equal(output["Worker/video"][:], values + 2)
        assert output["Main/value"][()] == 11
        assert output["Main/downstream"][()] == 12
        statuses = json.loads(output.attrs["pipeline_status"])
        assert statuses["shared"]["process_fixture"] == "failed"
        if not crash:
            assert output["Worker/independent"][0] == 17
            assert (
                manager.path_for(OutputType.PNG, "worker.png").read_bytes() == b"fixture artifact"
            )
    assert budget.cpu_used == 0
    assert sum(budget.used.values()) == 0


def test_resource_admission_shares_worker_memory_account():
    budget = ResourceBudget(cpu_threads=2, memory_bytes=100)
    worker = budget.memory_account()
    assert worker.reserve_memory(60)
    assert not budget.reserve_memory(50, "transient")
    assert budget.used["retained"] == 60
    worker.release_memory(60)
    assert budget.try_acquire(2, 90)
    assert not budget.try_acquire(1)
    budget.release(2, 90)
    assert sum(budget.used.values()) == 0


def test_coordinator_uses_remaining_cpu_slots_instead_of_waiting_for_all():
    from pipeline_engine.resources import thread_allowance

    budget = ResourceBudget(cpu_threads=4, memory_bytes=1024)
    assert budget.try_acquire(3)
    with budget.compute(4):
        assert thread_allowance() == 1
        assert budget.cpu_used == 4
    assert budget.cpu_used == 3
    budget.release(3)


def test_failed_final_serialization_does_not_report_a_missing_output(tmp_path, monkeypatch):
    from test_run_service import _write_input

    from pipeline_engine.run_service import execute_run, resolve_run_spec

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


def test_native_initialization_failure_releases_compute_reservations(monkeypatch):
    budget = ResourceBudget(cpu_threads=2, memory_bytes=1024)

    def fail_native(count):
        raise RuntimeError("native initialization failed")

    monkeypatch.setattr("pipeline_engine.resources.configure_native_threads", fail_native)
    with pytest.raises(RuntimeError, match="native initialization failed"), budget.compute(1, 256):
        pass
    assert budget.cpu_used == 0
    assert sum(budget.used.values()) == 0


def test_spawn_preparation_failure_keeps_independent_products(tmp_path, monkeypatch):
    from pipeline_engine.scheduler import CoordinatorScheduler

    def fail_launch(self, pipeline, threads, memory):
        raise OSError("cannot prepare worker")

    monkeypatch.setattr(CoordinatorScheduler, "_launch", fail_launch)
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
    assert budget.cpu_used == 0
    assert sum(budget.used.values()) == 0


@pytest.mark.parametrize("missing_artery", [False, True])
def test_production_displacement_uses_local_process_inputs_and_publishes_videos(
    tmp_path, missing_artery
):
    from contextlib import nullcontext

    from input_output.schema.doppler_view import DV_ARTERY_MASK_PATH, DV_VEIN_MASK_PATH
    from pipelines import load_pipeline_catalog

    pipelines, _ = load_pipeline_catalog()
    displacement = next(p for p in pipelines if p.name == "displacement_map")
    hd_path, dv_path = tmp_path / "HD.h5", tmp_path / "DV.h5"
    y, x = np.mgrid[:32, :32]
    image = np.exp(-((x - 16) ** 2 + (y - 16) ** 2) / 30).astype(np.float32)
    with h5py.File(hd_path, "w") as hd:
        hd["moment0"] = np.stack([image, np.roll(image, 1, axis=1)])
    with h5py.File(dv_path, "w") as dv:
        mask = np.zeros((32, 32), dtype=np.uint8)
        mask[12:20, 12:20] = 1
        if not missing_artery:
            dv[DV_ARTERY_MASK_PATH] = mask
        dv[DV_VEIN_MASK_PATH] = mask
    manager = OutputManager(HoloRunLayout.from_holo(tmp_path / "scan.holo"))
    with (
        ResourceBudget(cpu_threads=4, memory_bytes=512 * 2**20).activate(),
        pytest.raises(PipelineRunError) if missing_artery else nullcontext(),
    ):
        run_pipelines_to_output(
            output_manager=manager,
            pipelines=[displacement],
            holodoppler_h5=hd_path,
            doppler_vision_h5=dv_path,
            copy_hd_passthrough=False,
        )
    videos = list(manager.dir_for(OutputType.AVI).rglob("*.avi"))
    assert len(videos) == (1 if missing_artery else 2)
    assert all(video.stat().st_size > 0 for video in videos)
    with h5py.File(manager.path_for(OutputType.H5)) as output:
        statuses = json.loads(output.attrs["stage_status"])["displacement_map"]
        assert statuses == {
            "prepare": "succeeded",
            "artery": "failed" if missing_artery else "succeeded",
            "vein": "succeeded",
        }
