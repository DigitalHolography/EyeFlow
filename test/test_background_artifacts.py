"""Background exports overlap calculations and respect checkpoints/completion."""

import json
from concurrent.futures import ThreadPoolExecutor
from functools import partial
from threading import Event
from types import SimpleNamespace
from unittest.mock import patch

import h5py
import numpy as np
import pytest
from PIL import Image

from input_output.writers.h5_output import PipelineH5Output
from input_output.runs.layout import HoloRunLayout
from input_output.output_manager import OutputManager, OutputType
from input_output.schema import EyeFlowOutputPaths
from pipeline_engine.background_artifacts import ArtifactExportResult, BackgroundArtifacts
from pipeline_engine.base import PipelineDescriptor, pipeline, registerPipeline
from pipeline_engine.context import PipelineOutput
from pipeline_engine.runtime import PipelineRunError, run_pipelines_to_output
from pipelines.low_freq_spectra.artifacts import export_spectrogram_figures
from pipelines.low_freq_spectra.figure_export import (
    export_spectral_snapshot, snapshot_spectral_figures,
)
from utils.logger import Logger


@pytest.fixture
def threaded_exporter(monkeypatch):
    # These jobs write fixture bytes only. The real renderer uses a process;
    # test its spawn path separately, without running Matplotlib in threads.
    monkeypatch.setattr(
        "pipeline_engine.background_artifacts.ProcessPoolExecutor",
        lambda **kwargs: ThreadPoolExecutor(max_workers=1),
    )


def _manager(tmp_path):
    return OutputManager(HoloRunLayout.from_holo(tmp_path / "scan.holo"))


def _descriptor(name, run, **kwargs):
    return PipelineDescriptor(
        name, name, True,
        pipeline_factory=lambda: SimpleNamespace(name=name, run=run, source_path=__file__),
        **kwargs,
    )


def _run(manager, descriptors):
    return run_pipelines_to_output(
        output_manager=manager, pipelines=descriptors,
        holodoppler_h5=None, doppler_vision_h5=None,
    )


@pytest.mark.parametrize("failing_sibling", [False, True])
def test_export_overlaps_analysis_and_survives_unrelated_rollback(
    tmp_path, threaded_exporter, failing_sibling,
):
    manager = _manager(tmp_path)
    started, release, finished = Event(), Event(), Event()
    target = manager.path_for(OutputType.PNG, "background.png")

    def export(snapshot, staged):
        started.set()
        assert release.wait(10), "The next pipeline never ran while export was pending"
        path = staged.path_for(OutputType.PNG, "background.png")
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(snapshot.read_bytes())
        finished.set()
        return [path]

    def source(ctx):
        ctx.output.h5.write("Numerical/value", 42)
        ctx.defer_artifacts("figures", export, lambda path: path.write_bytes(b"snapshot"))
        assert not started.is_set(), "Exporter started before the checkpoint committed"

    def sibling(ctx):
        assert started.wait(10)
        assert not target.exists()
        release.set()
        assert finished.wait(10)
        # Publication and its completion log can now happen during a later
        # descriptor. Rollback must recognize those background-owned paths.
        ctx._background_artifacts.collect(wait=True)
        assert target.read_bytes() == b"snapshot"
        if failing_sibling:
            ctx.output.h5.write("Partial/value", -1)
            raise ValueError("injected unrelated failure")
        return {"Sibling/value": 7}

    def report(ctx):
        assert target.read_bytes() == b"snapshot"
        assert json.loads(ctx.output.h5.file.attrs["artifact_status"])["source:figures"] == "succeeded"
        return {"Report/value": 1}

    descriptors = [
        _descriptor("source", source), _descriptor("sibling", sibling),
        _descriptor("report", report, wait_for_artifacts=True),
    ]
    if failing_sibling:
        with pytest.raises(PipelineRunError, match="injected unrelated failure") as error:
            _run(manager, descriptors)
        path = error.value.output_h5_path
    else:
        path = _run(manager, descriptors)
    assert target.read_bytes() == b"snapshot"
    with h5py.File(path) as output:
        assert output["Numerical/value"][()] == 42
        assert output["Report/value"][()] == 1
        assert "Partial" not in output
        assert json.loads(output.attrs["artifact_status"]) == {"source:figures": "succeeded"}
        assert json.loads(output.attrs["artifact_failures"]) == {}


def test_completion_joins_export_without_an_artifact_consumer(tmp_path, threaded_exporter):
    manager = _manager(tmp_path)

    def export(snapshot, staged):
        path = staged.path_for(OutputType.EPS, "final.eps")
        path.parent.mkdir(parents=True)
        path.write_bytes(snapshot.read_bytes())
        return [path]

    def source(ctx):
        ctx.defer_artifacts("figures", export, lambda path: path.write_bytes(b"EPS fixture"))

    path = _run(manager, [_descriptor("source", source)])
    assert manager.path_for(OutputType.EPS, "final.eps").read_bytes() == b"EPS fixture"
    with h5py.File(path) as output:
        assert json.loads(output.attrs["artifact_status"])["source:figures"] == "succeeded"


@pytest.mark.parametrize("failure_stage", ["render", "submit", "render_result"])
def test_export_failure_preserves_numerics_continues_analysis_and_skips_report(
    tmp_path, threaded_exporter, monkeypatch, failure_stage,
):
    manager = _manager(tmp_path)

    def fail(*args, **kwargs):
        raise ValueError("injected export failure")

    if failure_stage == "submit":
        monkeypatch.setattr("pipeline_engine.background_artifacts.ProcessPoolExecutor", fail)

    def export(*args):
        if failure_stage == "render_result":
            return ArtifactExportResult([], ["[WARNING] fixture warning"], "injected export failure")
        return fail()

    def source(ctx):
        ctx.output.h5.write("Numerical/value", 42)
        ctx.defer_artifacts("figures", export, lambda path: path.write_bytes(b"snapshot"))

    def report(ctx):
        pytest.fail("An artifact consumer must not run after an export failure")

    with pytest.raises(PipelineRunError, match="injected export failure") as error:
        _run(manager, [
            _descriptor("source", source),
            _descriptor("sibling", lambda ctx: {"Sibling/value": 7}),
            _descriptor("report", report, wait_for_artifacts=True),
        ])
    with h5py.File(error.value.output_h5_path) as output:
        assert output["Numerical/value"][()] == 42
        assert output["Sibling/value"][()] == 7
        statuses = json.loads(output.attrs["pipeline_status"])["shared"]
        assert statuses == {"source": "succeeded", "sibling": "succeeded", "report": "skipped"}
        assert json.loads(output.attrs["artifact_status"]) == {"source:figures": "failed"}
        assert "injected export failure" in output.attrs["artifact_failures"]


def test_failed_producer_discards_snapshot_without_starting_export(tmp_path, threaded_exporter):
    manager = _manager(tmp_path)
    snapshots = []
    called = Event()

    def snapshotter(path):
        snapshots.append(path)
        path.write_bytes(b"snapshot")

    def source(ctx):
        ctx.defer_artifacts("figures", lambda *args: called.set(), snapshotter)
        raise ValueError("producer failed")

    with pytest.raises(PipelineRunError, match="producer failed"):
        _run(manager, [_descriptor("source", source)])
    assert not called.is_set()
    assert snapshots and all(not path.parent.exists() for path in snapshots)


@pytest.mark.parametrize("invalid_path", ["outside", "collision"])
def test_failed_publication_leaves_no_partial_files(tmp_path, threaded_exporter, invalid_path):
    manager = _manager(tmp_path)
    target = manager.path_for(OutputType.PNG, "second.png")
    outside = tmp_path / "outside.png"
    outside.write_bytes(b"outside")
    if invalid_path == "collision":
        target.parent.mkdir(parents=True)
        target.write_bytes(b"original")

    def export(snapshot, staged):
        first = staged.path_for(OutputType.PNG, "first.png")
        first.parent.mkdir(parents=True)
        first.write_bytes(b"first")
        return [first, outside if invalid_path == "outside" else first.with_name(target.name)]

    with BackgroundArtifacts() as jobs:
        job = jobs.prepare("test", manager, export, lambda path: path.write_bytes(b"fixture"))
        jobs.start(job)
        jobs.collect(wait=True)
        assert jobs.statuses == {"test": "failed"}
        assert jobs.failures
    assert not manager.path_for(OutputType.PNG, "first.png").exists()
    assert outside.read_bytes() == b"outside"
    if invalid_path == "collision":
        assert target.read_bytes() == b"original"
    assert not list(manager.layout.ef_dir.rglob("*.partial"))


def test_real_spawned_renderer_uses_independent_reduced_snapshot(tmp_path):
    schema = EyeFlowOutputPaths.active()
    manager = _manager(tmp_path / "background")
    synchronous = _manager(tmp_path / "synchronous")
    mask = np.array([[1, 0, 0], [0, 1, 1]], dtype=np.uint8)
    with h5py.File(tmp_path / "primary.h5", "w") as source:
        source.create_dataset(schema.segmentation.gaps, data=mask)
        # Missing calibration warns in the child, while the mask still exports.
        source.create_dataset(schema.spectrograms_masked, data=np.ones((2, 4)))
        source.create_dataset("UnrelatedLargeVideo", shape=(10, 64, 128, 128), dtype="f4")
        output = PipelineOutput(manager, PipelineH5Output(source))
        expected = export_spectrogram_figures(
            PipelineOutput(synchronous, PipelineH5Output(source)),
        )
        with BackgroundArtifacts() as jobs:
            job = jobs.prepare(
                "spectral", manager, export_spectral_snapshot,
                partial(snapshot_spectral_figures, output),
            )
            with h5py.File(job.snapshot) as snapshot:
                assert "UnrelatedLargeVideo" not in snapshot
                np.testing.assert_array_equal(snapshot[schema.segmentation.gaps][:], mask)
            # Neither subsequent writes nor a closed primary handle affect the worker.
            source[schema.segmentation.gaps][:] = 0
            source.close()
            jobs.start(job)
            with patch.object(Logger, "log") as log:
                jobs.collect(wait=True)
            assert jobs.statuses == {"spectral": "succeeded"}, jobs.failures
            assert any("FrequencyBinTable" in str(call) for call in log.call_args_list)
    for path in expected:
        actual = manager.layout.ef_dir / path.relative_to(synchronous.layout.ef_dir)
        assert actual.is_file()
        if path.suffix == ".png":
            np.testing.assert_array_equal(np.asarray(Image.open(actual)), np.asarray(Image.open(path)))
    assert not job.snapshot.parent.exists()


def test_decorators_preserve_artifact_consumer_metadata(monkeypatch):
    # Replace the registry locally so discovery state is untouched.
    registry = {}
    monkeypatch.setattr("pipeline_engine.base.PIPELINE_REGISTRY", registry)

    @pipeline("function_consumer", wait_for_artifacts=True)
    def run(ctx):
        pass

    @registerPipeline("class_consumer", wait_for_artifacts=True)
    class Consumer:
        pass

    for descriptor in registry.values():
        assert descriptor.wait_for_artifacts
        assert descriptor.instantiate().wait_for_artifacts
