"""Automatic writer capture, immutable snapshots, reports and batch overlap."""

import json
from concurrent.futures import ThreadPoolExecutor
from threading import Event
from types import SimpleNamespace

import h5py
import numpy as np
import pytest
from matplotlib import pyplot as plt
from matplotlib.figure import Figure
from matplotlib.ticker import FuncFormatter
from PIL import Image

from input_output.runs.layout import HoloRunLayout
from input_output.output_manager import OutputManager, OutputType
from input_output.writers.article import save_article_figure_pair
from input_output.writers.deferred import save_figure_file
from input_output.writers.eps import EpsArtifactWriter
from input_output.writers.png import FigureArtifactWriter
from pipeline_engine.base import PipelineDescriptor
from pipeline_engine.run_service import execute_run, resolve_run_spec
from pipeline_engine.runtime import PipelineRunError, run_pipelines_to_output
from test_run_service import _write_input
from utils.logger import Logger


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


def _figure():
    fig = Figure(figsize=(3, 2), layout="constrained")
    ax = fig.subplots()
    ax.plot([0, .5, 1], [1, 4, 2], label="signal")
    ax.fill_between([0, .5, 1], [0, 2, 1], alpha=.3)
    ax.set(xlabel="Time (s)", ylabel="Value", title="Example")
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{x:.2f}"))
    ax.legend()
    return fig


def _write_figures(output):
    paths = []
    paths.extend(save_article_figure_pair(output, _figure(), "article", dpi=120))
    paths.append(FigureArtifactWriter(output).save_figure(_figure(), "diagnostic.png", dpi=90))
    paths.append(EpsArtifactWriter(output).save_figure(_figure(), "vector.eps", dpi=90))
    fig = _figure()
    paths.append(save_figure_file(output.path_for(OutputType.PNG, "direct.png"), fig, dpi=90))
    # A queued snapshot must survive closing and later mutation of the source.
    fig.axes[0].lines[0].set_ydata([100, 200, 300])
    plt.close(fig)
    pixels = np.array([[np.nan, -1, 0], [2, 3, 5]], dtype=np.float32)
    paths.append(output.write_png(pixels, "array.png"))
    pixels[:] = 200
    return paths


def _eps_payload(path):
    return [line for line in path.read_bytes().splitlines()
            if not line.startswith((b"%%CreationDate:", b"%%Title:"))]


def test_shared_writers_capture_automatically_and_report_consumes_identical_files(tmp_path):
    original = OutputManager(HoloRunLayout.from_holo(tmp_path / "sync/scan.holo"))
    managed = OutputManager(HoloRunLayout.from_holo(tmp_path / "async/scan.holo"))
    expected = _write_figures(original)
    written = []

    def source(ctx):
        written.extend(_write_figures(ctx.output))
        assert all(not path.exists() for path in written)
        return {"Numerical/value": 42}

    def report(ctx):
        assert all(path.is_file() for path in written)
        # Report-owned intermediate figures remain synchronous for immediate use.
        temporary = ctx.output.write_png(np.ones((2, 2)), "report_intermediate.png")
        assert temporary.is_file()
        return {"Report/value": 1}

    path = _run(managed, [
        _descriptor("producer", source),
        _descriptor("report", report, wait_for_artifacts=True),
    ])
    for before, after in zip(expected, written):
        if before.suffix == ".png":
            with Image.open(before) as left, Image.open(after) as right:
                np.testing.assert_array_equal(np.asarray(left), np.asarray(right))
        else:
            assert _eps_payload(before) == _eps_payload(after)
    with h5py.File(path) as output:
        assert output["Numerical/value"][()] == 42
        assert output["Report/value"][()] == 1
        assert json.loads(output.attrs["artifact_status"]) == {"producer:figures": "succeeded"}


def test_unserializable_figure_is_an_export_failure_and_retains_numerics(tmp_path, monkeypatch):
    # Serialization fails before rendering, so these fixture workers never run Matplotlib.
    monkeypatch.setattr(
        "pipeline_engine.background_artifacts.ProcessPoolExecutor",
        lambda **kwargs: ThreadPoolExecutor(max_workers=1),
    )
    manager = OutputManager(HoloRunLayout.from_holo(tmp_path / "scan.holo"))

    def source(ctx):
        dataset = ctx.output.h5.file.create_dataset("Numerical/value", data=42)
        fig = _figure()
        fig.axes[0].yaxis.set_major_formatter(FuncFormatter(lambda x, pos: str(dataset[()])))
        FigureArtifactWriter(ctx.output).save_figure(fig, "unsupported.png")

    with pytest.raises(PipelineRunError, match="Cannot snapshot") as error:
        _run(manager, [_descriptor("producer", source)])
    with h5py.File(error.value.output_h5_path) as output:
        assert output["Numerical/value"][()] == 42
        assert json.loads(output.attrs["pipeline_status"])["shared"]["producer"] == "succeeded"
        assert json.loads(output.attrs["artifact_status"]) == {"producer:figures": "failed"}


def test_queued_writers_are_cancelled_when_numerical_persistence_fails(tmp_path, monkeypatch):
    def forbidden_pool(**kwargs):
        pytest.fail("A rolled-back producer must not start an exporter")

    monkeypatch.setattr("pipeline_engine.background_artifacts.ProcessPoolExecutor", forbidden_pool)
    manager = OutputManager(HoloRunLayout.from_holo(tmp_path / "scan.holo"))

    def source(ctx):
        _write_figures(ctx.output)
        return {"Invalid/value": object()}

    with pytest.raises(PipelineRunError):
        _run(manager, [_descriptor("producer", source)])
    assert not list(manager.layout.ef_dir.rglob("*.png"))
    assert not list(manager.layout.ef_dir.rglob("*.eps"))


def test_next_file_starts_before_exports_finish_and_completion_logs_are_asynchronous(tmp_path, monkeypatch):
    pools = []

    def fixture_pool(**kwargs):
        pools.append(kwargs)
        return ThreadPoolExecutor(max_workers=2)

    monkeypatch.setattr("pipeline_engine.background_artifacts.ProcessPoolExecutor", fixture_pool)
    first = _write_input(tmp_path / "inputs", "first")
    second = _write_input(tmp_path / "inputs", "second")
    render_started, second_started, release, first_complete = (Event() for _ in range(4))
    messages = []

    def on_log(message):
        messages.append(message)
        if "Completed run for first.holo:" in message:
            first_complete.set()

    def export(snapshot, manager):
        if manager.layout.stem == "first":
            render_started.set()
            assert release.wait(10)
        path = manager.path_for(OutputType.PNG, "fixture.png")
        path.parent.mkdir(parents=True)
        path.write_bytes(snapshot.read_bytes())
        return [path]

    def source(ctx):
        name = ctx.output.manager.layout.stem
        if name == "second":
            second_started.set()
            assert render_started.wait(10)
            assert not first_complete.is_set()
            first_path = spec.requests[0].output_manager.path_for(OutputType.H5)
            assert not first_path.exists(), "HDF5 publication waits for required exports"
            release.set()
            # The first file's final completion log must arrive while this
            # second descriptor is still executing, not at a later boundary.
            assert first_complete.wait(10)
        ctx.defer_artifacts("fixtures", export, lambda path: path.write_bytes(name.encode()))
        return {"Numerical/value": 42}

    spec = resolve_run_spec(
        input_paths=[first, second], target_names=["producer"],
        pipelines=[_descriptor("producer", source)],
    )
    previous = Logger._current
    Logger.configure(on_log=on_log)
    try:
        result = execute_run(spec)
    finally:
        Logger._current = previous
    assert result.succeeded, result.failures
    assert second_started.is_set()
    assert len(pools) == 1, "The same renderer pool should serve both inputs"
    assert len(result.outputs) == 2
    assert any("calculations complete; artifact exports continue" in message for message in messages)
    export_done = next(i for i, message in enumerate(messages) if "[EXPORT COMPLETE] first" in message)
    file_done = next(i for i, message in enumerate(messages) if "Completed run for first.holo:" in message)
    assert export_done < file_done
    for path in result.outputs:
        with h5py.File(path) as output:
            assert set(json.loads(output.attrs["artifact_status"]).values()) == {"succeeded"}
