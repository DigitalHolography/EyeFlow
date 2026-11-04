"""Result acceptance must preserve earlier work on collision or I/O failure."""

from pathlib import Path

import numpy as np
import pytest

from pipeline_engine.output_model import OutputModel, Record
from pipeline_engine.product_exchange.acceptance import accept_result
from pipeline_engine.protocol import InvocationResult
from pipeline_engine.storage import ProductStorage


def _result(tmp_path):
    root = tmp_path / "worker"
    root.mkdir()
    paths = tuple(root / name for name in ("first.png", "second.png"))
    for path in paths:
        path.write_bytes(path.name.encode())
    return InvocationResult(
        records={"": Record(attrs={"worker_attribute": 1}),
                 "Worker": Record(), "Worker/value": Record("dataset", np.array([7]))},
        state={"worker_path": paths[0], "updated": 10},
        artifacts=paths, artifact_root=root,
    )


def _accept(result, model, state, destination):
    accept_result(result, model=model, state=state, destination=destination,
                  name="worker", stages={"calculate": "succeeded"})


@pytest.mark.parametrize("collision", ["artifact", "dataset"])
def test_collision_leaves_prior_output_and_state_untouched(tmp_path, collision):
    result = _result(tmp_path)
    destination = tmp_path / "published"
    destination.mkdir()
    state = {"updated": 1}
    with ProductStorage() as storage:
        model = OutputModel(storage, destination / "future.h5")
        model.attrs["earlier"] = 2
        if collision == "artifact":
            (destination / "second.png").write_bytes(b"earlier")
        else:
            model.create_dataset("Worker/value", data=[3])
        with pytest.raises((ValueError, FileExistsError)):
            _accept(result, model, state, destination)
        assert state == {"updated": 1}
        assert dict(model.attrs) == {"earlier": 2}
        assert not (destination / "first.png").exists()
        if collision == "dataset":
            np.testing.assert_array_equal(model["Worker/value"][:], [3])
        else:
            assert "Worker" not in model
            assert (destination / "second.png").read_bytes() == b"earlier"


def test_copy_failure_rolls_back_model_state_and_files(tmp_path, monkeypatch):
    import pipeline_engine.product_exchange.acceptance as execution

    result = _result(tmp_path)
    destination = tmp_path / "published"
    state = {"updated": 1}
    copyfile = execution.shutil.copyfile

    def fail_second(source, target):
        if Path(source).name == "second.png":
            Path(target).write_bytes(b"partial")
            raise OSError("destination disconnected")
        return copyfile(source, target)

    monkeypatch.setattr(execution.shutil, "copyfile", fail_second)
    with ProductStorage() as storage:
        model = OutputModel(storage, destination / "future.h5")
        model.attrs["earlier"] = 2
        with pytest.raises(OSError, match="disconnected"):
            _accept(result, model, state, destination)
        assert state == {"updated": 1}
        assert dict(model.attrs) == {"earlier": 2}
        assert "Worker" not in model
        assert not list(destination.iterdir())


def test_log_callback_failure_does_not_invalidate_committed_records(tmp_path, monkeypatch):
    from dataclasses import replace
    result = replace(_result(tmp_path), messages=("presentation callback",))
    def fail_log(message):
        raise RuntimeError("presentation failure")
    monkeypatch.setattr("pipeline_engine.product_exchange.acceptance.Logger.log", fail_log)
    with ProductStorage() as storage:
        model = OutputModel(storage, tmp_path / "published/final.h5")
        state = {}
        _accept(result, model, state, tmp_path / "published")
        np.testing.assert_array_equal(model["Worker/value"][:], [7])
        assert state["updated"] == 10
