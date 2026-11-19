"""Fast, isolated filesystem defaults shared by the test suite."""

from __future__ import annotations

import os
from pathlib import Path

import h5py
import pytest

from input_output.writers.h5_publication import (
    CompletedLocalOutput,
    reconcile_attributes,
    write_record,
)
from pipeline_engine.io.manager import OutputService as _SpawnedOutputService


@pytest.fixture(scope="session", autouse=True)
def _local_scratch_directory(tmp_path_factory):
    previous = os.environ.get("EYEFLOW_SCRATCH_DIR")
    os.environ["EYEFLOW_SCRATCH_DIR"] = str(tmp_path_factory.mktemp("scratch"))
    yield
    if previous is None:
        os.environ.pop("EYEFLOW_SCRATCH_DIR", None)
    else:
        os.environ["EYEFLOW_SCRATCH_DIR"] = previous


class _InlineOutputService:
    """Materialize test outputs without testing the spawned writer again."""

    def __init__(self, model, *, budget=None, **_kwargs):
        self.model = model
        self.budget = budget if budget is not None else model.storage.budget
        self.metrics = model.storage.metrics
        self.local_path = model.storage.path(".h5")
        self._finished = self._closed = False
        self._completed = None
        model.bind_commit_listener(self.submit)

    def submit(self, _model):
        return None

    def check_health(self):
        if self._closed and not self._finished:
            raise OSError("Output serialization was aborted")

    def drain(self):
        self.check_health()

    def finish(self, snapshot=None):
        if self._finished:
            return self._completed
        self.check_health()
        snapshot = snapshot or self.model.finalization_snapshot()
        records = self.model.committed_records()
        with h5py.File(self.local_path, "w") as output:
            for path in snapshot.dataset_paths:
                write_record(output, path, records[path])
            for path, record in sorted(
                snapshot.groups,
                key=lambda item: (item[0].count("/"), item[0]),
            ):
                group = output.require_group(path) if path else output
                reconcile_attributes(group, record)
        self._completed = CompletedLocalOutput(
            self.local_path, self.local_path.stat().st_size
        )
        self._finished = True
        return self._completed

    def request_abort(self):
        return None

    def close(self):
        if self._closed:
            return
        self._closed = True
        self.model.detach_commit_listener(self.submit)
        Path(self.local_path).unlink(missing_ok=True)


@pytest.fixture
def inline_output_service(monkeypatch):
    """Use when a test's contract is above the output-process boundary."""

    monkeypatch.setattr(
        "pipeline_engine.io.manager.OutputService", _InlineOutputService
    )


@pytest.fixture
def real_output_service(monkeypatch, inline_output_service):
    """Restore the writer process for tests that own that boundary."""

    monkeypatch.setattr(
        "pipeline_engine.io.manager.OutputService", _SpawnedOutputService
    )
