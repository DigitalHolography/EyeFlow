"""Acquisition I/O ownership, with independent input and retained output lifetimes."""

from contextlib import ExitStack
import json
import multiprocessing
from pathlib import Path
from time import perf_counter
from threading import Event, RLock

from input_output.writers.h5 import open_h5, write_root_attributes
from input_output.writers.h5_publication import OutputPublisher
from pipeline_engine.execution.background_artifacts import BackgroundArtifacts
from pipeline_engine.output.output_model import ProductCheckpoint
from utils.logger import Logger

from .input.service import InputService
from .output.service import OutputService


class IOManager:
    """Own services and storage; the existing runtime still schedules science.

    The scientific ExitStack closes only input resources. This manager is then
    retained by batch finalization through artifact completion and publication.
    """

    def __init__(self, model, *, artifacts=None, owns_artifacts=False, budget=None,
                 context=None, output_policy=None):
        self.model = model
        self.artifacts = artifacts if artifacts is not None else BackgroundArtifacts()
        self.owns_artifacts = owns_artifacts or artifacts is None
        self.context = context or multiprocessing.get_context("spawn")
        self.budget = budget if budget is not None else model.storage.budget
        self._inputs = ExitStack()
        self._sources = None
        self.input = None
        self.files = {}
        self._input_closed = self._closed = self._finished = False
        self._lifecycle = RLock()
        self._cancelled = Event()
        try:
            self.output = OutputService(model, budget=self.budget, context=self.context, policy=output_policy)
        except BaseException:
            model.storage.close()
            raise

    def open_sources(self, hd_path, dv_path):
        if self._input_closed or self._sources is not None:
            raise RuntimeError("Acquisition inputs have already been opened or closed")
        with Logger.perf_stage("run/input_hdf5_open"):
            self._sources = {
                slot: self._inputs.enter_context(open_h5(path, "r")) if path is not None else None
                for slot, path in (("hd", hd_path), ("dv", dv_path))
            }
        return self._sources["hd"], self._sources["dv"]

    def start_input(self, policy=None):
        if self._input_closed or self._sources is None or self.input is not None:
            raise RuntimeError("Input service cannot be started in this acquisition state")
        self.input = InputService(self._sources, self.model.storage, self.context, policy=policy)
        self._inputs.callback(self.input.close)
        client = self.input.client()
        self._inputs.callback(client.close)
        self.files = self.input.files(client)
        return self.input

    def close_input(self):
        if not self._input_closed:
            self._input_closed = True
            self._inputs.close()

    def check_health(self):
        self.output.check_health()
        self.output.submit(self.model)
        self.output.check_health()

    def finish(self):
        with self._lifecycle:
            return self._finish()

    @property
    def closed(self):
        return self._closed

    def _finish(self):
        if self._finished:
            return Path(self.model.filename)
        if self._closed:
            raise RuntimeError("Cannot publish a closed acquisition")
        started = perf_counter()
        try:
            self.close_input()
            self.artifacts.collect(wait=True)
            write_root_attributes(self.model, {
                "artifact_status": json.dumps(self.artifacts.statuses, sort_keys=True),
                "artifact_failures": json.dumps(self.artifacts.failures, sort_keys=True),
                "artifact_timings": json.dumps(self.artifacts.timings, sort_keys=True),
            })
            if self.model._checkpoints:
                raise RuntimeError("Cannot finalize an active output transaction")
            # Final structural metadata may be written outside scientific stages;
            # this does not authorize any uncommitted scientific dataset.
            ProductCheckpoint(self.model).commit(tuple(
                path for path, record in self.model.records.items() if record.kind == "group"))
            local = self.output.finish()
            OutputPublisher(self.model.filename, self.model.storage.metrics,
                            cancelled=self._cancelled.is_set).publish(local)
            self._finished = True
            return Path(self.model.filename)
        finally:
            self.model.storage.metrics.add("output_finalization_seconds", perf_counter()-started)
            self.close()

    def close(self):
        self._cancelled.set()
        if not self._finished:
            self.output.request_abort()
        with self._lifecycle:
            if self._closed:
                return
            self._closed = True
            with ExitStack() as cleanup:
                cleanup.callback(self.model.storage.close)
                if self.owns_artifacts:
                    cleanup.callback(self.artifacts.__exit__, None, None, None)
                cleanup.callback(self.output.close)
                self.close_input()

    abort = close
