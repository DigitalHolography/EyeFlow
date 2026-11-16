"""One owner for acquisition exports, final metadata and atomic publication."""

from contextlib import ExitStack
import json
from pathlib import Path
from threading import Event, RLock
from time import perf_counter

from input_output.writers.h5 import write_root_attributes
from input_output.writers.h5_publication import OutputPublisher
from pipeline_engine.execution.background_artifacts import BackgroundArtifacts
from pipeline_engine.io.manager import IOManager
from utils.logger import Logger


class AcquisitionCompletion:
    """Retain accepted backing until exports and publication settle.

    Construction transfers ownership of the I/O lanes and per-acquisition jobs.
    Batch renderer pools remain borrowed through BackgroundArtifacts. Standalone
    publication explicitly borrows the caller's storage instead of releasing it.
    A failed constructor leaves all caller-owned storage untouched.
    """

    def __init__(self, model, *, artifacts=None, owns_storage=True,
                 budget=None, context=None, output_policy=None):
        self.model = model
        self.artifacts = artifacts if artifacts is not None else BackgroundArtifacts()
        self.owns_storage = owns_storage
        self._lifecycle, self._cancelled = RLock(), Event()
        self._closed = self._finished = False
        self.io = IOManager(
            model, budget=budget, context=context, output_policy=output_policy)

    @property
    def closed(self):
        return self._closed

    def collect_artifacts(self, *, wait=False):
        """Reconcile artifact diagnostics through a single metadata authority."""
        started = perf_counter()
        if wait:
            Logger.log_perf(f"START pipeline=run variant=shared stage=background_artifact_join at={started:.6f}")
        try:
            self.artifacts.collect(wait=wait)
            if wait or self.artifacts.statuses:
                write_root_attributes(self.model, {
                    "artifact_status": json.dumps(self.artifacts.statuses, sort_keys=True),
                    "artifact_failures": json.dumps(self.artifacts.failures, sort_keys=True),
                    "artifact_timings": json.dumps(self.artifacts.timings, sort_keys=True),
                })
        finally:
            if wait:
                Logger.log_perf(f"END pipeline=run variant=shared stage=background_artifact_join "
                                f"at={perf_counter():.6f} wall={perf_counter()-started:.3f}s")

    def finish(self):
        with self._lifecycle:
            if self._finished:
                return Path(self.model.filename)
            if self._closed:
                raise RuntimeError("Cannot publish a closed acquisition")
            started = perf_counter()
            try:
                self.io.close_input()
                self.collect_artifacts(wait=True)
                snapshot = self.model.finalization_snapshot()
                completed = self.io.output.finish(snapshot)
                OutputPublisher(self.model.filename, self.model.storage.metrics,
                                cancelled=self._cancelled.is_set).publish(completed)
                self._finished = True
                return Path(self.model.filename)
            finally:
                self.model.storage.metrics.add("output_finalization_seconds", perf_counter()-started)
                self.close()

    def close(self):
        # Wake a finalizer waiting on the writer before taking its lifecycle lock.
        self._cancelled.set()
        if not self._finished:
            self.io.output.request_abort()
        with self._lifecycle:
            if self._closed:
                return
            self._closed = True
            with ExitStack() as cleanup:
                if self.owns_storage:
                    cleanup.callback(self.model.storage.close)
                cleanup.callback(self.artifacts.__exit__, None, None, None)
                cleanup.callback(self.io.close)

    abort = close


def publish_output_model(model):
    """Standalone facade using the same completion path and borrowed storage."""
    model.commit_pending()
    completion = AcquisitionCompletion(model, owns_storage=False)
    return completion.finish()
