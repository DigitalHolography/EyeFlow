"""Isolated renderers, asynchronous publication, and rollback-safe artifact ownership."""

from __future__ import annotations

import multiprocessing
import os
import shutil
import traceback
from concurrent.futures import Future, ProcessPoolExecutor, ThreadPoolExecutor
from dataclasses import dataclass, replace
from pathlib import Path
from tempfile import TemporaryDirectory
from pipeline_engine.storage import temporary_directory
from threading import RLock
from uuid import uuid4

from input_output.output_manager import OutputManager
from input_output.artifact_ownership import ArtifactOwnership, current_artifact_transaction
from utils.logger import Logger


def _initialize_renderer() -> None:
    from runtime_limits import configure_numeric_threads

    configure_numeric_threads()
    import matplotlib

    matplotlib.use("Agg", force=True)


@dataclass
class ArtifactExportResult:
    paths: list[Path]
    messages: list[str]
    error: str | None = None


def _run_exporter(exporter, snapshot, manager, kwargs):
    messages = []
    with Logger.capture_messages(messages.append):
        try:
            result = exporter(snapshot, manager, **kwargs)
        except Exception:
            return ArtifactExportResult([], messages, traceback.format_exc())
    if isinstance(result, ArtifactExportResult):
        result.messages = messages + result.messages
        return result
    return ArtifactExportResult(list(result), messages)


class ArtifactWorkers:
    """Bounded rendering and publication pools, reused across a batch of inputs."""

    def __init__(self):
        self._renderer = None
        self._publisher = None
        self._admission = None
        self._lock = RLock()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        if self._admission is not None:
            self._admission.shutdown(wait=True)
        if self._renderer is not None:
            self._renderer.shutdown(wait=True, cancel_futures=True)
        if self._publisher is not None:
            self._publisher.shutdown(wait=True)

    def submit(self, job, complete):
        from runtime_limits import cap_parallel_jobs

        with self._lock:
            if self._renderer is None:
                self._renderer = ProcessPoolExecutor(
                    max_workers=cap_parallel_jobs(min(2, os.cpu_count() or 1)),
                    mp_context=multiprocessing.get_context("spawn"),
                    initializer=_initialize_renderer,
                )
                self._publisher = ThreadPoolExecutor(
                    max_workers=1, thread_name_prefix="EyeFlowArtifactPublisher",
                )
                self._admission = ThreadPoolExecutor(
                    max_workers=1, thread_name_prefix="EyeFlowExportAdmission",
                )
            from .resources import current_resources
            budget = current_resources()
            self._admission.submit(self._submit_admitted, job, complete, budget)

    def _submit_admitted(self, job, complete, budget):
        from time import monotonic
        memory = min(64 * 2**20, budget.memory_limit // 4) if budget else 0
        started = monotonic()
        if budget:
            while not budget.try_acquire(1, memory):
                with budget._condition:
                    budget._condition.wait(.1)
            if monotonic() - started > .1:
                Logger.log(f"[RESOURCES] Export {job.label} delayed {monotonic() - started:.2f}s.")
        try:
            rendered = self._renderer.submit(
                _run_exporter, job.exporter, job.snapshot, job.manager, job.kwargs,
            )
        except Exception as exc:
            if budget:
                budget.release(1, memory)
            rendered = Future()
            rendered.set_exception(exc)
            self._publisher.submit(complete, job, rendered)
            return

        def finished(future):
            # Publication only copies files; return compute slots before it.
            if budget:
                budget.release(1, memory)
            self._publisher.submit(complete, job, future)
        rendered.add_done_callback(finished)


@dataclass
class ArtifactJob:
    label: str
    temporary: TemporaryDirectory
    snapshot: Path
    manager: OutputManager
    destination: Path
    exporter: object
    kwargs: dict
    owner: object
    future: Future | None = None  # completes after publication and logging
    status: str = "prepared"


class BackgroundArtifacts:
    """Per-file artifact jobs over an optional shared batch worker pool.

    HDF5 handles stay on their owning thread. The publisher owns only staged
    files and in-memory job status. Every artifact has an explicit producing
    invocation, so rollback never touches another invocation's files.
    """

    def __init__(self, workers: ArtifactWorkers | None = None) -> None:
        self._workers = workers or ArtifactWorkers()
        self._owns_workers = workers is None
        self._jobs: list[ArtifactJob] = []
        self._lock = RLock()
        self.ownership = ArtifactOwnership()
        self._failures: dict[str, str] = {}

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.collect(wait=True)
        if self._owns_workers:
            self._workers.__exit__(*exc)
        for job in self._jobs:
            job.temporary.cleanup()

    @property
    def statuses(self) -> dict[str, str]:
        with self._lock:
            return {job.label: job.status for job in self._jobs}

    @property
    def failures(self) -> dict[str, str]:
        with self._lock:
            return dict(self._failures)

    @property
    def pending(self) -> bool:
        return any(status in {"prepared", "running"} for status in self.statuses.values())

    def prepare(self, label, manager, exporter, snapshotter, **kwargs) -> ArtifactJob:
        if label in self.statuses:
            raise ValueError(f"Duplicate background artifact job: {label}")
        temporary = temporary_directory(prefix="eyeflow-figures-")
        root = Path(temporary.name)
        staged = replace(manager, layout=replace(manager.layout, _root_dir=root / "render"))
        snapshot = root / "snapshot"
        try:
            snapshotter(snapshot)
        except BaseException:
            temporary.cleanup()
            raise
        job = ArtifactJob(
            label, temporary, snapshot, staged, manager.layout.ef_dir,
            exporter, kwargs,
            current_artifact_transaction() or self.ownership.transaction(manager.layout.ef_dir),
        )
        with self._lock:
            self._jobs.append(job)
        return job

    def cancel(self, job: ArtifactJob) -> None:
        with self._lock:
            if job.status != "prepared":
                raise RuntimeError("Only unsubmitted artifact jobs can be cancelled.")
            self._jobs.remove(job)
            job.status = "cancelled"
        job.temporary.cleanup()

    def start(self, job: ArtifactJob) -> None:
        job.future = Future()
        with self._lock:
            job.status = "running"
        Logger.log(f"[EXPORT] {job.manager.layout.stem} | {job.label}: generating artifacts in background.")
        try:
            self._workers.submit(job, self._complete)
        except Exception as exc:
            self._fail(job, exc)
            job.future.set_result(None)

    def collect(self, *, wait: bool = False) -> None:
        if wait:
            with self._lock:
                futures = [job.future for job in self._jobs if job.future is not None]
            for future in futures:
                future.result()

    def _complete(self, job, rendered):
        try:
            result = rendered.result()
            for message in result.messages:
                source = f"{job.manager.layout.stem} | {job.label}: "
                prefix, separator, body = message.partition("] ")
                if message.startswith("[") and separator:
                    Logger.log(f"{prefix}] {source}{body}")
                else:
                    Logger.log(f"{source}{message}")
            if result.error is not None:
                raise RuntimeError(result.error)
            self._publish(job, result.paths)
            with self._lock:
                job.status = "succeeded"
            Logger.log(
                f"[EXPORT COMPLETE] {job.manager.layout.stem} | {job.label}: "
                f"{len(result.paths)} artifact file(s) written."
            )
        except Exception as exc:
            self._fail(job, exc)
        finally:
            job.temporary.cleanup()
            job.future.set_result(None)

    def _fail(self, job: ArtifactJob, exc: Exception) -> None:
        message = str(exc)
        if exc.__cause__ is not None:
            message = "".join(traceback.format_exception(type(exc), exc, exc.__traceback__))
        with self._lock:
            job.status = "failed"
            self._failures[job.label] = message
        job.temporary.cleanup()
        Logger.log_error(f"[EXPORT FAILED] {job.manager.layout.stem} | {job.label}: {message}")

    def _publish(self, job, paths):
        staged_root = job.manager.layout.ef_dir.resolve()
        destination = job.destination.resolve()
        published = []
        try:
            for source in paths:
                source = Path(source).resolve()
                target = destination / source.relative_to(staged_root)
                if not target.resolve().is_relative_to(destination):
                    raise ValueError(f"Invalid artifact destination: {target}")
                target = target.resolve()
                # Keep the temporary basename short: repeating the destination
                # name can exceed Windows MAX_PATH even when the target fits.
                temporary = target.with_name(f".ef-{uuid4().hex[:16]}.partial")
                if target.exists():
                    raise FileExistsError(f"Artifact already exists: {target}")
                job.owner.reserve(target, writing=True)
                job.owner.reserve(temporary, writing=True)
                target.parent.mkdir(parents=True, exist_ok=True)
                try:
                    shutil.copyfile(source, temporary)
                    os.replace(temporary, target)
                    published.append(target)
                finally:
                    temporary.unlink(missing_ok=True)
        except BaseException:
            for path in published:
                path.unlink(missing_ok=True)
            raise
