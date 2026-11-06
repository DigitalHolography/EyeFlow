"""Isolated renderers, asynchronous publication, and rollback-safe artifact ownership."""

from __future__ import annotations

import multiprocessing
import os
import shutil
import traceback
from concurrent.futures import Future, ProcessPoolExecutor, ThreadPoolExecutor
from dataclasses import dataclass, field, replace
from pathlib import Path
from tempfile import TemporaryDirectory
from scratch_location import temporary_directory
from threading import RLock
from uuid import uuid4
from time import perf_counter

from input_output.output_manager import OutputManager
from input_output.artifact_ownership import ArtifactOwnership, current_artifact_transaction
from utils.logger import Logger


def _initialize_renderer() -> None:
    from runtime_limits import configure_native_threads, set_thread_allowance

    configure_native_threads(1)
    set_thread_allowance(1)
    import matplotlib

    matplotlib.use("Agg", force=True)


@dataclass
class ArtifactExportResult:
    paths: list[Path]
    messages: list[str]
    error: str | None = None
    rendering_seconds: float = 0.0
    started_at: float = 0.0


def _run_exporter(label, exporter, snapshot, manager, kwargs):
    started = perf_counter()
    messages = []
    with Logger.capture_messages(messages.append):
        Logger.log_perf(f"START artifact={label} stage=render at={started:.6f}")
        try:
            result = exporter(snapshot, manager, **kwargs)
        except Exception:
            duration = perf_counter() - started
            Logger.log_perf(
                f"END artifact={label} stage=render status=failed "
                f"at={perf_counter():.6f} wall={duration:.3f}s"
            )
            return ArtifactExportResult([], messages, traceback.format_exc(), duration, started)
        duration = perf_counter() - started
        Logger.log_perf(
            f"END artifact={label} stage=render status=succeeded "
            f"at={perf_counter():.6f} wall={duration:.3f}s"
        )
    if isinstance(result, ArtifactExportResult):
        result.messages = messages + result.messages
        result.rendering_seconds = duration
        result.started_at = started
        return result
    return ArtifactExportResult(list(result), messages, rendering_seconds=duration, started_at=started)


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
            from pipeline_engine.resources import current_resources
            budget = current_resources()
            self._admission.submit(self._submit_admitted, job, complete, budget)

    def _submit_admitted(self, job, complete, budget):
        admitted = False
        admission_started = perf_counter()
        try:
            resource_wait_started = perf_counter()
            request = budget.default_request(threads=1) if budget else None
            if budget:
                budget.acquire(request.threads, request.memory)
                admitted = True
            renderer_submitted_at = perf_counter()
            job.timings.update(
                admission_queue_seconds=admission_started - job.timings["submitted_at"],
                resource_admission_wait_seconds=renderer_submitted_at - resource_wait_started,
                admission_wait_seconds=renderer_submitted_at - job.timings["submitted_at"],
                renderer_submitted_at=renderer_submitted_at,
            )
            if job.timings["resource_admission_wait_seconds"] > .1:
                Logger.log_perf(
                    f"artifact={job.label} stage=resource_admission "
                    f"wait={job.timings['resource_admission_wait_seconds']:.3f}s"
                )
            rendered = self._renderer.submit(
                _run_exporter, job.label, job.exporter, job.snapshot, job.manager, job.kwargs,
            )
        except Exception as exc:
            job.timings["admission_wait_seconds"] = perf_counter() - job.timings["submitted_at"]
            job.timings["admission_queue_seconds"] = max(
                0.0, admission_started - job.timings["submitted_at"]
            )
            if admitted:
                budget.release(request.threads, request.memory)
            rendered = Future()
            rendered.set_exception(exc)
            self._schedule_completion(job, complete, rendered)
            return

        def finished(future):
            # Publication only copies files; return compute slots before it.
            try:
                if admitted:
                    budget.release(request.threads, request.memory)
            finally:
                self._schedule_completion(job, complete, future)
        rendered.add_done_callback(finished)

    def _schedule_completion(self, job, complete, rendered):
        def finish():
            try:
                complete(job, rendered)
            except BaseException as exc:
                if job.future is not None and not job.future.done():
                    job.future.set_exception(exc)

        try:
            self._publisher.submit(finish)
        except Exception:
            # Pool shutdown/submission failure must still settle the job.
            finish()


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
    timings: dict = field(default_factory=dict)


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
    def timings(self):
        with self._lock:
            return {job.label: dict(job.timings) for job in self._jobs}

    @property
    def pending(self) -> bool:
        return any(status in {"prepared", "running"} for status in self.statuses.values())

    def prepare(self, label, manager, exporter, snapshotter, **kwargs) -> ArtifactJob:
        started = perf_counter()
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
        job.timings.update(prepared_at=started, preparation_seconds=perf_counter() - started)
        with self._lock:
            self._jobs.append(job)
        Logger.log_perf(
            f"START artifact={label} stage=snapshot_preparation at={started:.6f}"
        )
        Logger.log_perf(
            f"END artifact={label} stage=snapshot_preparation "
            f"at={perf_counter():.6f} wall={job.timings['preparation_seconds']:.3f}s"
        )
        return job

    def cancel(self, job: ArtifactJob) -> None:
        with self._lock:
            if job.status != "prepared":
                raise RuntimeError("Only unsubmitted artifact jobs can be cancelled.")
            self._jobs.remove(job)
            job.status = "cancelled"
        job.temporary.cleanup()

    def start(self, job: ArtifactJob) -> None:
        job.timings["submitted_at"] = perf_counter()
        job.future = Future()
        with self._lock:
            job.status = "running"
        Logger.log(f"[EXPORT] {job.manager.layout.stem} | {job.label}: generating artifacts in background.")
        Logger.log_perf(
            f"START artifact={job.label} stage=queue_wait "
            f"at={job.timings['submitted_at']:.6f}"
        )
        try:
            self._workers.submit(job, self._complete)
        except Exception as exc:
            try:
                self._fail(job, exc)
            finally:
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
            job.timings.update(rendering_seconds=result.rendering_seconds,
                               render_started_at=result.started_at,
                               queue_seconds=max(
                                   0.0, result.started_at - job.timings["submitted_at"]
                               ),
                               renderer_queue_seconds=max(
                                   0.0, result.started_at - job.timings["renderer_submitted_at"]
                               ))
            for message in result.messages:
                source = f"{job.manager.layout.stem} | {job.label}: "
                prefix, separator, body = message.partition("] ")
                if message.startswith("[") and separator:
                    Logger.log(f"{prefix}] {source}{body}")
                else:
                    Logger.log(f"{source}{message}")
            if result.error is not None:
                raise RuntimeError(result.error)
            started = perf_counter()
            job.timings["publication_started_at"] = started
            Logger.log_perf(f"START artifact={job.label} stage=publication at={started:.6f}")
            self._publish(job, result.paths)
            job.timings["publication_seconds"] = perf_counter() - started
            Logger.log_perf(
                f"END artifact={job.label} stage=publication at={perf_counter():.6f} "
                f"wall={job.timings['publication_seconds']:.3f}s "
                f"file_copy={job.timings.get('copy_seconds', 0.0):.3f}s "
                f"atomic_replace={job.timings.get('atomic_replace_seconds', 0.0):.3f}s"
            )
            with self._lock:
                job.status = "succeeded"
            Logger.log(
                f"[EXPORT COMPLETE] {job.manager.layout.stem} | {job.label}: "
                f"{len(result.paths)} artifact file(s) written."
            )
        except Exception as exc:
            self._fail(job, exc)
        finally:
            job.timings["finished_at"] = perf_counter()
            Logger.log_perf(
                f"END artifact={job.label} status={job.status} "
                f"at={job.timings['finished_at']:.6f} "
                f"preparation={job.timings.get('preparation_seconds', 0.0):.3f}s "
                f"admission_queue={job.timings.get('admission_queue_seconds', 0.0):.3f}s "
                f"resource_wait={job.timings.get('resource_admission_wait_seconds', 0.0):.3f}s "
                f"renderer_queue={job.timings.get('renderer_queue_seconds', 0.0):.3f}s "
                f"render={job.timings.get('rendering_seconds', 0.0):.3f}s "
                f"publication={job.timings.get('publication_seconds', 0.0):.3f}s "
                f"file_copy={job.timings.get('copy_seconds', 0.0):.3f}s"
            )
            try:
                job.temporary.cleanup()
            finally:
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
                    copy_started = perf_counter()
                    shutil.copyfile(source, temporary)
                    job.timings["copy_seconds"] = (
                        job.timings.get("copy_seconds", 0.0) + perf_counter() - copy_started
                    )
                    replace_started = perf_counter()
                    os.replace(temporary, target)
                    job.timings["atomic_replace_seconds"] = (
                        job.timings.get("atomic_replace_seconds", 0.0)
                        + perf_counter() - replace_started
                    )
                    published.append(target)
                finally:
                    temporary.unlink(missing_ok=True)
        except BaseException:
            for path in published:
                path.unlink(missing_ok=True)
            raise
