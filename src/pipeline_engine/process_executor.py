"""Own spawned invocations and accept their results on the coordinator thread."""

import json
import shutil
from dataclasses import dataclass, fields, is_dataclass, replace
from pathlib import Path
from queue import Empty, Queue
from threading import Event, Thread
from time import monotonic
from uuid import uuid4

from utils.logger import Logger

from .protocol import InvocationSpec
from .transfer import dump_result, load_result
from .worker import run_worker


def published_paths(value, source_root, destination_root):
    """Translate artifact references; numerical scratch remains local until finish."""
    if isinstance(value, Path):
        path = value.resolve()
        return (
            destination_root / path.relative_to(source_root)
            if path.is_relative_to(source_root) else value
        )
    if isinstance(value, dict):
        return {
            key: published_paths(item, source_root, destination_root)
            for key, item in value.items()
        }
    if isinstance(value, list):
        return [published_paths(item, source_root, destination_root) for item in value]
    if isinstance(value, tuple) and type(value) is tuple:
        return tuple(published_paths(item, source_root, destination_root) for item in value)
    if is_dataclass(value) and not isinstance(value, type):
        return replace(value, **{
            field.name: published_paths(getattr(value, field.name), source_root, destination_root)
            for field in fields(value) if field.init
        })
    return value


def _result_attributes(model, root, name, stages, error):
    """Combine worker stage metadata without replacing other invocations' status."""
    stage_keys = {"stage_status", "stage_failures", "stage_skips", "stage_products"}
    attributes = {key: value for key, value in root.attrs.items() if key not in stage_keys}
    for key in stage_keys:
        incoming = json.loads(root.attrs.get(key, "{}"))
        if key == "stage_status":
            if stages:
                incoming.setdefault(name, {}).update(stages)
            if error and name in incoming:
                incoming[name] = {
                    stage: {"running": "failed", "pending": "skipped"}.get(status, status)
                    for stage, status in incoming[name].items()
                }
        if incoming:
            existing = json.loads(model.attrs.get(key, "{}"))
            existing.update(incoming)
            attributes[key] = json.dumps(existing, sort_keys=True)
    return attributes


def accept_result(result, *, model, state, destination, name, stages, error=None):
    """Accept records, bindings and files together, rolling back on any failure."""
    source_root, destination = result.artifact_root.resolve(), destination.resolve()
    records = {path: record for path, record in result.records.items() if path}
    model.validate_merge(records)
    attributes = _result_attributes(model, result.records[""], name, stages, error)
    changes = published_paths(result.state, source_root, destination)
    artifacts = []
    targets = set()
    for path in result.artifacts:
        source = Path(path).resolve()
        if not source.is_relative_to(source_root) or not source.is_file():
            raise ValueError(f"Invalid worker artifact: {source}")
        target = (destination / source.relative_to(source_root)).resolve()
        if not target.is_relative_to(destination):
            raise ValueError(f"Invalid artifact destination: {target}")
        if target.exists() or target in targets:
            raise FileExistsError(f"Pipeline artifact already exists: {target}")
        targets.add(target)
        artifacts.append((source, target))

    published = []
    previous_state = dict(state)
    try:
        with model.checkpoint():
            model.merge(records)
            for key, value in attributes.items():
                model.attrs[key] = value
            for source, target in artifacts:
                target.parent.mkdir(parents=True, exist_ok=True)
                incoming = target.with_name(f".ef-{uuid4().hex[:16]}{target.suffix}")
                try:
                    shutil.copyfile(source, incoming)
                    # Unlike replace(), rename refuses to overwrite a concurrently
                    # published artifact on Windows. Destinations were preflighted.
                    if target.exists():
                        raise FileExistsError(f"Pipeline artifact already exists: {target}")
                    incoming.rename(target)
                    published.append(target)
                finally:
                    incoming.unlink(missing_ok=True)
            state.update(changes)
    except BaseException:
        state.clear()
        state.update(previous_state)
        for path in published:
            path.unlink(missing_ok=True)
        raise
    for message in result.messages:
        Logger.log(message.replace(str(source_root), str(destination)))


@dataclass
class ActiveInvocation:
    process: object
    threads: int
    memory: int
    memory_account: object
    stages: dict
    last_result: str | None = None
    finished: bool = False
    exited_at: float | None = None


@dataclass(frozen=True)
class InvocationOutcome:
    name: str
    error: str | None = None


class ProcessExecutor:
    """A per-input process owner; scheduling decisions belong to the caller."""

    def __init__(self, *, context, input_service, inputs, model, manager, budget, state, settings):
        self.context, self.input_service, self.inputs = context, input_service, inputs
        self.model, self.manager, self.budget = model, manager, budget
        self.state, self.settings = state, settings
        self.events = context.Queue(maxsize=32)
        self._incoming = Queue()
        self._stop = Event()
        self._acks = {}
        self.active = {}
        self._receiver = Thread(target=self._receive, name="EyeFlowResults", daemon=True)
        self._receiver.start()

    @property
    def running(self):
        return any(not job.finished for job in self.active.values())

    @property
    def coordinator_ready(self):
        """Keep collecting worker completions instead of blocking on their slots."""
        if not self.running:
            return True
        request = self.budget.default_request(threads=1)
        if not self.budget.try_acquire(request.threads, request.memory):
            return False
        self.budget.release(request.threads, request.memory)
        return True

    def _receive(self):
        # Drain and acknowledge while coordinator-owned scientific code runs.
        while not self._stop.is_set():
            try:
                event = self.events.get(timeout=0.1)
            except Empty:
                continue
            if event.kind == "LogMessage":
                try:
                    Logger.log(event.message or "")
                except Exception:
                    # Presentation callbacks must not strand the worker protocol.
                    pass
            else:
                self._incoming.put(event)
            self._acks[event.invocation].put(True)

    def start(self, pipeline):
        request = self.budget.default_request(threads=min(4, max(1, self.budget.cpu_limit // 2)))
        if not self.budget.try_acquire(request.threads, request.memory):
            return False
        try:
            self.active[pipeline.name] = self._launch(pipeline, request.threads, request.memory)
        except BaseException:
            self.budget.release(request.threads, request.memory)
            raise
        Logger.log(f"[START] {pipeline.name} in process {self.active[pipeline.name].process.pid}, "
                   f"{request.threads} compute slots.")
        return True

    def _launch(self, pipeline, threads, memory):
        directory = self.model.storage.directory / f"worker-{uuid4().hex}"
        directory.mkdir()
        instance = pipeline.instantiate()
        entry = getattr(instance, "func", type(instance))
        initial = dump_result(self.state, self.model.storage)
        acknowledgements = self.context.Queue(maxsize=1)
        self._acks[pipeline.name] = acknowledgements
        account = self.budget.memory_account()
        spec = InvocationSpec(
            name=pipeline.name, module=entry.__module__, directory=directory,
            output_manager=self.manager,
            input_metadata={slot: self.input_service.metadata.get(slot) for slot in ("hd", "dv")},
            input_client=self.input_service.client(), settings=self.settings,
            initial_result=initial, threads=threads, memory_account=account,
            acknowledgements=acknowledgements, native_parallelism=pipeline.native_parallelism,
        )
        process = self.context.Process(
            target=run_worker, args=(spec, self.events), name=f"EyeFlow-{pipeline.name}",
        )
        try:
            process.start()
        except BaseException:
            process.close()
            account.reclaim()
            raise
        return ActiveInvocation(
            process=process, threads=threads, memory=memory, memory_account=account, stages={},
        )

    def _finish(self, name, error):
        active = self.active[name]
        if active.finished:
            return None
        active.finished = True
        try:
            if active.last_result is None:
                error = error or "Worker exited without results."
            else:
                result = load_result(active.last_result, self.model.storage, self.inputs)
                accept_result(
                    result, model=self.model, state=self.state,
                    destination=self.manager.layout.ef_dir, name=name,
                    stages=active.stages, error=error,
                )
        except Exception as exc:
            error = str(exc)
        finally:
            self.budget.release(active.threads, active.memory)
        return InvocationOutcome(name, error)

    def collect(self, *, wait=False):
        outcomes = []
        try:
            event = self._incoming.get(timeout=0.1 if wait else 0.001)
        except Empty:
            event = None
        while event is not None:
            active = self.active[event.invocation]
            status = {"StageStarted": "running", "StageFinished": "succeeded",
                      "StageFailed": "failed", "StageSkipped": "skipped"}.get(event.kind)
            if status:
                active.stages[event.stage] = status
            if event.result_path:
                active.last_result = event.result_path
            if event.kind == "InvocationFinished":
                outcome = self._finish(event.invocation, event.error)
                if outcome:
                    outcomes.append(outcome)
            elif event.kind in {"StageFailed", "StageSkipped"}:
                Logger.log_warning(f"{event.invocation}/{event.stage}: {event.error}")
            try:
                event = self._incoming.get_nowait()
            except Empty:
                event = None
        for name, active in self.active.items():
            if active.process.is_alive():
                continue
            active.process.join(timeout=0)
            active.memory_account.reclaim()
            if not active.finished:
                # The receiver can be between queue receipt and local delivery.
                active.exited_at = active.exited_at or monotonic()
                if monotonic() - active.exited_at > 0.2:
                    error = f"Worker exited with code {active.process.exitcode}."
                    outcomes.append(self._finish(name, error))
        return outcomes

    def close(self):
        for active in self.active.values():
            if active.process.is_alive() and not active.finished:
                active.process.terminate()
            active.process.join(timeout=5)
            if active.process.is_alive():
                active.process.kill()
                active.process.join()
            if not active.finished:
                self.budget.release(active.threads, active.memory)
            active.process.close()
            active.memory_account.reclaim()
        self._stop.set()
        self._receiver.join()
        for acknowledgements in self._acks.values():
            acknowledgements.close()
            acknowledgements.join_thread()
        self.events.close()
        self.events.join_thread()
