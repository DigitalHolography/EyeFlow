"""Dependency-ready scheduling of opt-in processes and sequential small jobs."""

import shutil
from dataclasses import dataclass, fields, is_dataclass, replace
from pathlib import Path
from queue import Empty, Queue
from threading import Event, Thread
from time import monotonic
from uuid import uuid4

from utils.logger import Logger

from .transfer import dump_result, load_result
from .worker import InvocationSpec, run_worker


def published_paths(value, source_root, destination_root):
    """Rewrite artifact catalog references, keeping local numerical scratch paths."""
    if isinstance(value, Path):
        path = value.resolve()
        return (
            destination_root / path.relative_to(source_root)
            if path.is_relative_to(source_root)
            else value
        )
    if isinstance(value, dict):
        return {
            key: published_paths(item, source_root, destination_root) for key, item in value.items()
        }
    if isinstance(value, list):
        return [published_paths(item, source_root, destination_root) for item in value]
    if isinstance(value, tuple) and type(value) is tuple:
        return tuple(published_paths(item, source_root, destination_root) for item in value)
    if is_dataclass(value) and not isinstance(value, type):
        return replace(
            value,
            **{
                field.name: published_paths(
                    getattr(value, field.name), source_root, destination_root
                )
                for field in fields(value)
                if field.init
            },
        )
    return value


@dataclass
class ActiveInvocation:
    process: object
    directory: Path
    threads: int
    memory: int
    memory_account: object
    acknowledgements: object
    stages: dict
    last_result: str | None = None
    finished: bool = False
    exited_at: float | None = None


class CoordinatorScheduler:
    def __init__(
        self,
        *,
        pipelines,
        dag,
        options,
        context,
        input_service,
        inputs,
        model,
        manager,
        budget,
        state,
        statuses,
        failures,
        skips,
        settings,
        on_start=None,
        on_success=None,
        on_progress=None,
    ):
        self.pending = list(pipelines)
        self.dag, self.options, self.context = dag, options, context
        self.input_service, self.inputs = input_service, inputs
        self.model, self.manager, self.budget = model, manager, budget
        self.state, self.statuses, self.failures = state, statuses, failures
        self.skips = skips
        self.settings = settings
        self.on_start, self.on_success, self.on_progress = on_start, on_success, on_progress
        self.events = context.Queue(maxsize=32)
        self._incoming = Queue()
        self._stop = Event()
        self._acks = {}
        self._receiver = Thread(target=self._receive, name="EyeFlowResults", daemon=True)
        self._receiver.start()
        self.active = {}
        self.completed = set()
        self.total = len(pipelines)
        self.started = 0

    def _receive(self):
        """Keep IPC draining while coordinator-owned calculations are running."""
        while not self._stop.is_set():
            try:
                event = self.events.get(timeout=0.1)
            except Empty:
                continue
            self._incoming.put(event)
            self._acks[event.invocation].put(True)

    def _dependencies(self, pipeline):
        return self.dag.dependencies_of(pipeline.name, pipeline_options=self.options)

    def _start(self, pipeline):
        # The first opt-in production worker is independent displacement. Other
        # process pipelines may consume completed shared producers, but variant
        # fan-out stays in the established coordinator path for this rollout.
        self.budget.refresh()
        threads = max(1, min(4, getattr(self.budget, "_dynamic_cpu", self.budget.cpu_limit) // 2))
        memory = min(64 * 2**20, self.budget.memory_limit // 4)
        if not self.budget.try_acquire(threads, memory):
            return False
        try:
            active = self._launch(pipeline, threads, memory)
        except Exception as exc:  # noqa: BLE001 - isolate spawning/serialization failures
            self.budget.release(threads, memory)
            self.statuses[pipeline.name] = "failed"
            self.failures[pipeline.name] = str(exc)
            self.pending.remove(pipeline)
            self.completed.add(pipeline.name)
            self.started += 1
            if self.on_start:
                self.on_start(pipeline.name, self.started, self.total)
            if self.on_progress:
                self.on_progress()
            return True
        except BaseException:
            self.budget.release(threads, memory)
            raise
        self.active[pipeline.name] = active
        self.statuses[pipeline.name] = "running"
        self.pending.remove(pipeline)
        self.started += 1
        if self.on_start:
            self.on_start(pipeline.name, self.started, self.total)
        Logger.log(
            f"[START] {pipeline.name} in process {active.process.pid}, {threads} compute slots."
        )
        return True

    def _launch(self, pipeline, threads, memory):
        directory = self.model.storage.directory / f"worker-{uuid4().hex}"
        directory.mkdir()
        instance = pipeline.instantiate()
        entry = getattr(instance, "func", type(instance))
        module = entry.__module__
        initial = dump_result(self.state, self.model.storage)
        client = self.input_service.client()
        metadata = {slot: self.input_service.metadata.get(slot) for slot in ("hd", "dv")}
        account = self.budget.memory_account()
        acknowledgements = self.context.Queue(maxsize=1)
        self._acks[pipeline.name] = acknowledgements
        spec = InvocationSpec(
            pipeline.name,
            module,
            str(directory),
            self.manager,
            metadata,
            client,
            self.settings,
            initial,
            threads,
            account,
            acknowledgements,
            native_parallelism=pipeline.native_parallelism,
        )
        process = self.context.Process(
            target=run_worker, args=(spec, self.events), name=f"EyeFlow-{pipeline.name}"
        )
        try:
            process.start()
        except BaseException:
            if not process.is_alive():
                process.close()
            raise
        return ActiveInvocation(process, directory, threads, memory, account, acknowledgements, {})

    def _merge(self, name, result_path, error):
        result = load_result(result_path, self.model.storage, self.inputs)
        records = result["records"]
        # Root statuses are aggregated by the coordinator, not overwritten by
        # whichever process finishes last.
        root = records.pop("")
        import json

        stage_attributes = {"stage_status", "stage_failures", "stage_skips", "stage_products"}
        for key, value in root.attrs.items():
            if key not in stage_attributes:
                self.model.attrs[key] = value
        for attribute in ("stage_status", "stage_failures", "stage_skips", "stage_products"):
            incoming = json.loads(root.attrs.get(attribute, "{}"))
            if attribute == "stage_status":
                stages = dict(self.active[name].stages)
                if error:
                    stages = {
                        stage: "failed" if status == "running" else status
                        for stage, status in stages.items()
                    }
                # Legacy invocation wrapper events use the pipeline name;
                # only actual graph stages belong in stage_status.
                stages.pop(name, None)
                if stages:
                    incoming.setdefault(name, {}).update(stages)
                if error and name in incoming:
                    incoming[name] = {
                        stage: "failed"
                        if status == "running"
                        else "skipped"
                        if status == "pending"
                        else status
                        for stage, status in incoming[name].items()
                    }
            if incoming:
                existing = json.loads(self.model.attrs.get(attribute, "{}"))
                existing.update(incoming)
                self.model.attrs[attribute] = json.dumps(existing, sort_keys=True)
        self.model.merge(records)
        staged_root = result["manager"].layout.ef_dir.resolve()
        self.state.update(
            published_paths(result["state"], staged_root, self.manager.layout.ef_dir.resolve())
        )
        for message in result.get("messages", ()):
            Logger.log(message.replace(str(staged_root), str(self.manager.layout.ef_dir.resolve())))
        for path in result["artifacts"]:
            source = Path(path).resolve()
            if not source.is_relative_to(staged_root):
                raise ValueError(f"Worker artifact outside its staging directory: {source}")
            relative = source.relative_to(staged_root)
            target = self.manager.layout.ef_dir / relative
            if target.exists():
                raise FileExistsError(f"Pipeline artifact already exists: {target}")
            target.parent.mkdir(parents=True, exist_ok=True)
            incoming = target.with_name(f".ef-{uuid4().hex[:16]}{target.suffix}")
            try:
                shutil.copyfile(source, incoming)
                incoming.replace(target)
            finally:
                incoming.unlink(missing_ok=True)
        self.statuses[name] = "failed" if error else "succeeded"
        if error:
            self.failures[name] = error
            Logger.log_warning(f"Failed {name}; keeping completed products: {error}")
        elif self.on_success:
            self.on_success(name)

    def _finish(self, name, error=None):
        active = self.active[name]
        if active.finished:
            return
        active.finished = True
        try:
            if active.last_result:
                self._merge(name, active.last_result, error)
            else:
                self.statuses[name] = "failed"
                self.failures[name] = error or "Worker exited without results."
        except Exception as exc:  # noqa: BLE001 - a bad result must not cancel independent pipelines
            self.statuses[name] = "failed"
            self.failures[name] = str(exc)
        finally:
            self.budget.release(active.threads, active.memory)
            self.completed.add(name)
            if self.on_progress:
                self.on_progress()

    def collect(self, *, wait=False):
        try:
            event = self._incoming.get(timeout=0.1 if wait else 0.001)
        except Empty:
            event = None
        while event is not None:
            active = self.active[event.invocation]
            status = {
                "StageStarted": "running",
                "StageFinished": "succeeded",
                "StageFailed": "failed",
                "StageSkipped": "skipped",
            }.get(event.kind)
            if status:
                active.stages[event.stage] = status
            if event.kind == "ProductCommit":
                active.last_result = event.detail["result"]
            elif event.kind == "InvocationFinished":
                active.last_result = event.detail["result"]
                self._finish(event.invocation, event.detail.get("error"))
            elif event.kind in {"StageFailed", "StageSkipped"}:
                Logger.log_warning(f"{event.invocation}/{event.stage}: {event.detail.get('error')}")
            try:
                event = self._incoming.get_nowait()
            except Empty:
                event = None
        for name, active in self.active.items():
            if not active.process.is_alive():
                active.process.join(timeout=0)
                active.memory_account.reclaim()
            if not active.process.is_alive() and not active.finished:
                # Queue feeder messages can arrive just after process exit.
                active.exited_at = active.exited_at or monotonic()
                if monotonic() - active.exited_at > 0.2:
                    self._finish(name, f"Worker exited with code {active.process.exitcode}.")

    def coordinator_jobs(self):
        while self.pending or any(not job.finished for job in self.active.values()):
            self.collect()
            for pipeline in tuple(self.pending):
                if pipeline.execution_mode != "process":
                    continue
                dependencies = self._dependencies(pipeline)
                if not set(dependencies) <= self.completed:
                    continue
                if any(self.statuses.get(name) in {"failed", "skipped"} for name in dependencies):
                    self.pending.remove(pipeline)
                    self.statuses[pipeline.name] = "skipped"
                    self.skips[pipeline.name] = "Unavailable dependencies: " + ", ".join(
                        dependencies
                    )
                    self.completed.add(pipeline.name)
                    self.started += 1
                    if self.on_start:
                        self.on_start(pipeline.name, self.started, self.total)
                    if self.on_progress:
                        self.on_progress()
                    continue
                self._start(pipeline)
            selected = next(
                (
                    pipeline
                    for pipeline in self.pending
                    if pipeline.execution_mode == "coordinator"
                    and set(self._dependencies(pipeline)) <= self.completed
                    and (
                        not pipeline.wait_for_artifacts
                        or (
                            not any(not p.wait_for_artifacts for p in self.pending)
                            and not any(not job.finished for job in self.active.values())
                        )
                    )
                ),
                None,
            )
            if selected is not None:
                self.pending.remove(selected)
                self.started += 1
                yield self.started, selected
                self.completed.add(selected.name)
            elif any(not job.finished for job in self.active.values()):
                self.collect(wait=True)
            elif self.pending:
                # Resource pressure is transient; allow monitoring to refresh.
                if any(
                    p.execution_mode == "process" and set(self._dependencies(p)) <= self.completed
                    for p in self.pending
                ):
                    self.collect(wait=True)
                else:
                    raise RuntimeError("No schedulable pipeline; inspect active dependencies.")

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
