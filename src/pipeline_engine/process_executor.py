"""Spawn, transport receipt, worker supervision and invocation outcomes."""

from dataclasses import dataclass
from queue import Empty, Queue
from threading import Event, Thread
from time import monotonic
from uuid import uuid4
import pickle

from utils.logger import Logger

from .protocol import InvocationSpec, TransportReceipt
from .transfer import dump_result
from .worker import run_worker
from .product_exchange import PipelineInvocationId


@dataclass
class ActiveInvocation:
    process: object
    threads: int
    memory: int
    memory_account: object
    last_result: str | None = None
    acceptance_error: str | None = None
    finished: bool = False
    exited_at: float | None = None


@dataclass(frozen=True)
class InvocationOutcome:
    name: str
    error: str | None = None


class ProcessExecutor:
    """A per-input process owner; scheduling decisions belong to the caller."""

    def __init__(self, *, context, input_service, inputs, model, manager, budget, state, settings,
                 acceptance, registry, run_id):
        self.context, self.input_service, self.inputs = context, input_service, inputs
        self.model, self.manager, self.budget = model, manager, budget
        self.state, self.settings = state, settings
        self.acceptance, self.registry, self.run_id = acceptance, registry, run_id
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
            self.model.storage.metrics.add("ipc_events", 1)
            self.model.storage.metrics.add("ipc_bytes", len(pickle.dumps(event, protocol=5)))
            self._acks[event.invocation].put(TransportReceipt(event.commit_id))

    def invocation_id(self, pipeline, variant=None):
        return PipelineInvocationId(self.run_id, pipeline.name, variant)

    def product_inputs(self, pipeline, variant=None):
        identity = self.invocation_id(pipeline, variant)
        return {item.alias: self.registry.get(item.resolve(identity)) for item in pipeline.product_requires}

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
            invocation_id=self.invocation_id(pipeline), product_inputs=self.product_inputs(pipeline),
            owner_id=str(self.model.storage.directory.resolve()),
        )
        process = self.context.Process(
            target=run_worker, args=(spec, self.events), name=f"EyeFlow-{pipeline.name}",
        )
        try:
            process.start()
            self.model.storage.metrics.add("processes_started", 1)
        except BaseException:
            process.close()
            account.reclaim()
            raise
        return ActiveInvocation(
            process=process, threads=threads, memory=memory, memory_account=account,
        )

    def _finish(self, name, error):
        active = self.active[name]
        if active.finished:
            return None
        active.finished = True
        try:
            error = active.acceptance_error or error
            if active.last_result is None:
                error = error or "Worker exited without results."
            self.acceptance.finalize(PipelineInvocationId(self.run_id, name), error)
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
            if active.finished:
                event = self._incoming.get_nowait() if not self._incoming.empty() else None
                continue
            identity = PipelineInvocationId(self.run_id, event.invocation)
            try:
                if event.invocation_id != identity:
                    raise ValueError("Worker event invocation identity mismatch.")
                if event.result_path:
                    self.acceptance.accept_event(event, identity)
                    active.last_result = event.result_path
                self.acceptance.observe(identity, event)
            except Exception as exc:
                active.acceptance_error = f"Coordinator acceptance failed: {exc}"
                if active.process.is_alive():
                    active.process.terminate()
                outcome = self._finish(event.invocation, active.acceptance_error)
                if outcome:
                    outcomes.append(outcome)
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
        for name, active in self.active.items():
            if active.process.is_alive() and not active.finished:
                active.process.terminate()
            active.process.join(timeout=5)
            if active.process.is_alive():
                active.process.kill()
                active.process.join()
            if not active.finished:
                self.acceptance.finalize(PipelineInvocationId(self.run_id, name), "Worker cancelled.")
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
