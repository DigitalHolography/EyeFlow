"""Bounded delivery of accepted records to one independently admitted writer."""

from collections import deque
from contextlib import nullcontext
from contextvars import copy_context
from hashlib import sha256
from math import prod
import multiprocessing
from pathlib import Path
import pickle
from queue import Empty
from threading import Condition, Event, Thread
from time import perf_counter

from pipeline_engine.resources import suspend_compute
from scratch_location import is_local_path
from utils.logger import Logger

from .protocol import OutputPolicy, OutputServiceError, WriteCommand, WriterConnection
from .writer import writer_main


class OutputService:
    """Materialize accepted datasets; logical acceptance remains with OutputModel.

    The coordinator writes small manifests, never queues numerical array bytes.
    A dispatcher owns command admission and acknowledgements, so a full backlog
    does not depend on the scientific coordinator draining a completion channel.
    """

    def __init__(self, model, *, budget=None, context=None, policy=None, _writer_target=None):
        self.model, self.metrics = model, model.storage.metrics
        self.budget = budget if budget is not None else model.storage.budget
        self.policy = policy or OutputPolicy()
        context = context or multiprocessing.get_context("spawn")
        if context.get_start_method() != 'spawn':
            raise ValueError('The output writer requires a multiprocessing spawn context')
        # Explicit ProductStorage directories and retained worker directories
        # bypass the default scratch selector. Verify the actual writer volume
        # before allocating IPC resources or starting the child.
        if not is_local_path(model.storage.directory):
            raise OSError(f"EyeFlow output scratch must be local: {model.storage.directory}")
        self.local_path = model.storage.path(".h5")
        self._output_name = Path(model.filename).name
        receiver, self.commands = context.Pipe(duplex=False)
        self.acknowledgements, sender = context.Pipe(duplex=False)
        self._condition, self._stop = Condition(), Event()
        self._ready = Event()
        self._pending = deque()
        self._jobs = self._bytes = self._sequence = self._completed = 0
        self._manifest_bytes = 0
        self._cursor, self._submitted, self._written = {}, set(), set()
        self._accepted_times = {}
        self._manifests = set()
        self._error = None
        self._closed = self._finished = False
        self.process = context.Process(target=_writer_target or writer_main,
            args=(WriterConnection(receiver), WriterConnection(sender), str(self.local_path), self.policy.hdf5_cache_bytes),
            name="EyeFlowOutputWriter")
        self._thread = Thread(target=copy_context().run, args=(self._dispatch,), name="EyeFlowOutputDelivery", daemon=True)
        started = perf_counter()
        self._started = started
        Logger.log_perf(f"START pipeline=run variant=shared stage=hdf5_serialization at={started:.6f} progressive=True output={self._output_name}")
        try:
            # Claim the single materialization listener before spawning. Failed
            # duplicate construction must leave the existing service attached.
            model.bind_commit_listener(self.submit)
            self.process.start()
            self._thread.start()
        except BaseException:
            self.close()
            raise
        finally:
            receiver.close()
            sender.close()
        self.metrics.add("output_process_start_seconds", perf_counter()-started)

    @property
    def pending_count(self):
        with self._condition:
            return self._jobs

    @property
    def completed_sequence(self):
        with self._condition:
            return self._completed

    @property
    def written_paths(self):
        with self._condition:
            return frozenset(self._written)

    def check_health(self):
        error = self._error or self.model.commit_listener_error
        if error is not None:
            if not isinstance(error, Exception):
                # Acceptance retains interrupts until a safe coordinator boundary;
                # they must still stop the batch rather than become run failures.
                raise error
            raise OutputServiceError(f"Output serialization failed: {error}") from error
        if self._closed and not self._finished:
            raise OutputServiceError("Output serialization was aborted")

    def _fail(self, error):
        with self._condition:
            if self._error is None:
                self._error = error
                self.metrics.add("output_failures", 1)
            self._condition.notify_all()

    def request_abort(self):
        """Wake delivery/finalization without deleting storage still in use."""
        self._fail(OutputServiceError("Output serialization was cancelled"))
        self._stop.set()
        with self._condition:
            self._condition.notify_all()

    @staticmethod
    def _size(record):
        value = record.value
        shape, dtype = getattr(value, "shape", None), getattr(value, "dtype", None)
        return prod(shape)*dtype.itemsize if shape is not None and dtype is not None else 0

    def submit(self, model):
        """Post-commit observer: failure must never unwind scientific acceptance."""
        if self._closed or self._finished or self._error is not None:
            return
        try:
            records, cursor = model.committed_delta(self._cursor)
            observed_at = perf_counter()
            for path, record in records.items():
                if record.kind == 'dataset':
                    self._accepted_times.setdefault(path, observed_at)
            batch, size = [], 0
            for path, record in records.items():
                if record.kind != "dataset" or path in self._submitted:
                    continue
                amount = self._size(record)
                if batch and (len(batch) >= 32 or size+amount > self.policy.max_pending_bytes):
                    if self._enqueue(batch, size, block=False) is None:
                        return
                    self._submitted.update(path for path, _ in batch)
                    batch, size = [], 0
                batch.append((path, record))
                size += amount
            if batch:
                if self._enqueue(batch, size, block=False) is None:
                    return
                self._submitted.update(path for path, _ in batch)
            self._cursor = cursor
        except Exception as exc:
            self._fail(exc)

    def _enqueue(self, contents, size=0, *, final=False, block=True):
        started = perf_counter()
        # Yield an enclosing scientific reservation while applying backpressure.
        # This matters on a single CPU: delivery must be able to admit the writer.
        with (suspend_compute() if block else nullcontext()), self._condition:
            while self._jobs >= self.policy.max_pending_jobs or (
                    self._jobs and self._bytes+size > self.policy.max_pending_bytes):
                self.check_health()
                if not block:
                    # Accepted records remain in the authoritative model ledger.
                    # A coordinator must remain free to reap scientific workers
                    # whose CPU reservations may currently prevent writer admission.
                    self.metrics.add("output_submission_deferred", 1)
                    return None
                if self._stop.is_set():
                    raise OutputServiceError("Output serialization was aborted")
                self._condition.wait(0.1)
            self.check_health()
            path = self.model.storage.path(".pickle")
            try:
                with path.open("wb") as stream:
                    pickle.dump(contents, stream, protocol=5)
                digest = sha256()
                with path.open("rb") as stream:
                    for block in iter(lambda: stream.read(2**20), b""):
                        digest.update(block)
            except BaseException:
                path.unlink(missing_ok=True)
                raise
            self._sequence += 1
            accepted_at = started if final else min(self._accepted_times[path] for path, _ in contents)
            command = WriteCommand(self._sequence, str(path), digest.hexdigest(), final, accepted_at)
            self._manifests.add(path)
            self._manifest_bytes += path.stat().st_size
            self._pending.append((command, size, accepted_at))
            self._jobs += 1
            self._bytes += size
            self.metrics.maximum("output_peak_pending_jobs", self._jobs)
            self.metrics.maximum("output_peak_pending_logical_bytes", self._bytes)
            self.metrics.maximum("output_peak_pending_manifest_bytes", self._manifest_bytes)
            self.metrics.add("output_command_manifest_bytes", path.stat().st_size)
            self.metrics.add("output_submission_seconds", perf_counter()-started)
            if not final:
                self.metrics.add("output_datasets_accepted", len(contents))
                self.metrics.add("output_logical_bytes_accepted", size)
            Logger.log_perf(f"output_commit output={self._output_name} sequence={command.sequence} datasets={0 if final else len(contents)} accepted_at={accepted_at:.6f} submitted_at={perf_counter():.6f}")
            self._condition.notify_all()
            return command.sequence

    def _dispatch(self):
        try:
            while not self._stop.is_set():
                try:
                    if not self.acknowledgements.poll(0.1):
                        raise Empty
                    startup = self.acknowledgements.recv()
                    break
                except Empty:
                    if not self.process.is_alive():
                        raise OutputServiceError(f"Output writer exited during startup with code {self.process.exitcode}") from None
                    if perf_counter()-self._started > self.policy.operation_timeout:
                        raise OutputServiceError("Output writer startup timed out") from None
            else:
                return
            if startup.error or startup.sequence != 0:
                raise OutputServiceError(startup.error or "Output writer startup acknowledgement differs")
            self.metrics.add("output_writer_startup_seconds", perf_counter()-self._started)
            self._ready.set()
            while not self._stop.is_set():
                with self._condition:
                    if not self._pending:
                        self._condition.wait(0.1)
                        if not self.process.is_alive() and not self._finished and not self._stop.is_set():
                            raise OutputServiceError("Output writer exited unexpectedly")
                        continue
                    command, size, committed = self._pending.popleft()
                reserved = False
                # Numerical payload stays mapped/copied by HDF5, while the small
                # manifest must be decoded into Python objects in the child.
                memory = (min(self.policy.admission_memory_bytes, self.budget.memory_limit//8)
                          + 2*Path(command.manifest).stat().st_size) if self.budget else 0
                admitted_at = perf_counter()
                try:
                    if self.budget:
                        while not self._stop.is_set():
                            if not self.process.is_alive():
                                raise OutputServiceError("Output writer exited before admission")
                            if self.budget.try_acquire(1, memory):
                                reserved = True
                                break
                            if perf_counter()-admitted_at > self.policy.operation_timeout:
                                raise OutputServiceError("Output writer resource admission timed out")
                            self._stop.wait(0.05)
                    if self._stop.is_set():
                        return
                    self.metrics.add("output_admission_wait_seconds", perf_counter()-admitted_at)
                    sent = perf_counter()
                    self.commands.send(command)
                    while True:
                        if self._stop.is_set():
                            return
                        try:
                            if not self.acknowledgements.poll(0.1):
                                raise Empty
                            ack = self.acknowledgements.recv()
                            break
                        except Empty:
                            if not self.process.is_alive():
                                raise OutputServiceError(f"Output writer exited with code {self.process.exitcode}") from None
                            if perf_counter()-sent > self.policy.operation_timeout:
                                raise OutputServiceError("Output writer exceeded its operation timeout") from None
                    if ack.sequence != command.sequence:
                        raise OutputServiceError("Output acknowledgement sequence differs")
                    if ack.error:
                        raise OutputServiceError(ack.error)
                    for name, value in ack.metrics.items():
                        if name == "output_writer_sampled_rss_bytes":
                            self.metrics.maximum(name, value)
                        else:
                            self.metrics.add(name, value)
                    for message in ack.messages:
                        Logger.log_perf(f"output={self._output_name} {message}")
                    self.metrics.add("output_commit_to_write_complete_seconds", perf_counter()-committed)
                    self.metrics.add("output_commands_completed", 1)
                    with self._condition:
                        self._jobs -= 1
                        self._bytes -= size
                        self._completed = command.sequence
                        self._written.update(ack.paths)
                        self._condition.notify_all()
                    manifest = Path(command.manifest)
                    with self._condition:
                        self._manifest_bytes -= manifest.stat().st_size
                    manifest.unlink(missing_ok=True)
                    self._manifests.discard(manifest)
                    if command.final:
                        return
                finally:
                    if reserved:
                        self.budget.release(1, memory)
        except BaseException as exc:
            if not self._stop.is_set():
                self._fail(OutputServiceError("Output writer exited or closed its IPC channel")
                           if isinstance(exc, (EOFError, BrokenPipeError)) else exc)

    def drain(self):
        """Wait for already submitted writes without finalizing mutable metadata."""
        with suspend_compute(), self._condition:
            while self._jobs:
                self.check_health()
                self._condition.wait(0.1)
            self.check_health()

    def finish(self):
        if self._finished:
            return self.local_path
        if self.model._checkpoints:
            raise RuntimeError("Cannot finalize an active output transaction")
        started = perf_counter()
        self.metrics.add("output_jobs_at_analysis_end", self.pending_count)
        records = self.model.committed_records()
        groups = {path: record for path, record in records.items() if record.kind == "group"}
        expected = tuple(path for path, record in records.items() if record.kind == "dataset")
        with suspend_compute():
            while not set(expected).issubset(self._submitted):
                self.submit(self.model)
                self.check_health()
                self.drain()
        self.check_health()
        self._enqueue((groups, expected), final=True)
        self.drain()
        self.process.join(self.policy.shutdown_timeout)
        if self.process.is_alive() or self.process.exitcode:
            self._fail(OutputServiceError("Output writer did not close successfully"))
        self.check_health()
        self._finished = True
        self.metrics.add("output_final_drain_seconds", perf_counter()-started)
        Logger.log_perf(f"END pipeline=run variant=shared stage=hdf5_serialization at={perf_counter():.6f} progressive=True wall={perf_counter()-self._started:.3f}s output={self._output_name}")
        Logger.log_perf(f"output_summary metrics={self.metrics.snapshot()}")
        return self.local_path

    def close(self):
        if self._closed:
            return
        self._closed = True
        self.model.detach_commit_listener(self.submit)
        self._stop.set()
        with self._condition:
            self._condition.notify_all()
        if self._thread.ident is not None:
            self._thread.join(self.policy.shutdown_timeout)
        if self.process.pid is not None:
            self.process.join(0 if not self._finished else self.policy.shutdown_timeout)
            if self.process.is_alive():
                self.process.terminate()
                self.process.join(self.policy.shutdown_timeout)
                if self.process.is_alive():
                    self.process.kill()
                    self.process.join()
        if self._thread.ident is not None:
            self._thread.join(self.policy.shutdown_timeout)
        self.commands.close()
        self.acknowledgements.close()
        if self.process.pid is not None:
            self.process.close()
        self.commands = self.acknowledgements = None
        for manifest in tuple(self._manifests):
            manifest.unlink(missing_ok=True)
        self._manifests.clear()
        with self._condition:
            self._pending.clear()
            self._jobs = self._bytes = self._manifest_bytes = 0
            self._condition.notify_all()
        self.local_path.unlink(missing_ok=True)

    abort = close
