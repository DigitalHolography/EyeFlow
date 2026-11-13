"""One source-I/O lane; immutable commands in, completion records out.

Only this thread opens native dataset handles or writes dispatched destinations.
It never mutates readiness, pending consumers, priorities, leases or reservations.
"""

from dataclasses import dataclass
from queue import Queue
from threading import Thread, get_ident
from contextvars import copy_context
from time import perf_counter

import numpy as np

from utils.logger import Logger


@dataclass(frozen=True)
class SourceJob:
    token: int
    key: tuple
    selection: object = None
    destination: object = None
    flush: object = None
    copy_from: object = None
    priority: int = 0
    queued: float = 0.


@dataclass(frozen=True)
class SourceCompletion:
    token: int
    value: object
    error: str | None
    reads: tuple
    started: float
    ended: float
    flush_seconds: float


class SourceReader:
    def __init__(self, sources, read, lookup, wake):
        self.sources, self.read, self.lookup, self.wake = sources, read, lookup, wake
        self.commands, self.completed = Queue(maxsize=1), Queue(maxsize=1)
        self.datasets = {}
        self.current = None
        self._reads = []
        self._identity = None
        self.thread = Thread(target=copy_context().run, args=(self._run,), name="EyeFlowSourceIO", daemon=True)

    def dataset(self, key):
        if get_ident() != self._identity:
            raise RuntimeError("Native input datasets belong to the source reader")
        if key not in self.datasets:
            self.datasets[key] = self.sources[key[0]][key[1]]
        return self.datasets[key]

    def read_direct(self, dataset, destination, selection):
        started = perf_counter()
        dataset.read_direct(destination, source_sel=selection)
        self._record(dataset, selection, destination.nbytes, started, perf_counter(), False)

    def _record(self, dataset, selection, size, started, ended, small):
        identity = f"{dataset.file.filename}:{dataset.name}"
        self._reads.append((identity, size, started, ended, small))
        command = self.current
        Logger.log_perf(f"input_source_hdf5_read dataset={dataset.name} region={selection} "
            f"priority={command.priority} reason={('demand', 'lookahead', 'active_prefetch', 'future_prefetch')[command.priority]} "
            f"enqueued={command.queued:.6f} at={started:.6f} end={ended:.6f} "
            f"read_seconds={ended-started:.6f} logical_bytes={size} job={command.token}")

    def _execute(self, command):
        value = None
        if command.copy_from is not None:
            np.copyto(command.destination, command.copy_from)
        elif command.destination is None:
            dataset = self.lookup(command.key)
            started = perf_counter()
            value = np.asarray(dataset[()])
            self._record(dataset, (), value.nbytes, started, perf_counter(), True)
        elif command.destination.flags.c_contiguous:
            self.read(self.lookup(command.key), command.destination, command.selection)
        else:
            tile = np.empty(command.destination.shape, dtype=command.destination.dtype)
            self.read(self.lookup(command.key), tile, command.selection)
            command.destination[...] = tile
        flushed = perf_counter()
        if command.flush is not None:
            command.flush.flush()
        return value, perf_counter() - flushed if command.flush is not None else 0.

    def _run(self):
        self._identity = get_ident()
        try:
            while True:
                command = self.commands.get()
                if command is None:
                    return
                self.current, self._reads = command, []
                started, value, flush_seconds, error = perf_counter(), None, 0., None
                try:
                    value, flush_seconds = self._execute(command)
                except BaseException as exc:
                    error = f"{type(exc).__name__}: {exc}"
                completion = SourceCompletion(command.token, value, error, tuple(self._reads),
                                              started, perf_counter(), flush_seconds)
                # Drop all borrowed ndarray/buffer views BEFORE acknowledging
                # completion: the coordinator may immediately retire the storage.
                command = self.current = None
                self._reads = []
                self.completed.put(completion)
                self.wake()
        finally:
            for dataset in self.datasets.values():
                dataset.id.close()
            self.datasets.clear()

    def close(self):
        if self.thread.is_alive():
            self.commands.put(None)
        self.thread.join()
