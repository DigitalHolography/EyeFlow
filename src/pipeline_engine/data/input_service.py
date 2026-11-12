"""One coordinator-owned NAS reader with reusable local read-only datasets."""

from __future__ import annotations

from dataclasses import dataclass, field
from collections import deque
from collections.abc import Mapping
from contextvars import copy_context
from threading import RLock, Thread, Event
from time import perf_counter
from queue import Empty, Full
from multiprocessing.shared_memory import SharedMemory
from pathlib import Path
import shutil

import numpy as np

from input_output.readers.nodes import is_dataset

from pipeline_engine.data.storage import ArrayReference, ValueReference
from pipeline_engine.stage_graph.model import InputAccess, StageGraph
from utils.logger import Logger


@dataclass
class SourceEntry:
    attrs: dict
    shape: tuple | None = None
    dtype: object = None
    identity: str | None = None
    chunks: tuple | None = None
    compression: str | None = None


class InputClient:
    def __init__(self, identity, requests, responses, stopped):
        self.identity, self.requests, self.responses = identity, requests, responses
        self.stopped = stopped
        self.cache = {}
        self._lock = RLock()
        self._streams = set()

    def __getstate__(self):
        if self._streams:
            raise RuntimeError("Cannot transfer an input client with active leases")
        # Spawn transports controls, never previously opened mappings/arrays.
        return {**self.__dict__, "_lock": None, "cache": {}, "_streams": set()}

    def __setstate__(self, state):
        self.__dict__.update(state)
        self._lock = RLock()

    def rpc(self, operation, *args):
        # A client has one response queue. Serialize threads as well as requests
        # so concurrent reads cannot consume each other's reply.
        with self._lock:
            started = perf_counter()
            if self.stopped.is_set():
                raise OSError("Input service has stopped")
            self.send(operation, *args)
            while True:
                try:
                    value, error = self.responses.get(timeout=0.1)
                    break
                except Empty:
                    if self.stopped.is_set():
                        raise OSError("Input service stopped while waiting for a response") from None
            if operation in ("read", "stream_next"):
                self.send("input_wait", started, perf_counter())
            if error:
                raise OSError(error)
            return value

    def send(self, operation, *args):
        while not self.stopped.is_set():
            try:
                self.requests.put((self.identity, operation, args, perf_counter()), timeout=0.1)
                return
            except Full:
                pass

    def read(self, slot, path):
        with self._lock:
            return self._read(slot, path)

    def _read(self, slot, path):
        key = (slot, path)
        if key not in self.cache:
            value = self.rpc("read", slot, path, None)
            array = value.open() if isinstance(value, (ArrayReference, ValueReference)) else value
            if isinstance(array, np.ndarray):
                array.flags.writeable = False
            self.cache[key] = array
        return self.cache[key]

    def read_slice(self, slot, path, shape, selection):
        with self._lock:
            return self._read_slice(slot, path, shape, selection)

    def _read_slice(self, slot, path, shape, selection):
        span, local = _selection_span(shape, selection)
        if span is None:
            return self.read(slot, path)[selection]
        first, stop = span
        if first == stop:
            value = np.empty((0, *shape[1:]), dtype=self.rpc("dtype", slot, path))
            value.flags.writeable = False
            return value[local]
        reference = self.rpc("read", slot, path, span)
        if not isinstance(reference, _RegionReference):
            key = (slot, path)
            if key not in self.cache:
                value = reference.open() if isinstance(reference, (ArrayReference, ValueReference)) else reference
                if isinstance(value, np.ndarray):
                    value.flags.writeable = False
                self.cache[key] = value
            return self.cache[key][selection]
        key = (reference.path, first, stop)
        if key not in self.cache:
            self.cache[key] = reference.open()
        return self.cache[key][local]

    def close(self):
        for stream in list(self._streams):
            stream.close()
        self.send("client_close")
        for value in self.cache.values():
            if isinstance(value, np.memmap) and not value._mmap.closed:
                value._mmap.close()
        self.cache.clear()


class InputNode:
    def __init__(self, file, path):
        self.file, self.path = file, path

    @property
    def name(self):
        return "/" + self.path

    @property
    def attrs(self):
        return dict(self.file.entries[self.path].attrs)

    @property
    def parent(self):
        return self.file[self.path.rpartition("/")[0]]


class InputDataset(InputNode):
    _eyeflow_dataset = True

    @property
    def id(self):
        # Stable within the run, including separate proxy objects/workflows.
        return self.file.slot, self.file.filename, self.file.entries[self.path].identity

    @property
    def shape(self):
        return self.file.entries[self.path].shape

    @property
    def dtype(self):
        return self.file.entries[self.path].dtype

    @property
    def chunks(self):
        return self.file.entries[self.path].chunks

    @property
    def ndim(self):
        return len(self.shape)

    @property
    def size(self):
        return int(np.prod(self.shape))

    @property
    def nbytes(self):
        return self.size * self.dtype.itemsize

    @property
    def _value(self):
        identity = self.file.entries[self.path].identity
        return self.file.client.read(self.file.slot, identity)

    def __getitem__(self, selection):
        return self.read_slice(selection)

    def read_slice(self, selection):
        return self.file.client.read_slice(self.file.slot, self.path, self.shape, selection)

    def read_direct(self, dest, source_sel=None, dest_sel=None):
        np.copyto(dest if dest_sel is None else dest[dest_sel],
                  self.read_slice(() if source_sel is None else source_sel), casting="unsafe")

    def iter_blocks(self, *, frames=None, writable=False, max_bytes=None):
        """Use as a context manager; each yielded lease owns one exclusive block."""
        return _BlockStream(self, frames, writable, max_bytes)

    def __array__(self, dtype=None, copy=None):
        array = np.asarray(self._value, dtype=dtype)
        return array.copy() if copy else array

    def astype(self, dtype):
        return np.asarray(self._value, dtype=dtype)


class InputFile(InputNode):
    _eyeflow_group = True

    def __init__(self, slot, filename, entries, client):
        super().__init__(self, "")
        self.slot, self.filename, self.entries, self.client = slot, filename, entries, client

    def get(self, path, default=None):
        path = str(path).strip("/")
        entry = self.entries.get(path)
        if entry is None:
            return default
        return (InputDataset if entry.shape is not None else InputGroup)(self, path)

    def __getitem__(self, path):
        found = self.get(path)
        if found is None:
            raise KeyError(path)
        return found

    def __contains__(self, path):
        return str(path).strip("/") in self.entries

    def keys(self):
        return InputGroup.keys(self)

    def items(self):
        return [(key, self[key]) for key in self.keys()]


class InputGroup(InputNode):
    _eyeflow_group = True

    def get(self, path, default=None):
        path = str(path)
        if path.startswith("/"):
            return self.file.get(path, default)
        return self.file.get("/".join((self.path, path.strip("/"))), default)

    def __getitem__(self, path):
        value = self.get(path)
        if value is None:
            raise KeyError(path)
        return value

    def keys(self):
        prefix = self.path + "/" if self.path else ""
        return sorted(
            {
                path[len(prefix) :].split("/")[0]
                for path in self.file.entries
                if path.startswith(prefix) and path != self.path
            }
        )

    def items(self):
        return [(key, self[key]) for key in self.keys()]


def _selection_span(shape, selection):
    """Bound basic first-axis indexing; advanced indexing uses full replay.

    Returning a mapping of only this extent prevents a caller's ndarray.base
    from exposing unfinished portions of the replay file.
    """
    items = selection if isinstance(selection, tuple) else (selection,)
    if not shape or not items:
        return None, selection
    first, *rest = items
    if any(not (isinstance(item, (slice, int, np.integer)) or item is Ellipsis or item is None)
           or isinstance(item, (bool, np.bool_)) for item in items):
        return None, selection
    if isinstance(first, (int, np.integer)):
        index = int(first)
        if not -shape[0] <= index < shape[0]:
            raise IndexError("Input frame index out of range")
        index %= shape[0]
        return (index, index + 1), (0, *rest)
    if isinstance(first, slice):
        start, stop, step = first.indices(shape[0])
        if step == 1:
            stop = max(start, stop)
            return (start, stop), (slice(None), *rest)
    return None, selection


def _union_intervals(intervals):
    merged = []
    for start, stop in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], stop))
        else:
            merged.append((start, stop))
    return merged


def _subtract_intervals(interval, excluded):
    first, stop = interval
    for left, right in _union_intervals(excluded):
        if first < min(stop, left):
            yield first, min(stop, left)
        first = max(first, right)
        if first >= stop:
            return
    if first < stop:
        yield first, stop


def interval_overlap(left, right):
    """Wall-clock intersection of unions; never double count concurrent spans."""
    a, b = _union_intervals(left), _union_intervals(right)
    i = j = 0
    total = 0.0
    while i < len(a) and j < len(b):
        total += max(0, min(a[i][1], b[j][1]) - max(a[i][0], b[j][0]))
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return total


@dataclass(frozen=True)
class InputPolicy:
    """All provisional I/O and capacity defaults, independently benchmarkable."""

    replay_frames: int = 32
    stream_frames: int = 16
    stream_buffers: int = 2
    max_io_bytes: int = 64 * 2**20
    max_streams: int = 2
    max_prefetch_jobs: int = 16
    lookahead: bool = True

    def __post_init__(self):
        for name in ("replay_frames", "stream_frames", "stream_buffers", "max_io_bytes", "max_streams",
                     "max_prefetch_jobs"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 1:
                raise ValueError(f"{name} must be a positive integer")


@dataclass(frozen=True)
class _RegionReference:
    path: str
    shape: tuple
    dtype: object
    offset: int

    def open(self):
        return np.memmap(self.path, dtype=self.dtype, mode="r", offset=self.offset, shape=self.shape)


@dataclass
class _Replay:
    array: np.memmap
    reference: ArrayReference
    frames: int
    states: list[str]
    offset: int
    created: float = field(default_factory=perf_counter)
    error: str | None = None
    demanded: set[int] = field(default_factory=set)
    served: set[int] = field(default_factory=set)


class _BlockLease:
    """Exclusive buffer ownership lasts until release; views must not escape it."""

    def __init__(self, stream, info):
        self.stream = stream
        self.index, name, shape, dtype, self.start, self.stop = info
        self.memory = SharedMemory(name=name)
        self.block = np.ndarray(shape, dtype=dtype, buffer=self.memory.buf)
        self.block.flags.writeable = stream.writable
        self.released = False

    def __enter__(self):
        if self.released:
            raise RuntimeError("Input lease has already been released")
        return self.block

    def __exit__(self, *exc):
        self.release()

    def release(self):
        if self.released:
            return
        self.released = True
        self.block = None
        self.memory.close()
        self.stream.client.send("stream_release", self.stream.identity, self.index)


class _BlockStream:
    def __init__(self, dataset, frames, writable, max_bytes):
        self.client = dataset.file.client
        self.identity = self.client.rpc("stream_open", dataset.file.slot, dataset.path, frames, max_bytes)
        self.writable = writable
        self.lease = None
        self.closed = False
        self.client._streams.add(self)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def __iter__(self):
        return self

    def __next__(self):
        if self.closed:
            raise StopIteration
        if self.lease is not None and not self.lease.released:
            raise RuntimeError("Release the current input lease before advancing")
        try:
            info = self.client.rpc("stream_next", self.identity)
            if info is None:
                self.close()
                raise StopIteration
            self.lease = _BlockLease(self, info)
            return self.lease
        except Exception:
            self.close()
            raise

    def close(self):
        if self.closed:
            return
        self.closed = True
        if self.lease is not None:
            self.lease.release()
        self.client.send("stream_close", self.identity)
        self.client._streams.discard(self)


class InputService:
    """Authoritative metadata, bounded source reads and input-buffer ownership.

    Only the reader thread mutates job, replay and stream state. Each source
    operation yields to request traffic at a slab boundary. Scientific products
    never enter this cache. Source files must outlive close().
    """

    def __init__(self, sources, storage, mp_context, *, policy=None):
        self.sources, self.storage = sources, storage
        self.policy = policy or InputPolicy()
        self.requests = mp_context.Queue(maxsize=32)
        self.responses = {}
        self._context = mp_context
        self._lock = RLock()
        self._closed = False
        self.cache, self.metadata, self.source_reads, self.aliases = {}, {}, {}, {}
        self.manifest, self.unresolved, self.pairs = {}, [], {}
        self._configured = False
        self._undeclared = set()
        self._replays, self._jobs, self._pending, self._streams = {}, {}, {}, {}
        self._datasets = {}
        self._serial = 0
        self._stream_serial = 0
        self._cancelled_clients = set()
        self._stop = Event()
        self.stopped = mp_context.Event()
        # Bound diagnostic memory even for unusually repetitive consumers.
        self.io_intervals = deque(maxlen=8192)
        self.wait_intervals = deque(maxlen=8192)
        self.compute_intervals = deque(maxlen=8192)
        self._intervals_dropped = False
        started = perf_counter()
        for slot, source in sources.items():
            if source is None:
                continue
            entries = {"": SourceEntry(dict(source.attrs))}
            identities = {}

            def walk(group, prefix="", ancestors=(), *, entries=entries,
                     identities=identities, slot=slot):
                for name in group:
                    node = group.get(name)
                    if node is None:
                        continue
                    path = "/".join(part for part in (prefix, name) if part)
                    identity = identities.setdefault(node.id, path)
                    dataset = is_dataset(node)
                    entries[path] = SourceEntry(
                        dict(node.attrs), tuple(node.shape) if dataset else None,
                        node.dtype if dataset else None, identity,
                        getattr(node, "chunks", None), getattr(node, "compression", None),
                    )
                    if dataset:
                        self.aliases[slot, path] = slot, identity
                    elif node.id not in ancestors:
                        walk(node, path, (*ancestors, node.id))

            walk(source, ancestors=(source.id,))
            self.metadata[slot] = (source.filename, entries)
        self._metric("input_metadata_seconds", perf_counter() - started)
        self.thread = Thread(target=copy_context().run, args=(self._serve,), name="EyeFlowInputs", daemon=True)
        self.thread.start()

    def _metric(self, name, value=1):
        self.storage.metrics.add(name, value)

    def _record_interval(self, intervals, interval):
        if len(intervals) == intervals.maxlen:
            self._intervals_dropped = True
            self._metric("input_interval_samples_dropped")
        intervals.append(interval)

    def configure(self, graphs: Mapping[str, StageGraph | None],
                  run_settings: Mapping[str, object] | None = None) -> None:
        """Resolve graph declarations without running callbacks or reading videos.

        Invalid/missing candidates are recorded, leaving branch validation and
        failure isolation with the scientific stages. No optional candidate is
        speculatively read merely because it appears here.
        """
        from input_output.schema.holodoppler import HolodopplerSource
        from input_output.schema.doppler_view import DopplerViewSource
        from input_output.readers.h5 import RawH5SourceReader
        from app_settings import selected_velocity_methods

        started = perf_counter()
        settings = run_settings or {}
        selection = settings.get("velocity_workflow_selection", "both")
        methods = set(selected_velocity_methods(selection))
        adapters = {slot: cls(RawH5SourceReader(h5file=self.sources[slot], label=slot.upper())) for slot, cls in
                    (("hd", HolodopplerSource), ("dv", DopplerViewSource))
                    if self.sources.get(slot) is not None}
        semantic = {}
        for pipeline, graph in graphs.items():
            if graph is None:
                continue
            for stage in graph.stages:
                for use in stage.inputs:
                    if use.when is not None and use.when not in methods:
                        continue
                    slot, name = use.source.split(".", 1)
                    try:
                        nodes = adapters[slot].input_datasets(name)
                        if not nodes:
                            raise KeyError(f"Missing {use.source}")
                    except (KeyError, ValueError, TypeError) as exc:
                        self.unresolved.append((pipeline, stage.id, use, str(exc)))
                        continue
                    for node in nodes:
                        key = self.aliases[slot, node.name.strip("/")]
                        semantic[use.source] = key
                        entry = self.metadata[slot][1][key[1]]
                        item = self.manifest.setdefault(key, {
                            "entry": entry, "consumers": [], "access": set(), "strategy": None})
                        item["consumers"].append((pipeline, stage.id, use))
                        item["access"].add(use.access)
        for item in self.manifest.values():
            access = item["access"]
            item["strategy"] = ("metadata" if access == {InputAccess.METADATA} else
                                "stream" if access <= {InputAccess.METADATA, InputAccess.STREAM} else
                                "small" if access <= {InputAccess.METADATA, InputAccess.SMALL} else "replay")
            if sum(use.access == InputAccess.STREAM for _, _, use in item["consumers"]) > 1:
                item["strategy"] = "replay"
        for left, right in (("hd.moment0", "hd.moment2"), ("hd.band_low", "hd.band_high")):
            if all(source in semantic and any(
                not use.optional and use.access == InputAccess.SEQUENTIAL
                for _, _, use in self.manifest[semantic[source]]["consumers"])
                   for source in (left, right)):
                self.pairs[semantic[left]] = semantic[right]
                self.pairs[semantic[right]] = semantic[left]
        self._metric("input_manifest_seconds", perf_counter() - started)
        self._configured = True
        Logger.log_perf(f"input_manifest datasets={len(self.manifest)} unresolved={len(self.unresolved)} "
                        f"strategies={[(key, item['strategy']) for key, item in self.manifest.items()]}")

    def activate(self, pipeline: str) -> None:
        """Called on actual dispatch, without changing scientific scheduling."""
        self.requests.put((None, "activate", (pipeline,), perf_counter()))

    def prefetch(self, key: tuple[str, str], *, first: int = 0) -> None:
        self.requests.put((None, "prefetch", (key, first), perf_counter()))

    def cancel_prefetch(self, key: tuple[str, str] | None = None) -> None:
        """Cancel queued speculation, preserving urgent demand and ready regions."""
        self.requests.put((None, "cancel_prefetch", (key,), perf_counter()))

    def release_client(self, identity: int) -> None:
        """Reclaim streams after a worker has exited, including hard crashes."""
        self.requests.put((None, "reclaim", (identity,), perf_counter()))

    def client(self) -> InputClient:
        with self._lock:
            identity = len(self.responses)
            response = self._context.Queue(maxsize=1)
            self.responses[identity] = response
        return InputClient(identity, self.requests, response, self.stopped)

    def files(self, client: InputClient) -> dict[str, InputFile | None]:
        return {slot: InputFile(slot, *self.metadata[slot], client) if slot in self.metadata else None
                for slot in self.sources}

    def _reply(self, identity, value=None, error=None):
        if identity is None or identity in self._cancelled_clients:
            return
        try:
            self.responses[identity].put_nowait((value, error))
        except Full:
            # An exited client may leave one unread reply; never block the
            # authoritative reader on an abandoned response queue.
            Logger.log_debug(f"Input client {identity} left an unread response")

    def _enqueue(self, key, priority, operation, *args):
        previous = self._jobs.get(key)
        if previous is not None and previous[0] <= priority:
            return True
        if previous is not None:
            self._metric("input_priority_changes")
        elif priority > 0 and operation == "replay":
            speculative = [(k, value) for k, value in self._jobs.items()
                           if value[0] > 0 and value[2] == "replay"]
            if len(speculative) >= self.policy.max_prefetch_jobs:
                worst, value = max(speculative, key=lambda item: item[1][:2])
                if value[0] <= priority:
                    return False
                self._discard_job(worst)
        self._serial += 1
        self._jobs[key] = (priority, self._serial, operation, args, perf_counter())
        return True

    def _discard_job(self, job_key):
        self._jobs.pop(job_key, None)
        if job_key[0] == "replay":
            _, key, region = job_key
            self._replays[key].states[region] = "not_requested"
            self._metric("input_speculative_jobs_cancelled")

    def _frames(self, dataset, frames):
        frame_bytes = max(1, int(np.prod(dataset.shape[1:])) * dataset.dtype.itemsize)
        return min(frames, max(1, self.policy.max_io_bytes // frame_bytes))

    def _replay(self, key):
        if key in self._replays:
            return self._replays[key]
        dataset = self._dataset(key)
        if dataset.dtype.hasobject or not dataset.shape or not dataset.size or dataset.nbytes < 4096:
            if key not in self.cache:
                self.cache[key] = self._stage_small(dataset)
                self.source_reads[key] = 1
            return None
        if shutil.disk_usage(self.storage.directory).free < dataset.nbytes + 4096:
            raise OSError("Insufficient local scratch disk space for input replay")
        path = self.storage.path(".npy")
        try:
            array = np.lib.format.open_memmap(path, mode="w+", dtype=dataset.dtype, shape=dataset.shape)
        except BaseException:
            path.unlink(missing_ok=True)
            raise
        frames = self._frames(dataset, dataset.shape[0] if self.manifest.get(key, {}).get("strategy") == "small"
                              else self.policy.replay_frames)
        replay = _Replay(array, ArrayReference(str(path), tuple(dataset.shape), dataset.dtype),
                         frames, ["not_requested"] * ((dataset.shape[0] + frames - 1) // frames), array.offset)
        self._replays[key] = replay
        return replay

    def _dataset(self, key):
        # Keep the native dataset handle (and its raw chunk cache) alive across
        # slabs. Reopening H5D for every block can defeat configured cache reuse.
        key = self.aliases[key]
        if key not in self._datasets:
            self._datasets[key] = self.sources[key[0]][key[1]]
        return self._datasets[key]

    def _queue_region(self, key, region, priority):
        replay = self._replay(key)
        if replay is None or replay.error or not 0 <= region < len(replay.states):
            return
        if replay.states[region] in ("ready", "loading"):
            return
        if self._enqueue(("replay", key, region), priority, "replay", key, region):
            replay.states[region] = "queued"

    def _lookahead(self, key, stop):
        if not self.policy.lookahead:
            return
        # Pair the current region and exactly one following window. Inactive
        # methods have no pair entries and are never pulled into this lookahead.
        windows = []
        for candidate in dict.fromkeys((key, self.pairs.get(key, key))):
            replay = self._replay(candidate)
            if replay is not None:
                region = max(0, (stop - 1) // replay.frames)
                windows.append((candidate, region))
                for job_key, job in list(self._jobs.items()):
                    if (job_key[0] == "replay" and job_key[1] == candidate and job[0] == 1
                            and job_key[2] not in (region, region + 1)):
                        self._discard_job(job_key)
        for offset in (0, 1):
            for candidate, region in windows:
                self._queue_region(candidate, region + offset, 1)

    def _handle(self, request):
        identity, operation, args, requested_at = request
        self._metric("input_queue_wait_seconds", max(0, perf_counter() - requested_at))
        try:
            if operation == "read":
                slot, path, span = args
                key = self.aliases[slot, path]
                if self._configured and key not in self.manifest and key not in self._undeclared:
                    Logger.log_debug(f"Undeclared source demand {slot}:{path}; using safe replay fallback")
                    self._undeclared.add(key)
                replay = self._replay(key)
                if replay is None:
                    self._reply(identity, self.cache[key])
                else:
                    first, stop = span if span is not None else (0, replay.array.shape[0])
                    regions = range(first // replay.frames, (stop - 1) // replay.frames + 1)
                    self._metric("input_prefetch_hits" if all(replay.states[r] == "ready" for r in regions)
                                 else "input_prefetch_misses")
                    replay.demanded.update(regions)
                    self._pending[identity] = (key, first, stop, span is None, requested_at)
                    # Only one missing slab per pending request is queued. Full
                    # materialization cannot create a volume-sized job queue.
                    self._fulfill()
                    self._lookahead(key, stop if span is not None else min(stop, replay.frames))
            elif operation == "dtype":
                slot, path = args
                self._reply(identity, self.metadata[slot][1][path].dtype)
            elif operation in ("input_wait", "compute"):
                intervals = self.wait_intervals if operation == "input_wait" else self.compute_intervals
                self._record_interval(intervals, (identity, *args))
                if operation == "input_wait":
                    self._metric("input_consumer_blocked_seconds", args[1] - args[0])
            elif operation == "shutdown":
                self._stop.set()
            elif operation == "cancel_prefetch":
                key = self.aliases[args[0]] if args[0] is not None else None
                for job_key, job in list(self._jobs.items()):
                    if job_key[0] == "replay" and job[0] > 0 and (key is None or job_key[1] == key):
                        self._discard_job(job_key)
            elif operation in ("activate", "prefetch"):
                if not self.policy.lookahead:
                    return
                if operation == "activate":
                    keys = [key for key, item in self.manifest.items() if item["strategy"] == "replay"
                            and any(p == args[0] and not use.optional
                                    for p, _, use in item["consumers"])]
                    first = 0
                else:
                    key, first = args
                    keys = [self.aliases[key]]
                for key in keys:
                    replay = self._replay(key)
                    if replay is not None:
                        self._queue_region(key, first // replay.frames, 2 if operation == "activate" else 3)
            elif operation == "stream_open":
                self._open_stream(identity, *args)
            elif operation == "stream_next":
                stream = self._owned_stream(identity, args[0])
                stream["waiting"] = requested_at
                self._stream_ready(args[0])
            elif operation == "stream_release":
                stream_id, index = args
                stream = self._owned_stream(identity, stream_id)
                buffer = stream["buffers"][index]
                if buffer["state"] != "in_use":
                    raise RuntimeError("Streaming buffer released without a lease")
                buffer["state"] = "free"
                self._fill_free(stream_id)
            elif operation == "stream_close":
                stream = self._streams.get(args[0])
                if stream is not None and stream["owner"] == identity:
                    self._close_stream(args[0])
            elif operation in ("client_close", "reclaim"):
                owner = identity if operation == "client_close" else args[0]
                self._cancelled_clients.add(owner)
                self._pending.pop(owner, None)
                for stream_id, stream in list(self._streams.items()):
                    if stream["owner"] == owner:
                        self._close_stream(stream_id)
        except Exception as exc:  # return source failures without killing unrelated consumers
            self._reply(identity, error=f"Input {operation} failed: {exc}")

    def _fulfill(self):
        for identity, (key, first, stop, full, started) in list(self._pending.items()):
            replay = self._replays[key]
            missing = next((r for r in range(first // replay.frames, (stop - 1) // replay.frames + 1)
                            if replay.states[r] != "ready"), None)
            if missing is not None:
                if replay.error:
                    self._reply(identity, error=replay.error)
                    del self._pending[identity]
                    continue
                self._queue_region(key, missing, 0)
                continue
            reference = replay.reference if full else _RegionReference(
                replay.reference.path, (stop - first, *replay.array.shape[1:]), replay.array.dtype,
                replay.offset + first * int(np.prod(replay.array.shape[1:])) * replay.array.dtype.itemsize)
            self._reply(identity, reference)
            self._metric("input_region_wait_seconds", perf_counter() - started)
            regions = range(first // replay.frames, (stop - 1) // replay.frames + 1)
            self._metric("input_region_reuses", sum(r in replay.served for r in regions))
            replay.served.update(regions)
            del self._pending[identity]

    def _read_direct(self, dataset, destination, selection):
        started = perf_counter()
        dataset.read_direct(destination, source_sel=selection)
        ended = perf_counter()
        self._record_interval(self.io_intervals, (started, ended))
        self._metric("input_hdf5_read_seconds", ended - started)
        self._metric("input_logical_bytes", destination.nbytes)
        self._metric("source_reads")
        if destination.nbytes >= 2**20:
            Logger.log_perf(
                f"input_read dataset={dataset.name} region={selection} "
                f"priority={getattr(self, '_priority', 0)} "
                f"queue_seconds={getattr(self, '_job_wait', 0):.6f} "
                f"at={started:.6f} end={ended:.6f} read_seconds={ended - started:.6f} "
                f"logical_bytes={destination.nbytes}"
            )

    def _run_job(self, job):
        priority, _, operation, args, queued = job
        self._priority, self._job_wait = priority, perf_counter() - queued
        self._metric("input_job_queue_seconds", self._job_wait)
        if operation == "stream":
            try:
                self._read_stream(*args)
            except Exception as exc:
                stream = self._streams.get(args[0])
                if stream is not None:
                    stream["error"] = f"Input stream failed: {exc}"
                    self._stream_ready(args[0])
            return
        key, region = args
        replay = self._replays[key]
        if replay.error:
            return
        first = region * replay.frames
        stop = min(first + replay.frames, replay.array.shape[0])
        working_bytes = replay.array[first:stop].nbytes
        budget = self.storage.budget
        reserved = budget is None or budget.reserve_memory(working_bytes, "input")
        if not reserved:
            if priority > 0:
                replay.states[region] = "not_requested"
                self._metric("input_speculation_declined_memory")
                return
            replay.error = f"Insufficient RAM budget for a {working_bytes:,}-byte input region"
            replay.states[region] = "failed"
            self._fulfill()
            return
        replay.states[region] = "loading"
        try:
            self._read_direct(self._dataset(key), replay.array[first:stop], np.s_[first:stop])
            flushed = perf_counter()
            replay.array.flush()
            self._metric("input_replay_flush_seconds", perf_counter() - flushed)
            # Readiness is published only after a successful complete read and
            # flush. No job ever writes this region again.
            replay.states[region] = "ready"
            self.source_reads.setdefault(key, 1)
            self._metric("input_replay_bytes", replay.array[first:stop].nbytes)
            if sum(state == "ready" for state in replay.states) == 1:
                self._metric("input_first_region_seconds", perf_counter() - replay.created)
            if all(state == "ready" for state in replay.states):
                self.cache[key] = replay.reference
                self._metric("input_replays_completed")
        except Exception as exc:
            replay.error = f"Could not read {key} region {region}: {exc}"
            replay.states[region] = "failed"
            self._metric("input_replays_failed")
        finally:
            if budget is not None:
                budget.release_memory(working_bytes, "input")
        self._fulfill()

    def _owned_stream(self, owner, stream_id):
        stream = self._streams[stream_id]
        if stream["owner"] != owner:
            raise RuntimeError("Input stream belongs to another client")
        return stream

    def _open_stream(self, owner, slot, path, frames, max_bytes):
        if len(self._streams) >= self.policy.max_streams:
            raise MemoryError("Input stream capacity exhausted; close another stream first")
        dataset = self._dataset((slot, path))
        if not dataset.shape or dataset.dtype.hasobject:
            raise TypeError("Streaming requires a nonscalar fixed-size dtype")
        frames = self.policy.stream_frames if frames is None else frames
        if isinstance(frames, bool) or not isinstance(frames, (int, np.integer)) or frames < 1:
            raise ValueError("frames must be a positive integer")
        frames = self._frames(dataset, int(frames))
        if max_bytes is not None:
            if isinstance(max_bytes, bool) or not isinstance(max_bytes, (int, np.integer)) or max_bytes < 1:
                raise ValueError("max_bytes must be a positive integer")
            frames = min(frames, max(1, max_bytes // max(1, int(np.prod(dataset.shape[1:])) * dataset.dtype.itemsize)))
        budget = self.storage.budget
        count = self.policy.stream_buffers
        while True:
            shape = (frames, *dataset.shape[1:])
            size = int(np.prod(shape)) * dataset.dtype.itemsize
            if size == 0:
                raise ValueError("Cannot stream empty spatial dimensions")
            reserved = size * count
            if budget is None or budget.reserve_memory(reserved, "input"):
                break
            # Preserve a bounded, useful serial stream under tighter budgets.
            # Fail explicitly if even a single source frame cannot be admitted.
            if count > 1:
                count -= 1
            elif frames > 1:
                frames = max(1, frames // 2)
            else:
                raise MemoryError(f"Input stream needs at least {reserved:,} reserved bytes; RAM budget exhausted")
            self._metric("input_stream_capacity_reductions")
        buffers = []
        try:
            for _ in range(count):
                memory = SharedMemory(create=True, size=size)
                buffers.append({"memory": memory, "state": "free", "start": None, "stop": None})
        except BaseException:
            for buffer in buffers:
                buffer["memory"].close()
                buffer["memory"].unlink()
            if budget is not None:
                budget.release_memory(reserved, "input")
            raise
        self._stream_serial += 1
        identity = self._stream_serial
        self._streams[identity] = dict(owner=owner, key=self.aliases[slot, path], buffers=buffers, frames=frames,
            shape=shape, dtype=dataset.dtype, reserved=reserved, next_fill=0, next_take=0,
            waiting=None, error=None, total=dataset.shape[0],
            replay=self.manifest.get(self.aliases[slot, path], {}).get("strategy") == "replay")
        self._metric("input_stream_reserved_bytes", reserved)
        self._metric("input_stream_active_bytes", reserved)
        self._metric("input_stream_active_buffers", len(buffers))
        self._metric("input_stream_buffers_created", len(buffers))
        self.storage.metrics.maximum("input_stream_peak_buffers", sum(len(item["buffers"]) for item in self._streams.values()))
        self.storage.metrics.maximum("input_stream_peak_bytes", sum(item["reserved"] for item in self._streams.values()))
        self._fill_free(identity)
        self._reply(owner, identity)

    def _fill_free(self, identity):
        stream = self._streams[identity]
        total = stream["total"]
        for index, buffer in enumerate(stream["buffers"]):
            if buffer["state"] == "free" and stream["next_fill"] < total and not stream["error"]:
                start = stream["next_fill"]
                stop = min(start + stream["frames"], total)
                buffer.update(state="filling", start=start, stop=stop)
                stream["next_fill"] = stop
                self._enqueue(("stream", identity, index), 1, "stream", identity, index)

    def _read_stream(self, identity, index):
        stream = self._streams.get(identity)
        if stream is None:
            return
        buffer = stream["buffers"][index]
        key = self.aliases[stream["key"]]
        replay = self._replay(key) if stream["replay"] or key in self._replays else None
        if replay is not None:
            regions = range(buffer["start"] // replay.frames, (buffer["stop"] - 1) // replay.frames + 1)
            missing = next((r for r in regions if replay.states[r] != "ready"), None)
            if missing is not None and not replay.error:
                self._queue_region(key, missing, 0)
                self._enqueue(("stream", identity, index), 2, "stream", identity, index)
                return
        block = np.ndarray((buffer["stop"] - buffer["start"], *stream["shape"][1:]),
                           dtype=stream["dtype"], buffer=buffer["memory"].buf)
        try:
            if replay is None:
                self._read_direct(self._dataset(stream["key"]), block,
                                  np.s_[buffer["start"]:buffer["stop"]])
            else:
                if missing is not None:
                    raise OSError(replay.error)
                replay.demanded.update(regions)
                self._metric("input_region_reuses", sum(r in replay.served for r in regions))
                replay.served.update(regions)
                np.copyto(block, replay.array[buffer["start"]:buffer["stop"]])
            buffer["state"] = "ready"
        except Exception as exc:
            buffer["state"] = "failed"
            stream["error"] = f"Input stream read failed: {exc}"
        finally:
            del block
        self._stream_ready(identity)

    def _stream_ready(self, identity):
        stream = self._streams[identity]
        if stream["waiting"] is None:
            return
        owner = stream["owner"]
        if stream["error"]:
            self._reply(owner, error=stream["error"])
            stream["waiting"] = None
            return
        total = stream["total"]
        if stream["next_take"] >= total:
            self._reply(owner)
            stream["waiting"] = None
            return
        for index, buffer in enumerate(stream["buffers"]):
            if buffer["start"] != stream["next_take"]:
                continue
            if buffer["state"] != "ready":
                self._enqueue(("stream", identity, index), 0, "stream", identity, index)
                return
            buffer["state"] = "in_use"
            self.storage.metrics.maximum("input_stream_peak_leases", sum(
                buf["state"] == "in_use" for item in self._streams.values() for buf in item["buffers"]))
            shape = (buffer["stop"] - buffer["start"], *stream["shape"][1:])
            self._reply(owner, (index, buffer["memory"].name, shape, stream["dtype"],
                                buffer["start"], buffer["stop"]))
            stream["next_take"] = buffer["stop"]
            self._metric("input_stream_wait_seconds", perf_counter() - stream["waiting"])
            stream["waiting"] = None
            return

    def _close_stream(self, identity):
        stream = self._streams.pop(identity)
        for index, buffer in enumerate(stream["buffers"]):
            self._jobs.pop(("stream", identity, index), None)
            buffer["memory"].close()
            buffer["memory"].unlink()
        if self.storage.budget is not None:
            self.storage.budget.release_memory(stream["reserved"], "input")
        self._metric("input_stream_active_bytes", -stream["reserved"])
        self._metric("input_stream_active_buffers", -len(stream["buffers"]))
        self._metric("input_streams_closed")

    def _serve(self):
        try:
            while not self._stop.is_set():
                # Drain bounded traffic before each I/O unit; demand can promote
                # an existing job instead of duplicating it.
                for _ in range(32):
                    try:
                        request = self.requests.get_nowait()
                    except Empty:
                        break
                    self._handle(request)
                if self._stop.is_set():
                    break
                if self._jobs:
                    key = min(self._jobs, key=lambda key: self._jobs[key][:2])
                    self._run_job(self._jobs.pop(key))
                else:
                    try:
                        self._handle(self.requests.get(timeout=0.02))
                    except Empty:
                        pass
        finally:
            self.stopped.set()
            for identity in list(self._pending):
                self._reply(identity, error="Input service cancelled")
            self._pending.clear()
            self._jobs.clear()
            for identity in list(self._streams):
                stream = self._streams[identity]
                if stream["waiting"] is not None:
                    self._reply(stream["owner"], error="Input stream cancelled")
                self._close_stream(identity)
            for replay in self._replays.values():
                unused = sum(min(replay.frames, replay.array.shape[0] - r * replay.frames)
                             for r, state in enumerate(replay.states)
                             if state == "ready" and r not in replay.demanded)
                self._metric("input_unused_prefetch_bytes", unused * int(np.prod(replay.array.shape[1:]))
                             * replay.array.dtype.itemsize)
                if not all(state == "ready" for state in replay.states):
                    self._metric("input_replays_cancelled")
                replay.array._mmap.close()
            for dataset in self._datasets.values():
                dataset.id.close()
            self._datasets.clear()
            useful = []
            for owner, first, last in self.compute_intervals:
                useful.extend(_subtract_intervals((first, last),
                    [(a, b) for client, a, b in self.wait_intervals if client == owner]))
            self._metric("input_overlap_valid", int(not self._intervals_dropped))
            if not self._intervals_dropped:
                self._metric("input_io_compute_overlap_seconds", interval_overlap(self.io_intervals, useful))
            Logger.log_perf(f"input_summary metrics={self.storage.metrics.snapshot()}")

    def _stage_small(self, dataset):
        started = perf_counter()
        value = np.asarray(dataset[()])
        ended = perf_counter()
        self._record_interval(self.io_intervals, (started, ended))
        self._metric("input_hdf5_small_read_seconds", ended - started)
        self._metric("source_reads")
        self._metric("input_logical_bytes", value.nbytes)
        if value.dtype.kind == "O":
            return self.storage.store_value(value)
        return self.storage.store(value, category="input")

    def close(self):
        if self._closed:
            return
        self._closed = True
        if self.thread.is_alive():
            self.requests.put((None, "shutdown", (), perf_counter()))
        self.thread.join()
        # Clients/workers must close first (runtime ExitStack guarantees this).
        for replay in self._replays.values():
            Path(replay.reference.path).unlink(missing_ok=True)
        for queue in [self.requests, *self.responses.values()]:
            queue.cancel_join_thread()
            queue.close()
