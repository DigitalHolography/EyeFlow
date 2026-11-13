"""Adaptive source access: worker-local HDF5, bounded transfers and replay."""

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
import os
import struct
from math import prod

import numpy as np

from input_output.readers.nodes import is_dataset

from pipeline_engine.data.storage import ArrayReference, ValueReference
from pipeline_engine.data.input_replay import ReplayMappings
from pipeline_engine.data.input_io import SourceReader, SourceJob
from pipeline_engine.stage_graph.model import InputAccess, StageGraph
from utils.logger import Logger
from pipeline_engine.data.input_location import storage_location


@dataclass
class SourceEntry:
    attrs: dict
    shape: tuple | None = None
    dtype: object = None
    identity: str | None = None
    chunks: tuple | None = None
    compression: str | None = None
    backend: str = "nas_direct"

    @property
    def size(self):
        return prod(self.shape)

    @property
    def nbytes(self):
        return self.size * self.dtype.itemsize


def _tiles(shape, limit, itemsize, origin=None):
    """Basic rectangular selections, including sub-frame tiles, bounded in bytes."""
    origin = (0,) * len(shape) if origin is None else origin
    if not prod(shape):
        return
    if prod(shape) * itemsize <= limit:
        yield tuple(slice(first, first + size) for first, size in zip(origin, shape))
        return
    axis = next((i for i, size in enumerate(shape) if size > 1), None)
    if axis is None:
        raise MemoryError("One source element exceeds max_io_bytes")
    step = max(1, min(shape[axis] // 2, limit // max(1, itemsize * prod(shape[axis + 1:]))))
    for start in range(0, shape[axis], step):
        sizes, offsets = list(shape), list(origin)
        sizes[axis] = min(step, shape[axis] - start)
        offsets[axis] += start
        yield from _tiles(tuple(sizes), limit, itemsize, tuple(offsets))


def _direct_extent(shape, selection):
    """Crop all basic axes before transport; NumPy still applies the selector."""
    items = selection if isinstance(selection, tuple) else (selection,)
    if not shape:
        return (), (), selection
    basic = all(item is None or item is Ellipsis or isinstance(item, (slice, int, np.integer))
                and not isinstance(item, (bool, np.bool_)) for item in items)
    consumed = sum(item is not None and item is not Ellipsis for item in items)
    ellipses = sum(item is Ellipsis for item in items)
    if not basic or consumed > len(shape) or ellipses > 1:
        return (0,) * len(shape), shape, selection
    expanded = []
    for item in items:
        expanded.extend([slice(None)] * (len(shape) - consumed) if item is Ellipsis else [item])
    if not ellipses:
        expanded.extend([slice(None)] * (len(shape) - consumed))
    origin, sizes, local = [], [], []
    axis = 0
    for item in expanded:
        if item is None:
            local.append(None)
            continue
        size = shape[axis]
        axis += 1
        if isinstance(item, (int, np.integer)):
            index = int(item)
            if not -size <= index < size:
                raise IndexError("Input index out of range")
            origin.append(index % size)
            sizes.append(1)
            local.append(0)
        else:
            selected = range(*item.indices(size))
            if not selected:
                origin.append(0)
                sizes.append(0)
                local.append(slice(0, 0))
            else:
                first, last = selected[0], selected[-1]
                origin.append(min(first, last))
                extent = abs(last - first) + 1
                sizes.append(extent)
                local.append(slice(0, extent, selected.step) if selected.step > 0 else
                             slice(extent - 1, None, selected.step))
    return tuple(origin), tuple(sizes), tuple(local)


def _native_cast(source, target):
    # HDF5 clips some narrowing float/integer conversions where NumPy wraps or
    # returns a sentinel. Keep those conversions on the existing NumPy path.
    target = np.dtype(target)
    return (source.kind in "biufc" and target.kind in "biufc"
            and np.can_cast(source, target, casting="safe"))


def _native_selection(selection):
    items = selection if isinstance(selection, tuple) else (selection,)
    return all(item is Ellipsis or isinstance(item, slice) and (item.step is None or item.step > 0)
               or isinstance(item, (int, np.integer)) and not isinstance(item, (bool, np.bool_))
               for item in items)


class _LocalReader:
    """One lazy read-only handle per source and process; guarded by client lock."""

    def __init__(self):
        self.pid = os.getpid()
        self.files = {}
        self.datasets = {}

    def dataset(self, filename, path):
        import h5py
        if self.pid != os.getpid():
            raise RuntimeError("HDF5 input readers must be created after process creation")
        if filename not in self.files:
            self.files[filename] = h5py.File(filename, "r")
        key = filename, path
        if key not in self.datasets:
            self.datasets[key] = self.files[filename][path]
        return self.datasets[key]

    def read(self, filename, path, selection, dtype=None):
        dataset = self.dataset(filename, path)
        reader = dataset if dtype is None else dataset.astype(dtype)
        # h5py supports fewer selectors than NumPy (negative strides, newaxis,
        # unordered fancy indices). Preserve the existing NumPy proxy contract.
        started = perf_counter()
        try:
            value = reader[selection]
            self.last_seconds = perf_counter() - started
            self.last_interval = started, started + self.last_seconds
            self.last_bytes = np.asarray(value).size * dataset.dtype.itemsize
            return value
        except (TypeError, ValueError):
            origin, shape, local = _direct_extent(dataset.shape, selection)
            extent = tuple(slice(offset, offset + size) for offset, size in zip(origin, shape))
            value = reader[extent]
            self.last_seconds = perf_counter() - started
            self.last_interval = started, started + self.last_seconds
            self.last_bytes = np.asarray(value).size * dataset.dtype.itemsize
            return value[local]

    def read_into(self, filename, path, destination, source_sel, dest_sel):
        dataset = self.dataset(filename, path)
        if (not isinstance(destination, np.ndarray) or not destination.flags.c_contiguous
                or not destination.flags.writeable or not destination.size
                or not dataset.shape or not dataset.size
                or not _native_cast(dataset.dtype, destination.dtype)
                or source_sel is not None and not _native_selection(source_sel)
                or dest_sel is not None and not _native_selection(dest_sel)):
            return False
        started = perf_counter()
        try:
            dataset.read_direct(destination, source_sel=source_sel, dest_sel=dest_sel)
        except (TypeError, ValueError, ZeroDivisionError):
            # NumPy selectors, noncontiguous destinations and casts retain the
            # proxy's previous semantics through read_slice + assignment.
            return False
        self.last_seconds = perf_counter() - started
        self.last_interval = started, started + self.last_seconds
        # Count selected source bytes, independent of destination dtype/broadcast.
        items = () if source_sel is None else source_sel if isinstance(source_sel, tuple) else (source_sel,)
        missing = len(dataset.shape) - sum(item is not Ellipsis for item in items)
        expanded = []
        for item in items:
            expanded.extend([slice(None)] * missing if item is Ellipsis else [item])
        expanded.extend([slice(None)] * (len(dataset.shape) - len(expanded)))
        self.last_bytes = prod(len(range(*item.indices(size))) if isinstance(item, slice) else 1
                               for size, item in zip(dataset.shape, expanded)) * dataset.dtype.itemsize
        return True

    def close(self):
        for dataset in self.datasets.values():
            dataset.id.close()
        self.datasets.clear()
        for source in self.files.values():
            source.close()
        self.files.clear()


class InputClient:
    def __init__(self, identity, requests, responses, stopped, policy):
        self.identity, self.requests, self.responses = identity, requests, responses
        self.stopped = stopped
        self.cache = {}
        self._lock = RLock()
        self._cache_lock = RLock()
        self._metric_lock = RLock()
        self._streams = set()
        self.policy = policy
        self._local = None
        self._metrics = {}
        self._local_io = deque(maxlen=8192)
        self._local_io_dropped = False
        self._closed = False
        self._budget = None
        self._replay_maps = ReplayMappings(self.metric)

    def __getstate__(self):
        if self._streams:
            raise RuntimeError("Cannot transfer an input client with active leases")
        # Spawn transports controls, never previously opened mappings/arrays.
        return {**self.__dict__, "_lock": None, "_cache_lock": None, "_metric_lock": None, "cache": {}, "_streams": set(),
                "_local": None, "_metrics": {}, "_budget": None, "_replay_maps": None,
                "_local_io": deque(maxlen=8192), "_local_io_dropped": False}

    def __setstate__(self, state):
        self.__dict__.update(state)
        self._lock = RLock()
        self._cache_lock = RLock()
        self._metric_lock = RLock()
        self._replay_maps = ReplayMappings(self.metric)

    def rpc(self, operation, *args):
        # A client has one response queue. Serialize threads as well as requests
        # so concurrent reads cannot consume each other's reply.
        with self._lock:
            started = perf_counter()
            if self._closed or self.stopped.is_set():
                raise OSError("Input service has stopped")
            self.send(operation, *args)
            while True:
                try:
                    value, error = self.responses.get(timeout=0.1)
                    break
                except Empty:
                    if self.stopped.is_set():
                        raise OSError("Input service stopped while waiting for a response") from None
            if operation in ("read", "stream_next", "direct_read"):
                ended = perf_counter()
                self.metric("input_rpc_roundtrip_seconds", ended - started)
                self.metric("input_read_rpcs", 1)
                self.send("input_wait", started, ended)
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

    def metric(self, name, value):
        with self._metric_lock:
            self._metrics[name] = self._metrics.get(name, 0) + value

    def _local_reader(self):
        if self._closed or self.stopped.is_set():
            raise OSError("Input client has closed")
        if self._local is None:
            self._local = _LocalReader()
        return self._local

    def _local_read_metrics(self):
        self.metric("input_local_hdf5_read_seconds", self._local.last_seconds)
        self.metric("input_local_logical_bytes", self._local.last_bytes)
        self.metric("input_local_source_reads", 1)
        self.metric("input_consumer_blocked_seconds", self._local.last_seconds)
        if len(self._local_io) == self._local_io.maxlen:
            self._local_io_dropped = True
        self._local_io.append(self._local.last_interval)

    def local_read_into(self, filename, path, destination, source_sel, dest_sel):
        with self._lock:
            if not self._local_reader().read_into(filename, path, destination, source_sel, dest_sel):
                return False
            self._local_read_metrics()
            self.metric("input_local_native_direct_reads", 1)
            return True

    def local_read(self, filename, path, selection, dtype=None):
        with self._lock:
            value = self._local_reader().read(filename, path, selection, dtype=dtype)
            self._local_read_metrics()
            if isinstance(value, np.ndarray):
                value.flags.writeable = False
            return value

    def direct_read(self, slot, path, shape, dtype, selection):
        # Keep the same NumPy indexing contract as replay. Complex indexing may
        # require a full worker array, but never a full coordinator allocation.
        with self._lock:
            if self._closed or self.stopped.is_set():
                raise OSError("Input client has closed")
            if dtype.hasobject:
                return self.read(slot, path)[selection]
            origin, target_shape, local = _direct_extent(shape, selection)
            result = np.empty(target_shape, dtype=dtype)
            for tile in _tiles(target_shape, self.policy.max_io_bytes, dtype.itemsize):
                source = tuple(slice(item.start + offset, item.stop + offset) for item, offset in zip(tile, origin))
                name, block_shape, token = self.rpc("direct_read", slot, path, source)
                memory = None
                try:
                    memory = SharedMemory(name=name)
                    started = perf_counter()
                    view = np.ndarray(block_shape, dtype=dtype, buffer=memory.buf)
                    result[tile] = view
                    del view
                    self.metric("input_direct_copy_seconds", perf_counter() - started)
                finally:
                    if memory is not None:
                        memory.close()
                    self.send("direct_release", token)
            result.flags.writeable = False
            value = result[local]
            if isinstance(value, np.ndarray):
                value.flags.writeable = False
            return value

    def read(self, slot, path):
        return self._read(slot, path)

    def _check_open(self):
        if self._closed or self.stopped.is_set():
            raise OSError("Input client has closed or stopped")

    def _read(self, slot, path):
        key = slot, path
        with self._cache_lock:
            self._check_open()
            if key in self.cache:
                return self.cache[key]
            cached = self._replay_maps.find(key)
            if cached is not None:
                return cached[0]
        reference = self.rpc("read", slot, path, None)
        with self._cache_lock:
            self._check_open()
            if isinstance(reference, ArrayReference):
                stop = reference.shape[0] if reference.shape else 1
                return self._replay_maps.authorize(key, 0, stop, True, reference)
            array = reference.open() if isinstance(reference, ValueReference) else reference
            if isinstance(array, np.ndarray):
                array.flags.writeable = False
            self.cache[key] = array
            return array

    def read_slice(self, slot, path, shape, selection, dtype=None):
        span, local = _selection_span(shape, selection)
        key = slot, path
        if span is None:
            return self._read(slot, path)[selection]
        first, stop = span
        with self._cache_lock:
            self._check_open()
            if key in self.cache:
                return self.cache[key][selection]
            cached = self._replay_maps.find(key, first, stop) if first != stop else None
        if first == stop:
            dtype = dtype if dtype is not None else self.rpc("dtype", slot, path)
            value = np.empty((0, *shape[1:]), dtype=dtype)
            value.flags.writeable = False
            return value[local]
        if cached is not None:
            array, origin = cached
        else:
            # The RPC lock serializes the protocol, never unrelated cache hits.
            reference = self.rpc("read", slot, path, span)
            with self._cache_lock:
                self._check_open()
                if not isinstance(reference, (_RegionReference, ArrayReference)):
                    value = reference.open() if isinstance(reference, ValueReference) else reference
                    if isinstance(value, np.ndarray):
                        value.flags.writeable = False
                    self.cache[key] = value
                    return value[selection]
                origin = first if isinstance(reference, _RegionReference) else 0
                end = stop if isinstance(reference, _RegionReference) else shape[0]
                array = self._replay_maps.authorize(key, origin, end, origin == 0 and end == shape[0], reference)
        started = perf_counter()
        temporal = first - origin if isinstance(local[0], int) else slice(first - origin, stop - origin)
        value = array[(temporal, *local[1:])]
        self.metric("input_replay_selection_seconds", perf_counter() - started)
        return value

    def close(self):
        with self._lock, self._cache_lock:
            self._close()

    def _close(self):
        if self._closed:
            return
        for stream in list(self._streams):
            stream.close()
        if self._local is not None:
            self._local.close()
        if self._local_io:
            # One compact, bounded binary batch, never a diagnostic object graph
            # or an IPC round trip per worker-local numerical read.
            packed = b"".join(struct.pack("dd", *interval) for interval in self._local_io)
            self.send("local_io", packed, self._local_io_dropped)
            self._local_io.clear()
        self.metric("input_replay_client_peak_retained_mappings", self._replay_maps.peak_retained)
        self.metric("input_replay_client_peak_live_mappings", self._replay_maps.peak_live)
        with self._metric_lock:
            metrics = dict(self._metrics)
        self.send("metrics", metrics)
        self.send("client_close")
        for value in self.cache.values():
            if isinstance(value, np.memmap) and not value._mmap.closed:
                value._mmap.close()
        self._replay_maps.close()
        self.cache.clear()
        self._closed = True


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
        if self.file.entries[self.path].backend != "replay":
            return self.read_slice(())
        identity = self.file.entries[self.path].identity
        return self.file.client.read(self.file.slot, identity)

    def __getitem__(self, selection):
        return self.read_slice(selection)

    def read_slice(self, selection):
        backend = self.file.entries[self.path].backend
        if backend == "local_direct":
            return self.file.client.local_read(self.file.filename, self.path, selection)
        if backend in ("nas_direct", "stream"):
            return self.file.client.direct_read(self.file.slot, self.path, self.shape, self.dtype, selection)
        value = self.file.client.read_slice(self.file.slot, self.file.entries[self.path].identity, self.shape, selection, self.dtype)
        if isinstance(value, np.ndarray):
            value.flags.writeable = False
        return value

    def read_direct(self, dest, source_sel=None, dest_sel=None):
        if self.file.entries[self.path].backend == "local_direct" and self.file.client.local_read_into(
            self.file.filename, self.path, dest, source_sel, dest_sel,
        ):
            return
        value = self.read_slice(() if source_sel is None else source_sel)
        if dest_sel is None:
            np.copyto(dest, value, casting="unsafe")
        else:
            dest[dest_sel] = value

    def iter_blocks(self, *, frames=None, writable=False, max_bytes=None):
        """Use as a context manager; each yielded lease owns one exclusive block."""
        stream = _LocalBlockStream if self.file.entries[self.path].backend == "local_direct" else _BlockStream
        return stream(self, frames, writable, max_bytes)

    def __array__(self, dtype=None, copy=None):
        array = self.astype(dtype) if dtype is not None else np.asarray(self._value)
        array.flags.writeable = False
        return array.copy() if copy else array

    def astype(self, dtype):
        if (dtype is not None and self.file.entries[self.path].backend == "local_direct"
                and _native_cast(self.dtype, dtype)):
            value = np.asarray(self.file.client.local_read(self.file.filename, self.path, (), dtype=np.dtype(dtype)))
            if np.dtype(dtype) != self.dtype:
                value.flags.writeable = True
            return value
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
    direct_buffers: int = 2
    max_prefetch_jobs: int = 16
    lookahead: bool = True
    future_prefetch: bool = field(default_factory=lambda: os.environ.get("EYEFLOW_FUTURE_PREFETCH", "1") != "0")
    prefetch_regions: int = 16
    max_prefetch_bytes: int = 1000 * 2**20
    bulk_bytes: int = 1024**2
    input_backend: str = "auto"
    unknown_transport: str = "coordinator"

    def __post_init__(self):
        if self.input_backend not in ("auto", "direct", "replay"):
            raise ValueError("input_backend must be auto, direct or replay")
        if self.unknown_transport not in ("coordinator", "worker"):
            raise ValueError("unknown_transport must be coordinator or worker")
        for name in ("replay_frames", "stream_frames", "stream_buffers", "max_io_bytes", "max_streams", "direct_buffers",
                     "max_prefetch_jobs", "prefetch_regions", "max_prefetch_bytes", "bulk_bytes"):
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
    tiles: dict = field(default_factory=dict)
    prefetched: set[int] = field(default_factory=set)


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


class _LocalLease:
    def __init__(self, block, start, stop, budget=None):
        self.block, self.start, self.stop = block, start, stop
        self.released = False
        self.budget = budget
        self.reserved = block.nbytes

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
        if self.budget is not None:
            self.budget.release_memory(self.reserved, "input")


class _LocalBlockStream:
    """Worker-owned blocks obey the same exclusive lease API without IPC."""
    def __init__(self, dataset, frames, writable, max_bytes):
        if not dataset.shape or dataset.dtype.hasobject:
            raise TypeError("Streaming requires a nonscalar fixed-size dtype")
        self.dataset, self.client, self.writable = dataset, dataset.file.client, writable
        frames = self.client.policy.stream_frames if frames is None else frames
        for name, value in (("frames", frames), ("max_bytes", max_bytes)):
            if value is not None and (isinstance(value, bool) or not isinstance(value, (int, np.integer)) or value < 1):
                raise ValueError(f"{name} must be a positive integer")
        limit = min(self.client.policy.max_io_bytes, max_bytes or self.client.policy.max_io_bytes)
        frame_bytes = prod(dataset.shape[1:]) * dataset.dtype.itemsize
        if not frame_bytes or frame_bytes > limit:
            raise MemoryError("One stream frame exceeds the input byte limit")
        self.frames = min(frames, limit // frame_bytes)
        self.next_frame, self.lease, self.closed = 0, None, False
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
        if self.next_frame >= self.dataset.shape[0]:
            self.close()
            raise StopIteration
        first, stop = self.next_frame, min(self.next_frame + self.frames, self.dataset.shape[0])
        budget = self.client._budget
        frame_bytes = prod(self.dataset.shape[1:]) * self.dataset.dtype.itemsize
        while budget is not None and not budget.reserve_memory((stop - first) * frame_bytes, "input"):
            if stop - first == 1:
                self.close()
                raise MemoryError("Insufficient RAM budget for one local input frame")
            stop = first + max(1, (stop - first) // 2)
        reserved = (stop - first) * frame_bytes
        try:
            block = self.dataset[first:stop]
            if self.writable:
                block.flags.writeable = True
            self.lease = _LocalLease(block, first, stop, budget)
            self.next_frame = stop
            return self.lease
        except Exception:
            if budget is not None:
                budget.release_memory(reserved, "input")
            self.close()
            raise

    def close(self):
        if not self.closed:
            self.closed = True
            if self.lease is not None:
                self.lease.release()
            self.client._streams.discard(self)


class InputService:
    """Authoritative metadata, bounded source reads and input-buffer ownership.

    The control thread alone mutates jobs, readiness and leases. SourceReader
    owns one independent blocking I/O lane and exclusive in-flight destinations.
    Scientific products never enter this cache. Sources must outlive close().
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
        self._plan_order = {}
        self._stage_order = {}
        self._viable = set()
        self._active = set()
        self._prefetch_cancelled = set()
        self._prefetch_storage_declined = set()
        self._prefetch_paused = False
        self._plan_retry = 0.0
        self._idle_reason = "exhausted"
        self._idle_cause = "exhausted"
        self._idle_detail = "not_configured"
        self._idle_period = None
        self.dataset_metrics = {}
        self._undeclared = set()
        self._replays, self._jobs, self._pending, self._streams = {}, {}, {}, {}
        self._datasets = {}
        self._inflight = None
        self._io_serial = 0
        self._small_pending = {}
        self._open_pending = {}
        self._transfers = {}
        self._direct_pending = {}
        self._transfer_serial = 0
        self.locations = {}
        self._serial = 0
        self._stream_serial = 0
        self._stream_planned = set()
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
            self.locations[slot] = storage_location(source.filename)
            Logger.log_perf(f"input_source slot={slot} storage={self.locations[slot]}")
        self._select_backends()
        self._metric("input_metadata_seconds", perf_counter() - started)
        self.reader = SourceReader(sources, lambda *args: self._read_direct(*args),
                                   lambda key: self._dataset(key), self._wake_completion)
        self._datasets = self.reader.datasets  # Reader-owned; retained for diagnostics.
        self.thread = Thread(target=copy_context().run, args=(self._serve,), name="EyeFlowInputs", daemon=True)
        with Logger.perf_stage("run/input_thread_start"):
            self.reader.thread.start()
            self.thread.start()

    def _metric(self, name, value=1):
        self.storage.metrics.add(name, value)

    def _record_interval(self, intervals, interval):
        if len(intervals) == intervals.maxlen:
            self._intervals_dropped = True
            self._metric("input_interval_samples_dropped")
        intervals.append(interval)

    def configure(self, graphs: Mapping[str, StageGraph | None],
                  run_settings: Mapping[str, object] | None = None, *,
                  expected_invocations: Mapping[str, tuple[str, ...]] | None = None,
                  execution_order: tuple[str, ...] | None = None) -> None:
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
            # Use the local DAG's topological order, including external products,
            # without evaluating any stage callbacks or conditions.
            produced = {product.id for stage in graph.stages for product in stage.products}
            external = {name for stage in graph.stages
                        for name in (*stage.requires, *stage.optional_requires)} - produced
            plan = graph.plan(external_products=external)
            self._stage_order[pipeline] = {name: i for i, name in enumerate(plan.order)}
            for stage_id in plan.order:
                stage = plan.stage_by_id[stage_id]
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
                            "entry": entry, "consumers": [], "access": set()})
                        item["consumers"].append((pipeline, stage.id, use))
                        item["access"].add(use.access)
        for item in self.manifest.values():
            item["potential_consumers"] = len({
                (pipeline, stage, invocation)
                for pipeline, stage, use in item["consumers"] if use.access != InputAccess.METADATA
                for invocation in (expected_invocations or {}).get(pipeline, (None,))
                if invocation is None or use.when is None or use.when == invocation
            })
            item["expected_consumers"] = len({
                (pipeline, stage, invocation)
                for pipeline, stage, use in item["consumers"]
                if not use.optional and use.access != InputAccess.METADATA
                for invocation in (expected_invocations or {}).get(pipeline, (None,))
                if invocation is None or use.when is None or use.when == invocation
            })
        self._configured = True
        self._select_backends()
        for left, right in (("hd.moment0", "hd.moment2"), ("hd.band_low", "hd.band_high")):
            if all(source in semantic and any(
                not use.optional and use.access == InputAccess.SEQUENTIAL
                for _, _, use in self.manifest[semantic[source]]["consumers"])
                   for source in (left, right)):
                self.pairs[semantic[left]] = semantic[right]
                self.pairs[semantic[right]] = semantic[left]
        self._metric("input_manifest_seconds", perf_counter() - started)
        Logger.log_perf(f"input_manifest datasets={len(self.manifest)} unresolved={len(self.unresolved)} "
                        f"backends={[(key, item['backend']) for key, item in self.manifest.items()]}")

        if execution_order is not None:
            self.requests.put((None, "plan", (tuple(execution_order),), perf_counter()))

    def update_prefetch(self, viable) -> None:
        """Scheduler-owned eligibility; computational readiness is not required."""
        self.requests.put((None, "viable", (tuple(viable),), perf_counter()))

    def _prefetch_idle(self, cause, detail):
        # Keep the existing inclusive resource metric; expose its ahead-limit and
        # genuine-constraint components separately in the trace and counters.
        self._idle_reason = "exhausted" if cause == "exhausted" else "resource"
        self._idle_cause, self._idle_detail = cause, detail

    def _set_idle(self, cause=None):
        """Log transitions, not each iteration of the existing request wait."""
        previous = self._idle_period
        if previous is not None and previous[0] == cause:
            return
        now = perf_counter()
        if previous is not None:
            Logger.log_perf(f"input_idle_end reason={previous[0]} at={now:.6f} "
                            f"wall_seconds={now - previous[1]:.6f}")
        self._idle_period = None if cause is None else (cause, now)
        if cause is not None:
            Logger.log_perf(f"input_idle_begin reason={cause} at={now:.6f} detail={self._idle_detail}")

    def _plan_prefetch(self):
        """Select one operation at a boundary, ordered by imminent first use.

        Pipeline/stage proximity precedes depth fairness, so an empty distant
        stream cannot take the allowance from the next regions of nearer inputs.
        Within that first-use group, balance ready regions and keep pairs adjacent.
        """
        blockers = []
        self._prefetch_idle("exhausted", "no_missing_required_regions")
        for stream in self._streams.values():
            if stream["next_fill"] < stream["total"] and all(
                    buffer["state"] != "free" for buffer in stream["buffers"]):
                blockers.append(("ahead_limit", stream["key"], "stream_buffers_full"))
        if not self.policy.lookahead or self._prefetch_paused or not self._viable:
            self._prefetch_idle("exhausted", "planning_disabled_cancelled_or_no_viable_consumers")
            return
        if perf_counter() < self._plan_retry:
            self._prefetch_idle("constraint", "memory_admission_backoff")
            return
        replay_bytes = sum(
            min(replay.frames, replay.array.shape[0] - r * replay.frames)
            * prod(replay.array.shape[1:]) * replay.array.dtype.itemsize
            for replay in self._replays.values()
            for r, state in enumerate(replay.states)
            if state == "ready" and r not in replay.served)
        stream_bytes = sum(stream["reserved"] for stream in self._streams.values() if stream["owner"] is None)
        unused_bytes = replay_bytes + stream_bytes
        source_order = {key: i for i, key in enumerate(self.manifest)}
        candidates = []
        for key, item in self.manifest.items():
            if key in self._prefetch_cancelled or item["backend"] == "local_direct":
                continue
            consumers = [(p, stage) for p, stage, use in item["consumers"]
                         if p in self._viable and not use.optional
                         and use.access != InputAccess.METADATA]
            if not consumers:
                continue
            active = [(p, stage) for p, stage in consumers if p in self._active]
            if not active and not self.policy.future_prefetch:
                continue
            if key in self._prefetch_storage_declined:
                blockers.append(("constraint", key, "scratch_admission_declined"))
                continue
            entry = item["entry"]
            if not entry.shape or entry.dtype.hasobject or prod(entry.shape) * entry.dtype.itemsize < 4096:
                continue
            first_use = min((self._plan_order[p], self._stage_order[p][stage]) for p, stage in (active or consumers))
            priority = 2 if active else 3
            pair_order = min(source_order[key], source_order.get(self.pairs.get(key), source_order[key]))
            if item["backend"] != "replay":
                if key in self._stream_planned or item["access"] - {InputAccess.METADATA} != {InputAccess.STREAM} or any(
                        stream["key"] == key for stream in self._streams.values()):
                    continue
                frames = self._frames(entry, self.policy.stream_frames)
                capacity = frames * prod(entry.shape[1:]) * entry.dtype.itemsize * min(
                    self.policy.stream_buffers, self.policy.max_prefetch_jobs)
                if unused_bytes + capacity > self.policy.max_prefetch_bytes:
                    blockers.append(("ahead_limit", key, "next_stream_pool_bytes", capacity))
                    continue
                if len(self._streams) >= self.policy.max_streams:
                    blockers.append(("ahead_limit", key, "stream_pool_count", self.policy.max_streams))
                    continue
                candidates.append(((priority, *first_use, 0, pair_order, source_order[key]), key, None))
                continue
            replay = self._replays.get(key)
            if replay is not None and replay.error:
                continue
            frames = replay.frames if replay else self._replay_frames(key, entry)
            states = replay.states if replay else ()
            region = next((r for r, state in enumerate(states) if state == "not_requested"), None) if replay else 0
            if region is None:
                continue
            unused = sum(state == "ready" and r not in replay.served for r, state in enumerate(states))
            size = min(frames, entry.shape[0] - region * frames) * prod(entry.shape[1:]) * entry.dtype.itemsize
            if unused >= self.policy.prefetch_regions:
                blockers.append(("ahead_limit", key, "next_region", region, "region_limit", self.policy.prefetch_regions))
                continue
            if unused_bytes + size > self.policy.max_prefetch_bytes:
                blockers.append(("ahead_limit", key, "next_region", region, "needed_bytes", size))
                continue
            candidates.append(((priority, *first_use, unused, pair_order, source_order[key]), key, region))
        queued_priority = min((job[0] for job in self._jobs.values()), default=4)
        for score, key, region in sorted(candidates):
            # Reconsider active inputs before an already queued future fill, but
            # do not accumulate peers behind it or reserve the same ahead budget
            # twice. Demand/lookahead remain entirely with the existing queue.
            if score[0] >= queued_priority:
                break
            try:
                if region is None:
                    self._open_stream(None, *key, None, None, priority=score[0])
                else:
                    self._queue_region(key, region, score[0])
            except (OSError, MemoryError) as exc:
                if region is None:
                    self._plan_retry = perf_counter() + 0.02
                    blockers.append(("constraint", key, "stream_admission", str(exc)))
                    self._metric("input_speculation_declined_memory")
                    continue
                self._prefetch_storage_declined.add(key)
                blockers.append(("constraint", key, "scratch_admission", str(exc)))
                self._metric("input_speculation_declined_storage")
                Logger.log_perf(f"input_prefetch_declined dataset={key} reason={exc}")
                continue
            if self._jobs:
                Logger.log_perf(f"input_plan_selection dataset={key} region={region} "
                                f"priority={score[0]} first_use={score[1:3]} unused_regions={score[3]} "
                                f"pair={self.pairs.get(key)} reason=imminent_first_use_then_balanced_depth "
                                f"replay_unused_bytes={replay_bytes} unopened_stream_bytes={stream_bytes}")
                return
        if blockers:
            cause = "constraint" if any(item[0] == "constraint" for item in blockers) else "ahead_limit"
            self._prefetch_idle(cause, f"replay_unused_bytes={replay_bytes} unopened_stream_bytes={stream_bytes} "
                                f"byte_limit={self.policy.max_prefetch_bytes} blockers={blockers}")

    def _select_backends(self):
        """Conservative access-cost policy; declaration counts are not pass counts."""
        for slot, (_, entries) in self.metadata.items():
            location = self.locations[slot]
            for path, entry in entries.items():
                if entry.shape is None:
                    continue
                item = self.manifest.get(self.aliases[slot, path], {})
                access = item.get("access", set()) - {InputAccess.METADATA}
                potential = item.get("potential_consumers", 0)
                size = prod(entry.shape) * entry.dtype.itemsize
                reusable = InputAccess.RANDOM in access or potential > 1
                backend = "local_direct" if (location == "local" or
                    location == "unknown" and self.policy.unknown_transport == "worker") else "nas_direct"
                reason = "local_storage" if location == "local" else "unknown_fallback" if location == "unknown" else "small_or_undeclared"
                if self.policy.input_backend == "replay":
                    backend, reason = "replay", "forced_replay"
                elif self.policy.input_backend == "direct":
                    reason = "forced_direct"
                elif location == "network" and access and size >= 4096 and not entry.dtype.hasobject:
                    if reusable:
                        backend, reason = "replay", "potential_repeated_regions"
                    elif access == {InputAccess.STREAM} and size >= self.policy.bulk_bytes:
                        backend, reason = "stream", "single_pass_bounded_prefetch"
                    elif InputAccess.SEQUENTIAL in access and size >= self.policy.bulk_bytes:
                        backend, reason = "replay", "sequential_prefetch"
                entry.backend = backend
                if item:
                    item.update(storage_location=location, expected_reuse=reusable, backend=backend, reason=reason)
                if path == entry.identity and self._configured:
                    Logger.log_perf(f"input_dataset name=/{path} source={location} bytes={size} "
                        f"access={sorted(use.value for use in access)} required_consumers={item.get('expected_consumers', 0)} "
                        f"potential_consumers={potential} reuse={'possible' if reusable else 'single_declared_use'} "
                        f"chunks={entry.chunks} compression={entry.compression} backend={backend} reason={reason}")

    def _backend(self, key):
        slot, path = self.aliases[key]
        return self.metadata[slot][1][path].backend

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
        client = InputClient(identity, self.requests, response, self.stopped, self.policy)
        client._budget = self.storage.budget
        return client

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
        if self._inflight is not None and self._inflight["job"][2:4] == (operation, args):
            job = self._inflight["job"]
            if priority < job[0]:
                self._inflight["job"] = (priority, *job[1:])
                self._metric("input_priority_changes")
            return True
        previous = self._jobs.get(key)
        if previous is not None and previous[0] <= priority:
            return True
        if previous is not None:
            self._metric("input_priority_changes")
            Logger.log_perf(f"input_promotion job={key} old={previous[0]} priority={priority} at={perf_counter():.6f}")
        elif priority > 0 and (operation == "replay" or operation == "stream"
                               and self._streams[args[0]]["owner"] is None):
            speculative = [(k, value) for k, value in self._jobs.items()
                           if value[0] > 0 and (value[2] == "replay" or value[2] == "stream"
                                               and self._streams[value[3][0]]["owner"] is None)]
            if len(speculative) >= self.policy.max_prefetch_jobs:
                worst, value = max(speculative, key=lambda item: item[1][:2])
                if value[0] <= priority:
                    return False
                self._discard_job(worst)
        self._serial += 1
        enqueued = perf_counter()
        self._jobs[key] = (priority, self._serial, operation, args, enqueued)
        Logger.log_perf(f"input_enqueue job={key} priority={priority} operation={operation} at={enqueued:.6f}")
        return True

    def _discard_job(self, job_key):
        self._jobs.pop(job_key, None)
        Logger.log_perf(f"input_cancellation job={job_key} at={perf_counter():.6f}")
        if job_key[0] == "replay":
            _, key, region = job_key
            self._replays[key].states[region] = "not_requested"
            self._replays[key].tiles.pop(region, None)
            self._metric("input_speculative_jobs_cancelled")
        elif job_key[0] == "stream":
            self._close_stream(job_key[1])

    def _frames(self, dataset, frames):
        frame_bytes = max(1, int(np.prod(dataset.shape[1:])) * dataset.dtype.itemsize)
        limit = max(1, self.policy.max_io_bytes // frame_bytes)
        frames = min(frames, limit)
        chunks = getattr(dataset, "chunks", None)
        if chunks and chunks[0] <= frames:
            frames = max(chunks[0], frames // chunks[0] * chunks[0])
        return frames

    def _replay_frames(self, key, dataset):
        access = self.manifest.get(key, {}).get("access", set())
        target = dataset.shape[0] if access and access <= {InputAccess.METADATA, InputAccess.SMALL} else self.policy.replay_frames
        return self._frames(dataset, target)

    def _replay(self, key):
        if key in self._replays:
            return self._replays[key]
        from scratch_location import is_local_path
        if not is_local_path(self.storage.directory):
            raise OSError(f"Input replay scratch must be local: {self.storage.directory}")
        dataset = self._entry(key)
        if dataset.dtype.itemsize > self.policy.max_io_bytes:
            raise MemoryError("One source element exceeds max_io_bytes")
        if dataset.dtype.hasobject or not dataset.shape or not dataset.size or dataset.nbytes < 4096:
            return None
        if shutil.disk_usage(self.storage.directory).free < dataset.nbytes + 4096:
            raise OSError("Insufficient local scratch disk space for input replay")
        path = self.storage.path(".npy")
        try:
            array = np.lib.format.open_memmap(path, mode="w+", dtype=dataset.dtype, shape=dataset.shape)
        except BaseException:
            path.unlink(missing_ok=True)
            raise
        frames = self._replay_frames(key, dataset)
        replay = _Replay(array, ArrayReference(str(path), tuple(dataset.shape), dataset.dtype),
                         frames, ["not_requested"] * ((dataset.shape[0] + frames - 1) // frames), array.offset)
        self._replays[key] = replay
        return replay

    def _entry(self, key):
        slot, path = self.aliases[key]
        return self.metadata[slot][1][path]

    def _dataset(self, key):
        return self.reader.dataset(self.aliases[key])

    @property
    def _priority(self):
        return self.reader.current.priority if self.reader.current is not None else 0

    def _queue_region(self, key, region, priority):
        replay = self._replay(key)
        if replay is None or replay.error or not 0 <= region < len(replay.states):
            return
        if replay.states[region] == "ready":
            return
        if self._inflight is not None and self._inflight["job"][2:4] == ("replay", (key, region)):
            job = self._inflight["job"]
            if priority < job[0]:
                self._inflight["job"] = (priority, *job[1:])
                self._metric("input_priority_changes")
                Logger.log_perf(f"input_promotion running_job={job[1]} dataset={key} region={region} priority={priority}")
            return
        if self._enqueue(("replay", key, region), priority, "replay", key, region):
            replay.states[region] = "queued"

    def _has_replay_demand(self, key, region):
        replay = self._replays[key]
        first, stop = region*replay.frames, (region+1)*replay.frames
        return any(k == key and a < stop and first < b for k, a, b, *_ in self._pending.values()) or any(
            stream["key"] == key and not stream.get("closing") and any(
                buffer["state"] == "filling" and buffer["start"] < stop and first < buffer["stop"]
                for buffer in stream["buffers"])
            for stream in self._streams.values())

    def _lookahead(self, key, stop):
        if not self.policy.lookahead:
            return
        # Pair the current region and exactly one following window. Inactive
        # methods have no pair entries and are never pulled into this lookahead.
        windows = []
        for candidate in dict.fromkeys((key, self.pairs.get(key, key))):
            if self._backend(candidate) != "replay":
                continue
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
        started = perf_counter()
        try:
            self._handle_request(request)
        finally:
            elapsed = perf_counter() - started
            self._metric("input_control_handling_seconds", elapsed)
            self.storage.metrics.maximum("input_control_max_handling_seconds", elapsed)

    def _handle_request(self, request):
        identity, operation, args, requested_at = request
        if identity is not None and identity in self._cancelled_clients:
            return
        self._metric("input_queue_wait_seconds", max(0, perf_counter() - requested_at))
        if self._stop.is_set() and operation in ("read", "direct_read", "stream_open", "stream_next"):
            self._reply(identity, error="Input service cancelled")
            return
        try:
            if operation == "read":
                slot, path, span = args
                key = self.aliases[slot, path]
                if self._configured and key not in self.manifest and key not in self._undeclared:
                    Logger.log_debug(f"Undeclared source demand {slot}:{path}; backend={self._backend(key)}")
                    self._undeclared.add(key)
                replay = self._replay(key)
                if replay is None:
                    if key in self.cache:
                        self._reply(identity, self.cache[key])
                    else:
                        self._small_pending.setdefault(key, set()).add(identity)
                        if self._inflight is None or self._inflight["job"][2:4] != ("small", (key,)):
                            self._enqueue(("small", key), 0, "small", key)
                else:
                    first, stop = span if span is not None else (0, replay.array.shape[0])
                    regions = range(first // replay.frames, (stop - 1) // replay.frames + 1)
                    ready = all(replay.states[r] == "ready" for r in regions)
                    self._metric("input_prefetch_hits" if ready else "input_prefetch_misses")
                    if ready:
                        self._metric("input_ready_reply_seconds", perf_counter() - requested_at)
                    replay.demanded.update(regions)
                    self._pending[identity] = (key, first, stop, span is None, requested_at)
                    # Only one missing slab per pending request is queued. Full
                    # materialization cannot create a volume-sized job queue.
                    self._fulfill()
                    self._lookahead(key, stop if span is not None else min(stop, replay.frames))
            elif operation == "direct_read":
                self._enqueue(("direct", identity), 0, "direct", identity, *args)
            elif operation == "direct_release":
                transfer = self._transfers.get(identity)
                if transfer is not None and transfer.get("token") == args[0]:
                    transfer["busy"] = False
                    self._wake_direct()
            elif operation == "metrics":
                for name, value in args[0].items():
                    self._metric(name, value)
            elif operation == "local_io":
                packed, dropped = args
                for first, last in struct.iter_unpack("dd", packed):
                    self._record_interval(self.io_intervals, (first, last))
                    self._record_interval(self.wait_intervals, (identity, first, last))
                if dropped:
                    self._intervals_dropped = True
                    self._metric("input_interval_samples_dropped")
            elif operation == "dtype":
                slot, path = args
                self._reply(identity, self.metadata[slot][1][path].dtype)
            elif operation in ("input_wait", "compute"):
                intervals = self.wait_intervals if operation == "input_wait" else self.compute_intervals
                self._record_interval(intervals, (identity, *args))
                if operation == "input_wait":
                    self._metric("input_consumer_blocked_seconds", args[1] - args[0])
            elif operation == "io_complete":
                self._drain_completions()
            elif operation == "shutdown":
                self._stop.set()
                self.stopped.set()
            elif operation == "plan":
                self._plan_order = {name: i for i, name in enumerate(args[0])}
                self._viable = set(args[0])
            elif operation == "viable":
                self._viable = set(args[0])
                self._release_planned_streams(lambda key: not any(
                    p in self._viable and not use.optional
                    for p, _, use in self.manifest.get(key, {}).get("consumers", ())))
                self._cancel_running_prefetch(lambda key: not any(
                    p in self._viable and not use.optional
                    for p, _, use in self.manifest.get(key, {}).get("consumers", ())))
                self._active.intersection_update(self._viable)
                for job_key, job in list(self._jobs.items()):
                    if job_key[0] == "replay" and job[0] > 0:
                        item = self.manifest.get(job_key[1], {})
                        if not any(p in self._viable and not use.optional
                                   for p, _, use in item.get("consumers", ())):
                            self._discard_job(job_key)
            elif operation == "cancel_prefetch":
                key = self.aliases[args[0]] if args[0] is not None else None
                self._release_planned_streams(lambda candidate: key is None or candidate == key)
                self._cancel_running_prefetch(lambda candidate: key is None or candidate == key)
                if key is None:
                    self._prefetch_paused = True
                else:
                    self._prefetch_cancelled.add(key)
                for job_key, job in list(self._jobs.items()):
                    if job_key[0] == "replay" and job[0] > 0 and (key is None or job_key[1] == key):
                        self._discard_job(job_key)
            elif operation in ("activate", "prefetch"):
                if not self.policy.lookahead:
                    return
                if operation == "activate":
                    self._active.add(args[0])
                    if not self._plan_order:
                        self._plan_order = {name: i for i, name in enumerate(self._stage_order)}
                        self._viable.add(args[0])
                    for job_key, job in list(self._jobs.items()):
                        if job_key[0] == "replay" and job[0] == 3 and any(
                            p == args[0] and not use.optional
                            for p, _, use in self.manifest.get(job_key[1], {}).get("consumers", ())
                        ):
                            self._queue_region(job_key[1], job_key[2], 2)
                else:
                    key, first = args
                    key = self.aliases[key]
                    if self._backend(key) == "replay":
                        replay = self._replay(key)
                        if replay is not None:
                            self._queue_region(key, first // replay.frames, 3)
            elif operation == "stream_open":
                try:
                    self._open_stream(identity, *args)
                except MemoryError:
                    if any(s.get("closing") for s in self._streams.values()):
                        self._open_pending[identity] = args
                    else:
                        raise
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
                if operation == "reclaim" and owner not in self._cancelled_clients:
                    # An ungraceful exit can lose the worker's batched local I/O
                    # intervals. Never report complete overlap for that run.
                    self._intervals_dropped = True
                    self._metric("input_metrics_clients_lost")
                self._cancelled_clients.add(owner)
                self._open_pending.pop(owner, None)
                for pending in self._small_pending.values():
                    pending.discard(owner)
                self._pending.pop(owner, None)
                self._close_transfer(owner)
                self._direct_pending.pop(owner, None)
                self._jobs.pop(("direct", owner), None)
                self._wake_direct()
                for stream_id, stream in list(self._streams.items()):
                    if stream["owner"] == owner:
                        self._close_stream(stream_id)
                for job_key, job in list(self._jobs.items()):
                    if job_key[0] == "replay" and job[0] == 0 and not self._has_replay_demand(*job_key[1:]):
                        self._discard_job(job_key)
                if self._inflight is not None:
                    job = self._inflight["job"]
                    if job[2] == "replay" and not self._has_replay_demand(*job[3]):
                        self._inflight["cancel_followup"] = True
        except Exception as exc:  # return source failures without killing unrelated consumers
            self._reply(identity, error=f"Input {operation} failed: {exc}")

    def _close_transfer(self, owner):
        transfer = self._transfers.get(owner)
        if transfer is not None and transfer.get("loading"):
            transfer["closing"] = True
            return
        transfer = self._transfers.pop(owner, None)
        if transfer is not None:
            transfer["memory"].close()
            transfer["memory"].unlink()
            if self.storage.budget is not None:
                self.storage.budget.release_memory(transfer["size"], "input")
            self._metric("input_direct_active_bytes", -transfer["size"])

    def _wake_direct(self):
        for owner, args in self._direct_pending.items():
            self._enqueue(("direct", owner), 0, "direct", owner, *args)
        self._direct_pending.clear()

    def _direct_transfer(self, owner, slot, path, selection, job):
        started = perf_counter()
        dataset = self._entry((slot, path))
        shape = tuple(item.stop - item.start for item in selection)
        size = prod(shape) * dataset.dtype.itemsize
        if dataset.dtype.hasobject or not 0 < size <= self.policy.max_io_bytes:
            raise ValueError("Invalid or oversized direct-transfer slab")
        transfer = self._transfers.get(owner)
        if transfer is not None and transfer["busy"]:
            raise RuntimeError("Release the current direct transfer before reading")
        if transfer is not None and transfer["size"] < size:
            self._close_transfer(owner)
            transfer = None
        if transfer is None:
            # Reclaim idle slabs from other clients instead of retaining a pool
            # proportional to every pipeline that has ever run.
            for other, slab in list(self._transfers.items()):
                if not slab["busy"]:
                    self._close_transfer(other)
            budget = self.storage.budget
            if budget is not None and not budget.reserve_memory(size, "input"):
                self._release_planned_streams()
                if not budget.reserve_memory(size, "input"):
                    raise MemoryError("Insufficient RAM budget for direct input transfer")
            try:
                memory = SharedMemory(create=True, size=size)
            except BaseException:
                if budget is not None:
                    budget.release_memory(size, "input")
                raise
            transfer = self._transfers[owner] = dict(memory=memory, size=size, busy=False)
            self._metric("input_direct_allocations")
            self._metric("input_direct_active_bytes", size)
            self.storage.metrics.maximum("input_direct_peak_bytes", sum(item["size"] for item in self._transfers.values()))
        else:
            self._metric("input_direct_buffer_reuses")
        block = np.ndarray(shape, dtype=dataset.dtype, buffer=transfer["memory"].buf)
        transfer.update(busy=True, loading=True)
        self._transfer_serial += 1
        transfer["token"] = self._transfer_serial
        self._submit_io(job, self.aliases[slot, path], selection, block,
                        info=dict(owner=owner, shape=shape, started=started))

    def _record_replay_consumption(self, key, regions):
        replay = self._replays[key]
        self._metric("input_region_reuses", sum(r in replay.served for r in regions))
        consumed = tuple(r for r in regions if r in replay.prefetched and r not in replay.served)
        self._metric("input_prefetched_regions_consumed", len(consumed))
        if consumed:
            Logger.log_perf(f"input_prefetch_consumed dataset={key} regions={consumed} at={perf_counter():.6f}")
        replay.served.update(regions)

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
            self._record_replay_consumption(key, regions)
            Logger.log_perf(f"input_local_replay_mapping dataset={key} frames={first}:{stop} at={perf_counter():.6f}")
            del self._pending[identity]

    def _read_direct(self, dataset, destination, selection):
        self.reader.read_direct(dataset, destination, selection)

    def _wake_completion(self):
        try:
            self.requests.put_nowait((None, "io_complete", (), perf_counter()))
        except Full:
            pass  # Existing request traffic already wakes the control loop.

    def _submit_io(self, job, key, selection=None, destination=None, *, flush=None, copy_from=None, reserved=0, info=None):
        if self._inflight is not None:
            raise RuntimeError("Only one source operation may be in flight")
        self._io_serial += 1
        self._inflight = dict(token=self._io_serial, job=job, reserved=reserved,
                              issued_priority=job[0], **(info or {}))
        self.reader.commands.put_nowait(SourceJob(self._io_serial, key, selection, destination,
                                                  flush, copy_from, job[0], job[4]))

    def _cancel_running_prefetch(self, predicate):
        if self._inflight is not None:
            job = self._inflight["job"]
            if job[2] == "replay" and job[0] > 0 and predicate(job[3][0]):
                self._inflight["cancel_followup"] = True

    def _drain_completions(self):
        while True:
            try:
                result = self.reader.completed.get_nowait()
            except Empty:
                return
            if self._inflight is None or result.token != self._inflight["token"]:
                self._metric("input_stale_completions")
                continue
            info, self._inflight = self._inflight, None
            job = info["job"]
            priority, _, operation, args, _ = job
            if info["reserved"] and self.storage.budget is not None:
                self.storage.budget.release_memory(info["reserved"], "input")
            self._metric("input_io_operation_seconds", result.ended - result.started)
            local_copy = operation == "stream" and info.get("replay_key") is not None
            self._metric("input_replay_stream_copy_seconds" if local_copy else "input_source_operation_seconds",
                         result.ended - result.started)
            self._metric("input_completion_handoff_seconds", perf_counter() - result.ended)
            self._metric("input_io_jobs_completed")
            for source, size, first, last, small in result.reads:
                self._record_interval(self.io_intervals, (first, last))
                self._metric("input_hdf5_small_read_seconds" if small else "input_hdf5_read_seconds", last-first)
                self._metric("source_reads")
                self._metric("input_logical_bytes", size)
                metrics = self.dataset_metrics.setdefault(source, {"bytes": 0, "seconds": 0., "reads": 0})
                metrics["bytes"] += size
                metrics["seconds"] += last-first
                metrics["reads"] += 1
            if result.error:
                self._metric("input_source_operation_failures")
                self._intervals_dropped = True  # An unsuccessful native operation has incomplete interval coverage.
            if operation == "direct":
                owner = args[0]
                transfer = self._transfers.get(owner)
                if transfer is not None:
                    transfer["loading"] = False
                    if self._stop.is_set() or transfer.get("closing") or owner in self._cancelled_clients or result.error:
                        self._close_transfer(owner)
                    if result.error:
                        self._reply(owner, error=f"Input direct_read failed: {result.error}")
                    elif owner in self._transfers and not self._stop.is_set():
                        self._metric("input_direct_transfer_seconds", perf_counter()-info["started"])
                        self._reply(owner, (transfer["memory"].name, info["shape"], transfer["token"]))
                self._wake_direct()
            elif operation == "stream":
                self._read_stream(*args, result, info)
            elif operation == "small":
                key = args[0]
                owners = self._small_pending.pop(key, set())
                error = result.error or ("Input service cancelled" if self._stop.is_set() else None)
                if error is None:
                    try:
                        value = result.value
                        self.cache[key] = self.storage.store_value(value) if value.dtype.hasobject else self.storage.store(value, category="input")
                        self.source_reads[key] = 1
                    except Exception as exc:
                        error = f"Small input storage failed: {exc}"
                for owner in owners:
                    self._reply(owner, self.cache.get(key), error)
            else:
                self._complete_replay(result, info)
            self._fulfill()
            for owner, arguments in list(self._open_pending.items()):
                if self._stop.is_set():
                    self._reply(owner, error="Input service cancelled")
                    self._open_pending.pop(owner, None)
                    continue
                if owner in self._cancelled_clients:
                    self._open_pending.pop(owner, None)
                    continue
                try:
                    self._open_stream(owner, *arguments)
                except MemoryError as exc:
                    if any(stream.get("closing") for stream in self._streams.values()):
                        continue
                    self._reply(owner, error=str(exc))
                except Exception as exc:
                    self._reply(owner, error=str(exc))
                self._open_pending.pop(owner, None)

    def _run_job(self, job):
        priority, _, operation, args, queued = job
        self._metric("input_job_queue_seconds", perf_counter() - queued)
        if operation == "small":
            if args[0] not in self.cache:
                self._submit_io(job, args[0])
            return
        if operation == "direct":
            owner, slot, path, selection = args
            if owner in self._cancelled_clients:
                return
            if sum(item["busy"] for item in self._transfers.values()) >= self.policy.direct_buffers:
                self._direct_pending[owner] = slot, path, selection
            else:
                try:
                    self._direct_transfer(owner, slot, path, selection, job)
                except Exception as exc:
                    self._reply(owner, error=f"Input direct_read failed: {exc}")
            return
        if operation == "stream":
            try:
                self._prepare_stream(*args, job)
            except Exception as exc:
                stream = self._streams.get(args[0])
                if stream is not None:
                    stream["error"] = f"Input stream failed: {exc}"
                    self._stream_ready(args[0])
            return
        key, region = args
        replay = self._replays[key]
        if replay.error or replay.states[region] == "ready":
            return
        first = region * replay.frames
        stop = min(first + replay.frames, replay.array.shape[0])
        selection, destination = np.s_[first:stop], replay.array[first:stop]
        if destination.nbytes > self.policy.max_io_bytes:
            if region not in replay.tiles:
                tiles = iter(_tiles(destination.shape, self.policy.max_io_bytes, destination.dtype.itemsize,
                                    (first, *(0 for _ in destination.shape[1:]))))
                replay.tiles[region] = next(tiles), tiles
            selection = replay.tiles[region][0]
            destination = replay.array[selection]
        size, budget = destination.nbytes, self.storage.budget
        reserved = budget is None or budget.reserve_memory(size, "input")
        if not reserved and priority == 0:
            self._release_planned_streams()
            reserved = budget.reserve_memory(size, "input")
        if not reserved:
            if priority > 0:
                replay.states[region] = "not_requested"
                self._plan_retry = perf_counter() + 0.02
                self._metric("input_speculation_declined_memory")
            else:
                replay.error = f"Insufficient RAM budget for a {size:,}-byte input region"
                replay.states[region] = "failed"
                self._fulfill()
            return
        replay.states[region] = "loading"
        self._submit_io(job, key, selection, destination, flush=replay.array, reserved=size,
                        info=dict(working_bytes=size, staged=perf_counter()))

    def _complete_replay(self, result, info):
        priority, _, _, (key, region), _ = info["job"]
        replay = self._replays[key]
        if result.error:
            replay.tiles.pop(region, None)
            replay.error = f"Could not read {key} region {region}: {result.error}"
            replay.states[region] = "failed"
            self._metric("input_replays_failed")
            return
        if self._stop.is_set():
            return
        self._metric("input_replay_flush_seconds", result.flush_seconds)
        self._metric("input_replay_stage_seconds", result.ended-info["staged"])
        self._metric("input_replay_bytes", info["working_bytes"])
        if info["issued_priority"] > 0:
            self._metric("input_speculative_source_bytes", info["working_bytes"])
        if region in replay.tiles:
            _, tiles = replay.tiles[region]
            following = next(tiles, None)
            if following is not None:
                if info.get("cancel_followup") and (priority > 0 or not self._has_replay_demand(key, region)):
                    replay.tiles.pop(region)
                    replay.states[region] = "not_requested"
                    return
                replay.tiles[region] = following, tiles
                replay.states[region] = "queued"
                self._enqueue(("replay", key, region), priority, "replay", key, region)
                return
            del replay.tiles[region]
        replay.states[region] = "ready"
        if info["issued_priority"] > 0:
            replay.prefetched.add(region)
        following = next((r for r, state in enumerate(replay.states) if state == "not_requested"), None)
        Logger.log_perf(f"input_region_ready dataset={key} region={region} at={perf_counter():.6f} "
                        f"next_missing_region={following} eligibility=recheck_at_operation_boundary job={result.token}")
        self.source_reads.setdefault(key, 1)
        if sum(state == "ready" for state in replay.states) == 1:
            self._metric("input_first_region_seconds", perf_counter()-replay.created)
        if all(state == "ready" for state in replay.states):
            self.cache[key] = replay.reference
            self._metric("input_replays_completed")

    def _owned_stream(self, owner, stream_id):
        stream = self._streams[stream_id]
        if stream["owner"] != owner:
            raise RuntimeError("Input stream belongs to another client")
        return stream

    def _release_planned_streams(self, predicate=lambda key: True):
        for identity, stream in list(self._streams.items()):
            if stream["owner"] is None and predicate(stream["key"]):
                Logger.log_perf(f"input_cancellation stream={identity} dataset={stream['key']}")
                self._close_stream(identity)

    def _open_stream(self, owner, slot, path, frames, max_bytes, *, priority=1):
        dataset = self._entry((slot, path))
        if not dataset.shape or dataset.dtype.hasobject:
            raise TypeError("Streaming requires a nonscalar fixed-size dtype")
        frames = self.policy.stream_frames if frames is None else frames
        if isinstance(frames, bool) or not isinstance(frames, (int, np.integer)) or frames < 1:
            raise ValueError("frames must be a positive integer")
        frames = self._frames(dataset, int(frames))
        frame_bytes = prod(dataset.shape[1:]) * dataset.dtype.itemsize
        limit = self.policy.max_io_bytes if max_bytes is None else min(self.policy.max_io_bytes, max_bytes)
        if max_bytes is not None:
            if isinstance(max_bytes, bool) or not isinstance(max_bytes, (int, np.integer)) or max_bytes < 1:
                raise ValueError("max_bytes must be a positive integer")
            frames = min(frames, max(1, max_bytes // max(1, int(np.prod(dataset.shape[1:])) * dataset.dtype.itemsize)))
        if frame_bytes > limit:
            raise MemoryError("One stream frame exceeds the input byte limit")
        key = self.aliases[slot, path]
        if owner is not None:
            for identity, stream in list(self._streams.items()):
                if stream["owner"] is not None or stream["key"] != key:
                    continue
                # A pool may have shrunk under memory pressure. Smaller blocks
                # still meet the caller's upper bound and preserve lease semantics.
                if stream["frames"] <= frames and stream["frames"] * frame_bytes <= limit:
                    stream.update(owner=owner, priority=1)
                    for index, buffer in enumerate(stream["buffers"]):
                        if buffer["state"] == "filling":
                            self._enqueue(("stream", identity, index), 1, "stream", identity, index)
                    self._metric("input_planned_streams_adopted")
                    self._reply(owner, identity)
                    return
                self._close_stream(identity)
            if len(self._streams) >= self.policy.max_streams:
                self._release_planned_streams()
        if len(self._streams) >= self.policy.max_streams:
            raise MemoryError("Input stream capacity exhausted; close another stream first")
        budget = self.storage.budget
        count = min(self.policy.stream_buffers, self.policy.max_prefetch_jobs) if owner is None else self.policy.stream_buffers
        while True:
            shape = (frames, *dataset.shape[1:])
            size = int(np.prod(shape)) * dataset.dtype.itemsize
            if size == 0:
                raise ValueError("Cannot stream empty spatial dimensions")
            reserved = size * count
            if budget is None or budget.reserve_memory(reserved, "input"):
                break
            if owner is not None and any(stream["owner"] is None for stream in self._streams.values()):
                self._release_planned_streams()
                continue
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
            replay=self._backend((slot, path)) == "replay", priority=priority)
        if owner is None:
            self._stream_planned.add(key)
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
                self._enqueue(("stream", identity, index), stream["priority"], "stream", identity, index)

    def _prepare_stream(self, identity, index, job):
        stream = self._streams.get(identity)
        if stream is None or stream.get("closing") or stream.get("error"):
            return
        buffer = stream["buffers"][index]
        if buffer["state"] in ("ready", "in_use"):
            return
        key = stream["key"]
        replay = self._replay(key) if stream["replay"] or key in self._replays else None
        regions, source = (), None
        if replay is not None:
            regions = range(buffer["start"] // replay.frames, (buffer["stop"]-1) // replay.frames+1)
            missing = next((r for r in regions if replay.states[r] != "ready"), None)
            if missing is not None:
                if replay.error:
                    raise OSError(replay.error)
                self._queue_region(key, missing, 0)
                self._enqueue(("stream", identity, index), 2, "stream", identity, index)
                return
            source = replay.array[buffer["start"]:buffer["stop"]]
        block = np.ndarray((buffer["stop"]-buffer["start"], *stream["shape"][1:]),
                           dtype=stream["dtype"], buffer=buffer["memory"].buf)
        self._submit_io(job, key, np.s_[buffer["start"]:buffer["stop"]], block, copy_from=source,
                        info=dict(replay_key=key if replay is not None else None, regions=regions))

    def _read_stream(self, identity, index, result, info):
        stream = self._streams.get(identity)
        if stream is None:
            return
        if stream.get("closing") or self._stop.is_set():
            self._close_stream(identity)
            return
        buffer = stream["buffers"][index]
        if result.error:
            buffer["state"] = "failed"
            stream["error"] = f"Input stream read failed: {result.error}"
        else:
            if info["replay_key"] is not None:
                replay = self._replays[info["replay_key"]]
                replay.demanded.update(info["regions"])
                self._record_replay_consumption(info["replay_key"], info["regions"])
            buffer.update(state="ready", speculative=info["issued_priority"] >= 2, consumed=False)
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
            buffer["consumed"] = True
            if buffer.get("speculative"):
                Logger.log_perf(f"input_prefetch_consumed dataset={stream['key']} "
                                f"stream_frames={buffer['start']}:{buffer['stop']} at={perf_counter():.6f}")
                self._metric("input_buffer_prefetch_hits")
                self._metric("input_prefetched_bytes_consumed", (buffer["stop"] - buffer["start"])
                             * prod(stream["shape"][1:]) * stream["dtype"].itemsize)
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
        stream = self._streams[identity]
        if stream["waiting"] is not None:
            self._reply(stream["owner"], error="Input stream cancelled")
            stream["waiting"] = None
        stream["closing"] = True
        for index in range(len(stream["buffers"])):
            self._jobs.pop(("stream", identity, index), None)
        if self._inflight is not None and self._inflight["job"][2] == "stream" and self._inflight["job"][3][0] == identity:
            return  # Reader still owns a destination; completion retires the pool.
        self._streams.pop(identity)
        unused = sum((buffer["stop"] - buffer["start"]) * prod(stream["shape"][1:]) * stream["dtype"].itemsize
                     for buffer in stream["buffers"] if buffer.get("speculative") and not buffer.get("consumed"))
        self._metric("input_unused_prefetch_bytes", unused)
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
            while True:
                self._drain_completions()
                for _ in range(32):
                    try:
                        request = self.requests.get_nowait()
                    except Empty:
                        break
                    self._handle(request)
                self._drain_completions()
                if self._stop.is_set() and self._inflight is None:
                    break
                if not self._stop.is_set() and self._inflight is None:
                    if not self._jobs or min(job[0] for job in self._jobs.values()) > 2:
                        self._plan_prefetch()
                    if self._jobs:
                        self._set_idle()
                        key = min(self._jobs, key=lambda key: self._jobs[key][:2])
                        self._run_job(self._jobs.pop(key))
                        if self._inflight is None:
                            continue
                busy = self._inflight is not None
                if busy:
                    self._set_idle()
                else:
                    self._set_idle(self._idle_cause)
                idle = perf_counter()
                reason, cause = self._idle_reason, self._idle_cause
                try:
                    request = self.requests.get(timeout=0.1 if busy else 0.02)
                except Empty:
                    request = None
                if not busy:
                    waited = perf_counter()-idle
                    self._metric(f"input_idle_{reason}_seconds", waited)
                    if cause != reason:
                        self._metric(f"input_idle_{cause}_seconds", waited)
                if request is not None:
                    self._handle(request)
                if busy and not self.reader.thread.is_alive() and self.reader.completed.empty():
                    raise OSError("Input source reader exited without a completion")
        finally:
            self._set_idle()
            self.stopped.set()
            self.reader.close()  # No destination may be retired while it is being written.
            if self._inflight is not None:
                if self._inflight["reserved"] and self.storage.budget is not None:
                    self.storage.budget.release_memory(self._inflight["reserved"], "input")
                self._inflight = None
            for transfer in self._transfers.values():
                transfer["loading"] = False
            for owners in self._small_pending.values():
                for owner in owners:
                    self._reply(owner, error="Input service cancelled")
            self._small_pending.clear()
            for owner in self._open_pending:
                self._reply(owner, error="Input service cancelled")
            self._open_pending.clear()
            for identity in list(self._pending):
                self._reply(identity, error="Input service cancelled")
            self._pending.clear()
            self._jobs.clear()
            for owner in self._direct_pending:
                self._reply(owner, error="Input service cancelled")
            self._direct_pending.clear()
            for owner in list(self._transfers):
                self._close_transfer(owner)
            for identity in list(self._streams):
                stream = self._streams[identity]
                if stream["waiting"] is not None:
                    self._reply(stream["owner"], error="Input stream cancelled")
                self._close_stream(identity)
            for key, replay in self._replays.items():
                region_bytes = lambda r: (min(replay.frames, replay.array.shape[0] - r * replay.frames)
                                          * prod(replay.array.shape[1:]) * replay.array.dtype.itemsize)
                self._metric("input_prefetched_bytes_consumed", sum(region_bytes(r)
                             for r in replay.prefetched & replay.served))
                self._metric("input_prefetched_regions_unused", len(replay.prefetched - replay.served))
                Logger.log_perf(f"input_replay_usage dataset={key} prefetched={len(replay.prefetched)} "
                                f"consumed={len(replay.prefetched & replay.served)}")
                unused = sum(min(replay.frames, replay.array.shape[0] - r * replay.frames)
                             for r, state in enumerate(replay.states)
                             if state == "ready" and r not in replay.served)
                self._metric("input_unused_prefetch_bytes", unused * int(np.prod(replay.array.shape[1:]))
                             * replay.array.dtype.itemsize)
                if not all(state == "ready" for state in replay.states):
                    self._metric("input_replays_cancelled")
                replay.array._mmap.close()
            useful = []
            for owner, first, last in self.compute_intervals:
                useful.extend(_subtract_intervals((first, last),
                    [(a, b) for client, a, b in self.wait_intervals if client == owner]))
            self._metric("input_overlap_valid", int(not self._intervals_dropped))
            if not self._intervals_dropped:
                self._metric("input_io_compute_overlap_seconds", interval_overlap(self.io_intervals, useful))
            Logger.log_perf(f"input_datasets metrics={self.dataset_metrics}")
            Logger.log_perf(f"input_summary metrics={self.storage.metrics.snapshot()}")

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
