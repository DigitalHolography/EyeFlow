"""One coordinator-owned NAS reader with reusable local read-only datasets."""

from __future__ import annotations

from dataclasses import dataclass
from threading import RLock, Thread

import numpy as np

from input_output.readers.nodes import is_dataset

from .storage import ArrayReference, ValueReference


@dataclass
class SourceEntry:
    attrs: dict
    shape: tuple | None = None
    dtype: object = None
    identity: str | None = None


class InputClient:
    def __init__(self, identity, requests, responses):
        self.identity, self.requests, self.responses = identity, requests, responses
        self.cache = {}

    def read(self, slot, path):
        key = (slot, path)
        if key not in self.cache:
            self.requests.put((self.identity, slot, path))
            value, error = self.responses.get()
            if error:
                raise OSError(error)
            array = value.open() if isinstance(value, (ArrayReference, ValueReference)) else value
            if isinstance(array, np.ndarray):
                array.flags.writeable = False
            self.cache[key] = array
        return self.cache[key]

    def close(self):
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
        return self.file.client.read(self.file.slot, self.path)

    def __getitem__(self, selection):
        return self._value[selection]

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
        return self.file.get("/".join((self.path, str(path).strip("/"))), default)

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


class InputService:
    """Stage only requested datasets, in bounded reads, once per input run.

    The complete requested dataset is replayable locally. This avoids repeated
    NAS reads by pipelines with different slicing patterns without loading full
    videos into RAM. Unrequested HD datasets are never copied.
    """

    def __init__(self, sources, storage, mp_context):
        self.sources, self.storage = sources, storage
        self.requests = mp_context.Queue(maxsize=32)
        self.responses = {}
        self._context = mp_context
        self._lock = RLock()
        self.cache = {}
        self.metadata = {}
        self.source_reads = {}
        self.aliases = {}
        for slot, source in sources.items():
            if source is None:
                continue
            entries = {"": SourceEntry(dict(source.attrs))}
            identities = {}

            def walk(
                group, prefix="", ancestors=(), *, entries=entries, identities=identities, slot=slot
            ):
                for name in group:
                    node = group.get(name)
                    if node is None:
                        continue
                    path = "/".join(part for part in (prefix, name) if part)
                    identity = identities.setdefault(node.id, path)
                    dataset = is_dataset(node)
                    entries[path] = SourceEntry(
                        dict(node.attrs),
                        tuple(node.shape) if dataset else None,
                        node.dtype if dataset else None,
                        identity,
                    )
                    if dataset:
                        self.aliases[slot, path] = slot, identity
                    elif node.id not in ancestors:
                        walk(node, path, (*ancestors, node.id))

            walk(source, ancestors=(source.id,))
            self.metadata[slot] = (source.filename, entries)
        self.thread = Thread(target=self._serve, name="EyeFlowInputs", daemon=True)
        self.thread.start()

    def client(self):
        with self._lock:
            identity = len(self.responses)
            response = self._context.Queue(maxsize=1)
            self.responses[identity] = response
        return InputClient(identity, self.requests, response)

    def files(self, client):
        return {
            slot: InputFile(slot, *self.metadata[slot], client) if slot in self.metadata else None
            for slot in self.sources
        }

    def _serve(self):
        while True:
            request = self.requests.get()
            if request is None:
                return
            identity, slot, path = request
            try:
                key = self.aliases[slot, path]
                if key not in self.cache:
                    dataset = self.sources[slot][path]
                    if dataset.dtype.kind in "biufc" and dataset.nbytes >= 64 * 1024:
                        # store_array slices the source rather than materializing it.
                        value = self.storage.store_array(dataset)
                    else:
                        value = self.storage.store(dataset[()], category="input")
                    if isinstance(value, np.ndarray) and (
                        value.nbytes >= 4096 or value.dtype.kind == "O"
                    ):
                        if value.dtype.kind != "O":
                            value = self.storage.store_array(value)
                        else:
                            value = self.storage.store_value(value)
                    self.cache[key] = value
                    self.source_reads[key] = self.source_reads.get(key, 0) + 1
                response = self.cache[key], None
            except Exception as exc:  # noqa: BLE001 - return source failures to the requesting invocation
                response = None, f"Could not stage {slot}:{path}: {exc}"
            self.responses[identity].put(response)

    def close(self):
        self.requests.put(None)
        self.thread.join()
        for queue in [self.requests, *self.responses.values()]:
            queue.close()
            queue.join_thread()
