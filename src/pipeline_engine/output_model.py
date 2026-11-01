"""In-memory output structure with local immutable payloads and stage journals.

Checkpoints track edited metadata without copying dataset bytes. A failed stage
restores only its touched paths; the I/O writers own final serialization.
"""

from __future__ import annotations

from collections.abc import MutableMapping
from contextlib import contextmanager
from copy import deepcopy
from dataclasses import dataclass, field, replace
import numpy as np

from .storage import ArrayReference, ValueReference


@dataclass
class Record:
    kind: str = "group"
    value: object = None
    attrs: dict = field(default_factory=dict)
    attr_dtypes: dict = field(default_factory=dict)
    options: dict = field(default_factory=dict)

    def snapshot(self):
        return replace(
            self,
            attrs=deepcopy(self.attrs),
            attr_dtypes=dict(self.attr_dtypes),
            options=dict(self.options),
        )


class Attributes(MutableMapping):
    def __init__(self, node):
        self.node = node

    @property
    def record(self):
        return self.node.file.records[self.node.path]

    def __iter__(self):
        return iter(self.record.attrs)

    def __len__(self):
        return len(self.record.attrs)

    def __getitem__(self, key):
        return deepcopy(self.record.attrs[key])

    def __setitem__(self, key, value):
        self.node.file._touch(self.node.path)
        self.record.attrs[key] = deepcopy(value)
        self.record.attr_dtypes.pop(key, None)

    def __delitem__(self, key):
        self.node.file._touch(self.node.path)
        del self.record.attrs[key]
        self.record.attr_dtypes.pop(key, None)

    def create(self, key, value, dtype=None):
        self[key] = value
        if dtype is not None:
            self.record.attr_dtypes[key] = dtype


class OutputNode:
    def __init__(self, file, path):
        self.file, self.path = file, path

    @property
    def name(self):
        return "/" + self.path

    @property
    def parent(self):
        return self.file[self.path.rpartition("/")[0]]

    @property
    def attrs(self):
        return Attributes(self)


class OutputDataset(OutputNode):
    _eyeflow_dataset = True

    @property
    def _value(self):
        return self.file.storage.read(self.file.records[self.path].value)

    @property
    def shape(self):
        return self.file.records[self.path].value.shape

    @property
    def dtype(self):
        return self.file.records[self.path].value.dtype

    @property
    def ndim(self):
        return len(self.shape)

    @property
    def size(self):
        return int(np.prod(self.shape))

    def __getitem__(self, key):
        return self._value[key]

    def __array__(self, dtype=None, copy=None):
        array = np.asarray(self._value, dtype=dtype)
        return array.copy() if copy else array


class OutputGroup(OutputNode):
    _eyeflow_group = True

    def _path(self, key):
        key = str(key).replace("\\", "/")
        return (
            key.strip("/")
            if key.startswith("/")
            else "/".join(part for part in (self.path, key.strip("/")) if part)
        )

    def get(self, key, default=None):
        path = self._path(key)
        record = self.file.records.get(path)
        if record is None:
            return default
        return (OutputDataset if record.kind == "dataset" else OutputGroup)(self.file, path)

    def __getitem__(self, key):
        value = self.get(key)
        if value is None:
            raise KeyError(key)
        return value

    def __contains__(self, key):
        return self.get(key) is not None

    def __delitem__(self, key):
        path = self._path(key)
        for child in list(self.file.records):
            if child == path or child.startswith(path + "/"):
                self.file._touch(child)
                del self.file.records[child]

    def keys(self):
        prefix = self.path + "/" if self.path else ""
        return sorted(
            {
                path[len(prefix) :].split("/")[0]
                for path in self.file.records
                if path.startswith(prefix) and path != self.path
            }
        )

    def items(self):
        return [(key, self[key]) for key in self.keys()]

    def __iter__(self):
        return iter(self.keys())

    def require_group(self, name):
        path = self._path(name)
        for index in range(1, len(path.split("/")) + 1):
            parent = "/".join(path.split("/")[:index])
            if parent not in self.file.records:
                self.file._touch(parent)
                self.file.records[parent] = Record()
            elif self.file.records[parent].kind != "group":
                raise ValueError(f"Dataset blocks group {parent}")
        return self.file[path]

    def create_group(self, name):
        if name in self:
            raise ValueError(f"Output group already exists: {name}")
        return self.require_group(name)

    def create_dataset(self, name, data=None, shape=None, dtype=None, **options):
        if data is None:
            data = np.full(shape, options.get("fillvalue", 0), dtype=dtype)
        elif dtype is not None:
            data = np.asarray(data, dtype=dtype)
        return self.store_dataset(name, data, options)

    def store_dataset(self, name, data, options=None):
        path = self._path(name)
        if not isinstance(data, (ArrayReference, ValueReference)):
            array = np.asarray(data)
            if array.dtype.kind == "O" and not all(
                isinstance(item, (str, bytes)) for item in array.flat
            ):
                raise TypeError(
                    f"Could not write HDF5 dataset '{name}': unsupported object payload"
                )
            if array.dtype.kind not in "biufcSUO":
                raise TypeError(f"Unsupported HDF5 dtype: {array.dtype}")
        self.file.require_group(path.rpartition("/")[0])
        value = self.file.storage.store(data)
        self.file._touch(path)
        self.file.records[path] = Record("dataset", value, options=dict(options or {}))
        return self.file[path]


class OutputModel(OutputGroup):
    def __init__(self, storage, filename):
        super().__init__(self, "")
        self.storage = storage
        self.filename = str(filename)
        self.records = {"": Record()}
        self._checkpoints = []

    def _touch(self, path):
        for checkpoint in self._checkpoints:
            if path not in checkpoint.before:
                old = self.records.get(path)
                checkpoint.before[path] = old.snapshot() if old is not None else None

    def flush(self):
        pass

    @contextmanager
    def checkpoint(self):
        checkpoint = ProductCheckpoint(self)
        self._checkpoints.append(checkpoint)
        try:
            yield checkpoint
        except BaseException:
            checkpoint.rollback()
            raise
        finally:
            self._checkpoints.remove(checkpoint)

    def visititems(self, callback):
        for path in sorted(self.records):
            if path:
                callback(path, self[path])

    def visit(self, callback):
        for path in sorted(self.records):
            if path:
                callback(path)

    def copy(self, source, target, name=None):
        """Materialize a private export snapshot without opening the final file."""
        from input_output.writers.h5_publication import copy_output_node

        copy_output_node(self, source, target, name)

    def validate_merge(self, records):
        """Reject conflicts before accepting any invocation records."""
        for path, record in records.items():
            if path != path.strip("/") or "\\" in path:
                raise ValueError(f"Invalid output path: {path!r}")
            if record.kind not in {"group", "dataset"}:
                raise ValueError(f"Invalid output record kind: {record.kind!r}")
            current = self.records.get(path)
            if current and (record.kind == "dataset" or current.kind != record.kind):
                raise ValueError(f"Duplicate pipeline output: {path}")
            parent = path.rpartition("/")[0]
            while path and parent:
                ancestor = records.get(parent, self.records.get(parent))
                if ancestor is None or ancestor.kind != "group":
                    raise ValueError(f"Output parent must be a group: {parent}")
                parent = parent.rpartition("/")[0]

    def merge(self, records):
        """Accept one invocation atomically; only group attributes can be shared."""
        self.validate_merge(records)
        with self.checkpoint():
            for path, record in records.items():
                current = self.records.get(path)
                self._touch(path)
                if current:
                    current.attrs.update(deepcopy(record.attrs))
                    for key in record.attrs:
                        current.attr_dtypes.pop(key, None)
                    current.attr_dtypes.update(record.attr_dtypes)
                else:
                    adopted = record.snapshot()
                    if adopted.kind == "dataset":
                        adopted.value = self.storage.store(adopted.value)
                    self.records[path] = adopted

    def committed_records(self):
        """Snapshot accepted records without exposing unfinished checkpoint edits."""
        records = {path: record.snapshot() for path, record in self.records.items()}
        for checkpoint in reversed(self._checkpoints):
            for path, old in checkpoint.before.items():
                if old is None:
                    records.pop(path, None)
                else:
                    records[path] = old.snapshot()
        return records

    def publish(self):
        """Publish through the shared I/O policy (also used by standalone callers)."""
        from input_output.writers.h5_publication import publish_output_model

        publish_output_model(self)


class ProductCheckpoint:
    def __init__(self, model):
        self.model = model
        self.before = {}

    def commit(self, paths=()):
        committed = {""}
        for path in paths:
            parts = path.strip("/").split("/")
            committed.update("/".join(parts[:i]) for i in range(1, len(parts) + 1))
        for path in committed:
            self.before.pop(path, None)

    def rollback(self):
        for path, old in self.before.items():
            if old is None:
                self.model.records.pop(path, None)
            else:
                self.model.records[path] = old
