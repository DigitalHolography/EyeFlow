"""In-memory output structure with local immutable payloads and stage journals.

Checkpoints track edited metadata without copying dataset bytes. A failed stage
restores only its touched paths; the I/O writers own final serialization.
"""

from __future__ import annotations

from collections.abc import MutableMapping
from contextlib import contextmanager
from copy import deepcopy
from dataclasses import dataclass, field, replace
from types import MappingProxyType
import numpy as np

from data_manager.handles import FrozenMapping, StringArray
from pipeline_engine.data.storage import ArrayReference, ValueReference, H5Array, H5ArrayReference, _H5ArrayBuilder


class OutputError(ValueError):
    """A rejected managed-output operation."""


class OutputImmutableError(OutputError):
    pass


class OutputOwnershipError(OutputError):
    pass


class OutputPathConflictError(OutputError):
    pass


def _metadata(value):
    if isinstance(value, dict):
        return FrozenMapping({key: _metadata(item) for key, item in value.items()})
    if isinstance(value, (list, tuple)):
        return tuple(_metadata(item) for item in value)
    if isinstance(value, np.ndarray):
        if value.dtype.hasobject:
            return StringArray(tuple(value.flat), value.shape, value.dtype)
        return np.frombuffer(value.tobytes(), dtype=value.dtype).reshape(value.shape)
    return value


@dataclass(frozen=True)
class Record:
    kind: str = "group"
    value: object = None
    attrs: dict = field(default_factory=dict)
    attr_dtypes: dict = field(default_factory=dict)
    options: dict = field(default_factory=dict)
    owner: str | None = None
    committed: bool = False
    metadata_owner: str | None = None

    def __post_init__(self):
        for key in ("attrs", "attr_dtypes", "options"):
            value = getattr(self, key)
            if not isinstance(value, FrozenMapping):
                object.__setattr__(self, key, _metadata(dict(value)))

    def snapshot(self):
        return self


class Attributes(MutableMapping):
    def __init__(self, node):
        self.node = node

    @property
    def record(self):
        return self.node.file._records[self.node.path]

    def __iter__(self):
        return iter(self.record.attrs)

    def __len__(self):
        return len(self.record.attrs)

    def __getitem__(self, key):
        value = self.record.attrs[key]
        return np.asarray(value).copy() if isinstance(value, StringArray) else deepcopy(value)

    def __setitem__(self, key, value):
        self.node.file._check_attributes(self.node.path, key)
        attrs, dtypes = dict(self.record.attrs), dict(self.record.attr_dtypes)
        attrs[key] = value
        dtypes.pop(key, None)
        self.node.file._touch(self.node.path)
        self.node.file._records[self.node.path] = replace(self.record, attrs=attrs, attr_dtypes=dtypes)

    def __delitem__(self, key):
        self.node.file._check_attributes(self.node.path, key)
        attrs, dtypes = dict(self.record.attrs), dict(self.record.attr_dtypes)
        del attrs[key]
        dtypes.pop(key, None)
        self.node.file._touch(self.node.path)
        self.node.file._records[self.node.path] = replace(self.record, attrs=attrs, attr_dtypes=dtypes)

    def create(self, key, value, dtype=None):
        self[key] = value
        if dtype is not None:
            self.node.file._records[self.node.path] = replace(
                self.record, attr_dtypes={**self.record.attr_dtypes, key: dtype})


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

    def __setitem__(self, key, value):
        self.file._check_dataset_write(self.path)
        record = self.file._records[self.path]
        if isinstance(record.value, H5ArrayReference):
            backing = self.file.storage.updated_h5_slice(record.value, key, value)
        else:
            array = np.array(self._value, copy=True)
            array[key] = value
            backing = self.file.storage.store(array)
        self.file._touch(self.path)
        self.file._records[self.path] = replace(record, value=backing)

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
        if not path:
            raise OutputPathConflictError("Cannot delete the output root.")
        record = self.file._records.get(path)
        if record is None:
            raise KeyError(path)
        children = [child for child in self.file._records if child == path or child.startswith(path + "/")]
        if record.kind == "group" and any(self.file._records[child].kind == "dataset"
                                           and self.file._records[child].committed for child in children):
            raise OutputImmutableError(f"Cannot delete group '{path}' containing committed datasets.")
        for child in children:
            if self.file._records[child].kind == "dataset":
                self.file._check_dataset_write(child)
        if record.kind == "group":
            self.file._check_attributes(path)
        for child in children:
            self.file._touch(child)
            del self.file._records[child]
            self.file._versions.pop(child, None)

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
        if not path:
            return self.file
        parents = ["/".join(path.split("/")[:index]) for index in range(1, len(path.split("/")) + 1)]
        for parent in parents:
            if parent in self.file._records and self.file._records[parent].kind != "group":
                raise OutputPathConflictError(f"Output path '{path}' conflicts with dataset '{parent}'.")
        for parent in parents:
            if parent not in self.file._records:
                self.file._touch(parent)
                self.file._records[parent] = Record(metadata_owner=self.file.current_owner)
        return self.file[path]

    def create_group(self, name):
        if name in self:
            raise OutputPathConflictError(f"Output group already exists: {name}")
        return self.require_group(name)

    def create_dataset(self, name, data=None, shape=None, dtype=None, **options):
        if data is None:
            fill = options.get("fillvalue", 0)
            dtype = np.dtype(dtype if dtype is not None else np.asarray(fill).dtype)
            if dtype.kind in "buifc":
                with self.file.build_arrays() as allocate:
                    data = allocate(shape, dtype=dtype, fillvalue=fill, h5_options=options)
            else:
                data = np.full(shape, fill, dtype=dtype)
        elif dtype is not None:
            data = np.asarray(data, dtype=dtype)
        return self.store_dataset(name, data, options)

    def store_dataset(self, name, data, options=None):
        path = self._path(name)
        self.file._check_dataset_write(path)
        if not isinstance(data, (ArrayReference, ValueReference, H5ArrayReference, H5Array)):
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
        self.file._records[path] = Record("dataset", value, options=dict(options or {}), owner=self.file.current_owner)
        return self.file[path]


class OutputModel(OutputGroup):
    def __init__(self, storage, filename):
        super().__init__(self, "")
        self.storage = storage
        self.filename = str(filename)
        self._records = {"": Record(committed=True, metadata_owner="runtime")}
        self.current_owner = "runtime"
        self._checkpoints = []
        self._versions = {"": 0}
        self._generation = 0

    @contextmanager
    def build_arrays(self):
        """Construct bounded HDF5 payloads; acceptance still uses our checkpoint.

        No output record exists until the finished value is handed off. Failure
        closes and removes every unfinished payload, without publishing a path.
        """
        from input_output.writers.h5_storage import _DEFAULT_GZIP_LEVEL
        from pipeline_engine.resources import ensure_memory_available

        arrays = []

        def allocate(shape, *, h5_options=None, dtype=np.float32, fillvalue=np.nan):
            ensure_memory_available(16 * 2**20)
            shape = (int(shape),) if np.isscalar(shape) else tuple(shape)
            options = dict(h5_options or {})
            options.pop("sparse_nan_chunks", None)
            options.setdefault("chunks", True)
            options.setdefault("compression", "gzip")
            if options["compression"] == "gzip":
                options.setdefault("compression_opts", _DEFAULT_GZIP_LEVEL)
            options.setdefault("shuffle", True)
            options["fillvalue"] = np.asarray(fillvalue, dtype=dtype)[()]
            options["dtype"] = dtype
            if not shape:
                for key in ("chunks", "compression", "compression_opts", "shuffle"):
                    options.pop(key, None)
            array = _H5ArrayBuilder(self.storage, tuple(shape), options)
            arrays.append(array)
            return array

        try:
            yield allocate
            for array in arrays:
                array.close(success=True)
        except BaseException:
            for array in arrays:
                array.close(success=False)
            raise

    @property
    def records(self):
        return MappingProxyType(self._records)

    @contextmanager
    def writing_as(self, owner):
        previous = self.current_owner
        self.current_owner = owner
        try:
            yield self
        finally:
            self.current_owner = previous

    def _check_dataset_write(self, path):
        record = self._records.get(path)
        if record is None:
            return
        if record.kind != "dataset":
            raise OutputPathConflictError(f"Cannot create dataset at existing group '{path}'.")
        if record.owner != self.current_owner:
            raise OutputOwnershipError(
                f"Pipeline '{self.current_owner}' cannot write dataset '{path}' owned by pipeline '{record.owner}'.")
        if record.committed:
            raise OutputImmutableError(f"Cannot modify committed dataset '{path}'.")

    def _check_attributes(self, path, key=None):
        record = self._records[path]
        if record.kind == "dataset":
            self._check_dataset_write(path)
            return
        owner = ("stage_journal" if not path and key in
                 {"stage_status", "stage_failures", "stage_skips", "stage_products", "stage_timings"} else record.metadata_owner)
        if owner != self.current_owner:
            raise OutputOwnershipError(
                f"Pipeline '{self.current_owner}' cannot modify metadata of group '{path or '/'}' owned by '{owner}'.")

    def designate_metadata(self, path, owner):
        if self.current_owner != "runtime":
            raise OutputOwnershipError("Only the runtime may designate group metadata ownership.")
        group = self.require_group(path)
        self._touch(group.path)
        self._records[group.path] = replace(self._records[group.path], metadata_owner=owner)

    def _touch(self, path):
        for checkpoint in self._checkpoints:
            if path not in checkpoint.before:
                old = self.records.get(path)
                checkpoint.before[path] = old.snapshot() if old is not None else None
                checkpoint.versions[path] = self._versions.get(path)
        self._generation += 1
        self._versions[path] = self._generation

    def flush(self):
        pass

    @contextmanager
    def checkpoint(self, *, owner=None):
        if owner is not None:
            with self.writing_as(owner), self.checkpoint() as checkpoint:
                yield checkpoint
            return
        checkpoint = ProductCheckpoint(self)
        self._checkpoints.append(checkpoint)
        try:
            yield checkpoint
        except BaseException:
            checkpoint.rollback()
            raise
        else:
            if len(self._checkpoints) == 1:
                checkpoint.commit(tuple(checkpoint.before))
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

    def validate_merge(self, records, *, owner=None):
        """Reject conflicts before accepting any invocation records."""
        if owner is not None:
            with self.writing_as(owner):
                return self.validate_merge(records)
        for path, record in records.items():
            if path != path.strip("/") or "\\" in path:
                raise ValueError(f"Invalid output path: {path!r}")
            if record.kind not in {"group", "dataset"}:
                raise ValueError(f"Invalid output record kind: {record.kind!r}")
            current = self.records.get(path)
            if current and current.kind != record.kind:
                raise OutputPathConflictError(f"Output path '{path}' has a conflicting dataset/group kind.")
            if record.kind == "dataset":
                self._check_dataset_write(path)
                if record.owner is not None and record.owner != self.current_owner:
                    raise OutputOwnershipError(f"Invalid producing pipeline for dataset '{path}'.")
            elif current and record.attrs:
                for key in record.attrs:
                    self._check_attributes(path, key)
            parent = path.rpartition("/")[0]
            while path and parent:
                ancestor = records.get(parent, self.records.get(parent))
                if ancestor is not None and ancestor.kind != "group":
                    raise OutputPathConflictError(f"Output path '{path}' conflicts with dataset '{parent}'.")
                parent = parent.rpartition("/")[0]

    def merge(self, records, *, owner=None):
        """Accept one invocation atomically; only group attributes can be shared."""
        if owner is not None:
            with self.writing_as(owner):
                return self.merge(records)
        self.validate_merge(records)
        with self.checkpoint():
            for path, record in records.items():
                current = self.records.get(path)
                if current and record.kind == "group" and not record.attrs:
                    # An existing structural parent is shared, not rewritten.
                    continue
                self._touch(path)
                if current and record.kind == "group":
                    attrs, dtypes = dict(current.attrs), dict(current.attr_dtypes)
                    attrs.update(record.attrs)
                    for key in record.attrs:
                        dtypes.pop(key, None)
                    dtypes.update(record.attr_dtypes)
                    self._records[path] = replace(current, attrs=attrs, attr_dtypes=dtypes)
                else:
                    adopted = record.snapshot()
                    if adopted.kind == "dataset":
                        self.require_group(path.rpartition("/")[0])
                        adopted = replace(adopted, value=self.storage.store(adopted.value),
                                          owner=self.current_owner, committed=False)
                    else:
                        adopted = replace(adopted, metadata_owner=self.current_owner, committed=False)
                    self._records[path] = adopted

    def committed_records(self):
        """Snapshot accepted records without exposing unfinished checkpoint edits."""
        records = {path: record.snapshot() for path, record in self.records.items()}
        for checkpoint in reversed(self._checkpoints):
            for path, old in checkpoint.before.items():
                if old is None:
                    records.pop(path, None)
                else:
                    records[path] = old.snapshot()
        return {path: record for path, record in records.items() if record.committed}

    def committed_delta(self, cursor):
        """Return only changed accepted records and an invocation-owned cursor."""
        versions = dict(self._versions)
        overrides = {}
        for checkpoint in reversed(self._checkpoints):
            for path, old in checkpoint.before.items():
                version = checkpoint.versions[path]
                if old is None:
                    versions.pop(path, None)
                else:
                    versions[path] = version
                overrides[path] = old
        records = {
            path: overrides.get(path, self.records.get(path)).snapshot()
            for path, version in versions.items() if cursor.get(path) != version
            and overrides.get(path, self.records.get(path)).committed
        }
        # The root transports this invocation's diagnostics even for an empty
        # numerical delta; it never contains another invocation's worker tree.
        records.setdefault("", overrides.get("", self.records[""]).snapshot())
        return records, {path: version for path, version in versions.items()
                         if overrides.get(path, self.records.get(path)).committed}

    def publish(self):
        """Publish through the shared I/O policy (also used by standalone callers)."""
        from input_output.writers.h5_publication import publish_output_model

        if self._checkpoints:
            raise RuntimeError("Cannot publish an active output transaction.")
        ProductCheckpoint(self).commit(tuple(self._records))
        publish_output_model(self)


class ProductCheckpoint:
    def __init__(self, model):
        self.model = model
        self.before = {}
        self.versions = {}

    def commit(self, paths=()):
        committed = {""}
        for path in paths:
            parts = path.strip("/").split("/")
            committed.update("/".join(parts[:i]) for i in range(1, len(parts) + 1))
        for path in committed:
            record = self.model._records.get(path)
            if record is not None and not record.committed:
                self.model._records[path] = replace(record, committed=True)
            self.before.pop(path, None)
            self.versions.pop(path, None)

    def rollback(self):
        for path, old in self.before.items():
            if old is None:
                self.model._records.pop(path, None)
                self.model._versions.pop(path, None)
            else:
                self.model._records[path] = old
                self.model._versions[path] = self.versions[path]
