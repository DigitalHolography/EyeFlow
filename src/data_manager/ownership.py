"""Snapshot once, share read-only representations, and own local resources.

ProductStorage implements physical persistence. This manager gives structured
scientific values and reusable handles that policy without consumer copies.
It knows nothing about graph or invocation readiness.
"""

from collections.abc import Mapping
from dataclasses import fields, is_dataclass
from pathlib import Path
from uuid import uuid4

import numpy as np

from .handles import FrozenMapping, InvocationResource, ProductHandle, ResourceReference, ScientificValue, StringArray


def _envelope(value_type, values):
    """Rebuild existing data, without rerunning scientific __post_init__ work."""
    result = object.__new__(value_type)
    for name, value in values.items():
        object.__setattr__(result, name, value)
    return result


class DataManager:
    def __init__(self, storage, *, owner_id=None):
        self.storage = storage
        self.owner_id = owner_id or str(storage.directory.resolve())
        self.owner_root = Path(self.owner_id).resolve()
        self._values = {}
        self._backing = {}
        self._resources = {}

    def accept(self, value, *, transportable=False):
        snapshot = self._snapshot(value, {})
        if transportable and isinstance(snapshot, (Mapping, tuple, frozenset, ScientificValue, StringArray)):
            # A structured scientific value can contain local Python types or
            # extensive metadata. IPC carries one stable metadata reference.
            snapshot = self.storage.store_value(snapshot)
        handle = ProductHandle(uuid4().hex, snapshot, self.owner_id)
        return handle

    def _snapshot(self, value, memo):
        from pipeline_engine.data.storage import ArrayReference, ValueReference

        if isinstance(value, InvocationResource):
            key = uuid4().hex
            self._resources[key] = value.value
            return ResourceReference(key)
        if isinstance(value, (ArrayReference, ValueReference)):
            self.validate_storage(value)
            return value
        if isinstance(value, np.ndarray):
            cached = memo.get(id(value))
            if cached is not None and cached[0] is value:
                return cached[1]
            if value.dtype.hasobject:
                if not all(isinstance(item, (str, bytes)) for item in value.flat):
                    raise TypeError("Object arrays cannot be immutable scientific products.")
                return StringArray(tuple(value.flat), tuple(value.shape), value.dtype)
            # ProductStorage is the single snapshot implementation used by
            # handles, managed output records and legacy invocation transport.
            snapshot = self.storage.store(value)
            memo[id(value)] = value, snapshot
            return snapshot
        if isinstance(value, Mapping):
            return FrozenMapping({key: self._snapshot(item, memo) for key, item in value.items()})
        if isinstance(value, (list, tuple)):
            return tuple(self._snapshot(item, memo) for item in value)
        if isinstance(value, (set, frozenset)):
            return frozenset(self._snapshot(item, memo) for item in value)
        if is_dataclass(value) and not isinstance(value, type):
            return ScientificValue(type(value), FrozenMapping({f.name: self._snapshot(getattr(value, f.name), memo)
                                                               for f in fields(value)}))
        if value is None or isinstance(value, (str, bytes, int, float, complex, bool, np.generic, np.dtype, Path)):
            return value
        raise TypeError(f"Use InvocationResource for invocation-local {type(value).__name__} values.")

    def read(self, handle):
        if handle.id not in self._values:
            self._values[handle.id] = self._resolve(handle.storage)
        if isinstance(handle.storage, ResourceReference):
            return self._values[handle.id]
        # Mutable scientific dataclass envelopes are private to a reader; all
        # underlying arrays remain shared. No consumer can edit the canonical value.
        return self._consumer_view(self._values[handle.id])

    def _consumer_view(self, value):
        if isinstance(value, Mapping):
            return FrozenMapping({key: self._consumer_view(item) for key, item in value.items()})
        if isinstance(value, tuple):
            return tuple(self._consumer_view(item) for item in value)
        if is_dataclass(value) and not isinstance(value, type):
            return _envelope(type(value), {f.name: self._consumer_view(getattr(value, f.name))
                                           for f in fields(value)})
        return value

    def _resolve(self, value):
        from pipeline_engine.data.storage import ArrayReference, ValueReference

        if isinstance(value, ResourceReference):
            return self._resources[value.id]
        if isinstance(value, ScientificValue):
            return _envelope(value.type, {key: self._resolve(item) for key, item in value.fields.items()})
        if isinstance(value, (ArrayReference, ValueReference)):
            if value not in self._backing:
                loaded = self.storage.read(value)
                self._backing[value] = self._resolve(loaded) if isinstance(value, ValueReference) else loaded
            return self._backing[value]
        if isinstance(value, Mapping):
            return FrozenMapping({key: self._resolve(item) for key, item in value.items()})
        if isinstance(value, tuple):
            return tuple(self._resolve(item) for item in value)
        if is_dataclass(value) and not isinstance(value, type):
            return _envelope(type(value), {f.name: self._resolve(getattr(value, f.name))
                                           for f in fields(value)})
        return value

    def writable(self, handle):
        def copy(value):
            if isinstance(value, np.ndarray):
                self.storage.metrics.add("array_copies", 1)
                self.storage.metrics.add("array_copied_bytes", value.nbytes)
                return value.copy()
            if isinstance(value, Mapping):
                return {key: copy(item) for key, item in value.items()}
            if isinstance(value, tuple):
                return tuple(copy(item) for item in value)
            if is_dataclass(value) and not isinstance(value, type):
                return _envelope(type(value), {f.name: copy(getattr(value, f.name))
                                               for f in fields(value)})
            return value
        return copy(self.read(handle))

    def acquire(self, handle):
        """Validate run-owned backing before global availability, without copying."""
        if handle.owner_id != self.owner_id:
            raise ValueError("Product belongs to a different run storage owner.")
        self.validate_storage(handle.storage, product_paths=True)
        return handle

    def validate_storage(self, value, *, product_paths=False):
        from pipeline_engine.data.storage import ArrayReference, ValueReference

        if isinstance(value, ResourceReference):
            raise ValueError("Invocation-local resources cannot cross process boundaries.")
        if isinstance(value, (ArrayReference, ValueReference)):
            path = Path(value.path).resolve()
            if not path.is_relative_to(self.owner_root) or not path.is_file():
                raise ValueError(f"Invalid run-owned storage reference: {path}")
            if isinstance(value, ArrayReference):
                array = value.open()
                try:
                    if tuple(array.shape) != value.shape or array.dtype != value.dtype:
                        raise ValueError("Product storage shape/dtype mismatch.")
                finally:
                    array._mmap.close()
            else:
                self.validate_storage(value.open(), product_paths=product_paths)
        elif isinstance(value, Path) and product_paths:
            path = value.resolve()
            if not path.is_relative_to(self.owner_root) or not path.exists():
                raise ValueError("Product file data must use stable run-owned storage references.")
        elif isinstance(value, Mapping):
            for item in value.values():
                self.validate_storage(item, product_paths=product_paths)
        elif isinstance(value, (list, tuple, set, frozenset)):
            for item in value:
                self.validate_storage(item, product_paths=product_paths)
        elif is_dataclass(value) and not isinstance(value, type):
            for f in fields(value):
                self.validate_storage(getattr(value, f.name), product_paths=product_paths)

    def close_resources(self):
        errors = []
        try:
            for value in reversed(tuple(self._resources.values())):
                close = getattr(value, "close", None)
                if close is not None:
                    try:
                        close()
                    except Exception as exc:
                        errors.append(exc)
        finally:
            self._resources.clear()
            self._values.clear()
            self._backing.clear()
        if errors:
            raise errors[0]
