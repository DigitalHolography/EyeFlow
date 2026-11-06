"""Stable product identities and explicitly invocation-local resources."""

from dataclasses import dataclass
from collections.abc import Mapping
from types import MappingProxyType


@dataclass(frozen=True, slots=True)
class ProductHandle:
    id: str
    storage: object
    owner_id: str


@dataclass(frozen=True, slots=True)
class InvocationResource:
    """Own a local reader/replay/working resource, not a scientific value.

    Reads retain resource identity. Scientific results belong in ordinary
    handoffs; this wrapper does not provide a value snapshot or writable copy.
    Never eligible for IPC transport.
    """

    value: object


@dataclass(frozen=True, slots=True)
class ResourceReference:
    id: str


@dataclass(frozen=True, slots=True)
class ScientificValue:
    type: object
    fields: object


@dataclass(frozen=True, slots=True)
class StringArray:
    """Immutable string metadata; materialize only a private object-array view."""

    values: tuple
    shape: tuple
    dtype: object

    def __array__(self, dtype=None, copy=None):
        import numpy as np

        array = np.asarray(self.values, dtype=self.dtype if dtype is None else dtype).reshape(self.shape)
        array.flags.writeable = False
        return array


class FrozenMapping(Mapping):
    """Picklable read-only mapping, including recursively frozen metadata."""

    __slots__ = ("_data",)

    def __init__(self, values=()):
        object.__setattr__(self, "_data", MappingProxyType(dict(values)))

    def __setattr__(self, name, value):
        raise TypeError("Accepted product mappings are read-only.")

    def __iter__(self):
        return iter(self._data)

    def __len__(self):
        return len(self._data)

    def __getitem__(self, key):
        return self._data[key]

    def __reduce__(self):
        return type(self), (dict(self),)

    def __deepcopy__(self, memo):
        return self
