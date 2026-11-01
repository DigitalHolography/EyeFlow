"""Immutable payload storage with explicit mapping lifetime."""

from __future__ import annotations

import shutil
import weakref
from dataclasses import dataclass
from pathlib import Path
from threading import RLock
from uuid import uuid4

import numpy as np
from scratch_location import temporary_directory


@dataclass(frozen=True)
class ArrayReference:
    path: str
    shape: tuple
    dtype: np.dtype

    def open(self):
        return np.load(self.path, mmap_mode="r", allow_pickle=False)


@dataclass(frozen=True)
class ValueReference:
    path: str
    shape: tuple | None = None
    dtype: object = None

    def open(self):
        import pickle

        with Path(self.path).open("rb") as stream:
            return pickle.load(stream)


class ProductStorage:
    """Own persisted payloads and their mappings until final publication.

    Read-only arrays are immutable published products and may reuse bytes.
    Writable arrays always receive a fresh snapshot.
    """

    def __init__(self, *, directory=None, budget=None):
        self.temporary = temporary_directory("eyeflow-products-") if directory is None else None
        self.directory = Path(self.temporary.name if self.temporary else directory)
        self.directory.mkdir(parents=True, exist_ok=True)
        self.budget = budget
        self._resident = {}
        self._resident_lock = RLock()
        self._closed = False
        self._maps = []
        self._immutable_arrays = {}
        self._array_lock = RLock()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def path(self, suffix):
        if self._closed:
            raise RuntimeError("Product storage is closed.")
        return self.directory / f"{uuid4().hex}{suffix}"

    def store(self, value, *, category="retained"):
        if isinstance(value, (ArrayReference, ValueReference)):
            return value
        array = np.asarray(value)
        # Object/string payloads remain small metadata; numeric arrays spill.
        numeric = array.dtype.kind in "biufcSU"
        if numeric and array.nbytes >= 64 * 1024:
            return self.store_array(array)
        if self.budget and not self.budget.reserve_memory(array.nbytes, category):
            if numeric:
                return self.store_array(array)
            return self.store_value(array)
        with self._resident_lock:
            self._resident[category] = self._resident.get(category, 0) + array.nbytes
        result = array.copy()
        result.flags.writeable = False
        weakref.finalize(result, self._release_resident, array.nbytes, category)
        return result

    def store_value(self, value):
        import pickle

        path = self.path(".pickle")
        with path.open("wb") as stream:
            pickle.dump(value, stream, protocol=5)
        return ValueReference(str(path), tuple(value.shape), value.dtype)

    def _release_resident(self, size, category):
        with self._resident_lock:
            if not self._closed:
                self._resident[category] -= size
                if self.budget:
                    self.budget.release_memory(size, category)

    def store_array(self, array):
        immutable = _is_immutable_array(array)
        if immutable:
            with self._array_lock:
                cached = self._immutable_arrays.get(id(array))
                if cached is not None and cached[0]() is array:
                    return cached[1]
                reference = self._write_array(array)
                self._remember_array(array, reference)
                return reference
        return self._write_array(array)

    def _remember_array(self, array, reference):
        identity = id(array)

        def forget(expired):
            with self._array_lock:
                cached = self._immutable_arrays.get(identity)
                if cached is not None and cached[0] is expired:
                    del self._immutable_arrays[identity]

        self._immutable_arrays[identity] = weakref.ref(array, forget), reference

    def _write_array(self, array):
        if shutil.disk_usage(self.directory).free < array.nbytes:
            raise OSError("Insufficient local scratch disk space for an array.")
        path = self.path(".npy")
        output = None
        try:
            output = np.lib.format.open_memmap(path, mode="w+", dtype=array.dtype, shape=array.shape)
            if array.ndim and array.shape[0]:
                stride = max(
                    1, 16 * 2**20 // max(1, (array.dtype.itemsize * int(np.prod(array.shape[1:]))))
                )
                for first in range(0, array.shape[0], stride):
                    output[first : first + stride] = array[first : first + stride]
            else:
                output[...] = array
            output.flush()
        except BaseException:
            if output is not None:
                output._mmap.close()
            path.unlink(missing_ok=True)
            raise
        finally:
            if output is not None and not output._mmap.closed:
                output._mmap.close()
        return ArrayReference(str(path), tuple(array.shape), array.dtype)

    def read(self, value):
        if isinstance(value, (ArrayReference, ValueReference)):
            array = value.open()
            if isinstance(array, np.memmap):
                self._maps.append(weakref.ref(array))
                with self._array_lock:
                    self._remember_array(array, value)
            return array
        return value

    def close(self):
        if self._closed:
            return
        for ref in self._maps:
            array = ref()
            if array is not None and not array._mmap.closed:
                array._mmap.close()
        self._maps.clear()
        self._immutable_arrays.clear()
        with self._resident_lock:
            self._closed = True
            if self.budget:
                for category, size in self._resident.items():
                    self.budget.release_memory(size, category)
            self._resident.clear()
        if self.temporary:
            self.temporary.cleanup()
            self.temporary = None


def _is_immutable_array(array):
    """A read-only view of a writable base is still mutable through that base."""
    if not isinstance(array, np.ndarray):
        return False
    current = array
    while isinstance(current, np.ndarray):
        if current.flags.writeable:
            return False
        current = current.base
    if current is None:
        return True
    try:
        return memoryview(current).readonly
    except TypeError:
        return False
