"""Immutable payload storage with explicit mapping lifetime."""

from __future__ import annotations

import shutil
import weakref
import mmap
from dataclasses import dataclass
from pathlib import Path
from threading import RLock
from uuid import uuid4

import numpy as np
from scratch_location import temporary_directory
from pipeline_engine.data.metrics import ExecutionMetrics


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
        self.metrics = ExecutionMetrics()
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
        # Reuse owned resident snapshots as well as managed file references.
        if not array.flags.writeable:
            with self._array_lock:
                cached = self._immutable_arrays.get(id(array))
                if cached is not None and cached[0]() is array:
                    reference = cached[1]
                    if reference is None:
                        return array
                    if reference.shape == array.shape and reference.dtype == array.dtype:
                        return reference
        # Object/string payloads remain small metadata; numeric arrays spill.
        numeric = array.dtype.kind in "biufcSU"
        if array.dtype.hasobject:
            return self.store_value(array)
        if numeric and (array.nbytes >= 64 * 1024 or isinstance(value, np.memmap)):
            return self.store_array(array)
        if self.budget and not self.budget.reserve_memory(array.nbytes, category):
            if numeric:
                return self.store_array(array)
            return self.store_value(array)
        with self._resident_lock:
            self._resident[category] = self._resident.get(category, 0) + array.nbytes
        result = np.frombuffer(array.tobytes(), dtype=array.dtype).reshape(array.shape)
        self.metrics.add("array_copies", 1)
        self.metrics.add("array_copied_bytes", array.nbytes)
        result.flags.writeable = False
        with self._array_lock:
            self._remember_array(result, None)
        weakref.finalize(result, self._release_resident, array.nbytes, category)
        return result

    def store_value(self, value):
        import cloudpickle

        path = self.path(".pickle")
        try:
            with path.open("wb") as stream:
                cloudpickle.dump(value, stream, protocol=5)
        except BaseException:
            path.unlink(missing_ok=True)
            raise
        self.metrics.add("value_serializations", 1)
        self.metrics.add("value_bytes", path.stat().st_size)
        shape = getattr(value, "shape", None)
        return ValueReference(str(path), tuple(shape) if shape is not None else None,
                              getattr(value, "dtype", None))

    def _release_resident(self, size, category):
        with self._resident_lock:
            if not self._closed:
                self._resident[category] -= size
                if self.budget:
                    self.budget.release_memory(size, category)

    def store_array(self, array):
        if isinstance(array, np.ndarray) and not array.flags.writeable:
            with self._array_lock:
                current = array
                while isinstance(current, np.ndarray):
                    cached = self._immutable_arrays.get(id(current))
                    if (cached is not None and cached[1] is not None and cached[0]() is current
                            and current.shape == array.shape and current.dtype == array.dtype
                            and current.strides == array.strides
                            and current.ctypes.data == array.ctypes.data):
                        return cached[1]
                    current = current.base
        if _is_immutable_array(array):
            with self._array_lock:
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
            self.metrics.add("array_copies", 1)
            self.metrics.add("array_copied_bytes", array.nbytes)
            self.metrics.add("npy_writes", 1)
            self.metrics.add("npy_bytes", path.stat().st_size)
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
        # An owning ndarray can re-enable writes later. Only immutable backing
        # (bytes or a managed read-only map) establishes reusable ownership.
        return False
    if isinstance(current, mmap.mmap):
        # A second writer can mutate an arbitrary file's read-only mapping.
        # Managed maps are reused through the ownership cache above instead.
        return False
    try:
        return memoryview(current).readonly
    except TypeError:
        return False
