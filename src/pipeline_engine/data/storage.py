"""Immutable payload storage with explicit mapping lifetime."""

from __future__ import annotations

import shutil
import weakref
import mmap
from dataclasses import dataclass
from pathlib import Path
from threading import RLock
from uuid import uuid4
from math import prod
from threading import get_ident

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


@dataclass(frozen=True)
class H5ArrayReference:
    """Closed, run-owned numerical backing; transfer metadata, never array bytes."""

    path: str
    shape: tuple
    dtype: np.dtype

    def open(self):
        return H5Array(self)


class H5Array:
    """Immutable array-like reads with no persistent HDF5 handles."""

    def __init__(self, reference):
        self._reference = reference

    @property
    def reference(self):
        return self._reference

    @property
    def shape(self):
        return self._reference.shape

    @property
    def dtype(self):
        return self._reference.dtype

    @property
    def ndim(self):
        return len(self.shape)

    @property
    def size(self):
        return prod(self.shape)

    @property
    def nbytes(self):
        return self.size * self.dtype.itemsize

    def __getitem__(self, key):
        import h5py

        selectors = key if isinstance(key, tuple) else (key,)
        if any(item is Ellipsis for item in selectors):
            axis = next(i for i, item in enumerate(selectors) if item is Ellipsis)
            selectors = (*selectors[:axis], *(slice(None),) * (self.ndim - len(selectors) + 1), *selectors[axis + 1:])
        selectors = (*selectors, *(slice(None),) * (self.ndim - len(selectors)))
        sizes = [len(range(*item.indices(size))) for item, size in zip(selectors, self.shape) if isinstance(item, slice)]
        if prod(sizes) * self.dtype.itemsize > 64 * 2**20:
            raise MemoryError("Read bounded slices of a large HDF5-backed product.")
        with h5py.File(self.reference.path, "r") as source:
            return source["data"][key]

    def __array__(self, dtype=None, copy=None):
        if self.nbytes > 64 * 2**20:
            raise MemoryError("Read bounded slices of a large HDF5-backed product, not its full dense array.")
        return np.asarray(self[()], dtype=dtype)

    def copy_to(self, group, name):
        import h5py

        # HDF5 copies allocated chunks, including their filters and fill value.
        # Missing chunks remain unallocated; Python never receives dense data.
        with h5py.File(self.reference.path, "r") as source:
            source.copy("data", group, name=name)
        return group[name]


class _H5ArrayBuilder(H5Array):
    """A physical payload under construction, owned by OutputModel.build_arrays."""

    def __init__(self, storage, shape, options):
        import h5py

        options = dict(options)
        dtype = np.dtype(options.pop("dtype", np.float32))
        reference = H5ArrayReference(str(storage.path(".h5")), tuple(shape), dtype)
        super().__init__(reference)
        self._storage = storage
        self._thread = get_ident()
        self._finished = False
        self._dirty = False
        self._file = None
        try:
            self._file = h5py.File(reference.path, "w", rdcc_nbytes=2**20)
            self._file.create_dataset("data", shape=shape, dtype=dtype, **options)
        except BaseException:
            self.close(success=False)
            raise
        storage.metrics.add("h5_arrays", 1)
        storage.metrics.add("h5_logical_bytes", self.nbytes)

    @property
    def reference(self):
        if not self._finished:
            raise RuntimeError("Unfinished HDF5 array cannot be accepted as a scientific product.")
        return self._reference

    def __getitem__(self, key):
        if self._file is not None:
            return self._file["data"][key]
        return super().__getitem__(key)

    def __setitem__(self, key, value):
        from pipeline_engine.resources import ensure_memory_available

        if self._file is None:
            raise TypeError("Finished scientific array is immutable.")
        if get_ident() != self._thread:
            raise RuntimeError("Only the constructing thread may write an HDF5 array.")
        values = np.asarray(value, dtype=self.dtype)
        ensure_memory_available(min(values.nbytes, 64 * 2**20))
        if shutil.disk_usage(Path(self._reference.path).parent).free < values.nbytes + 16 * 2**20:
            raise OSError("Insufficient local scratch space for the next scientific output chunk.")
        self._file["data"][key] = values
        self._dirty = True
        self._storage.metrics.add("h5_written_bytes", values.nbytes)

    def fill(self, value):
        # Initialization is the HDF5 fill-value declaration, not a page scan.
        if self._file is None or self._dirty or not np.isnan(value):
            raise ValueError("Only initial NaN fill is supported by an incremental array.")

    def close(self, *, success):
        if self._file is not None:
            try:
                self._file.close()
            finally:
                self._file = None
        self._finished = success
        path = Path(self._reference.path)
        if success:
            self._storage.metrics.add("h5_backing_bytes", path.stat().st_size)
        else:
            path.unlink(missing_ok=True)


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
        if isinstance(value, H5Array):
            return value.reference
        if isinstance(value, (ArrayReference, ValueReference, H5ArrayReference)):
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

    def updated_h5_slice(self, reference, key, value):
        """Private chunk-preserving edit; checkpoint rollback retains old backing."""
        import h5py
        from pipeline_engine.resources import ensure_memory_available

        ensure_memory_available(16 * 2**20)
        if shutil.disk_usage(self.directory).free < Path(reference.path).stat().st_size + 16 * 2**20:
            raise OSError("Insufficient local scratch space for a private dataset edit.")
        path = self.path(".h5")
        try:
            with h5py.File(path, "w") as target:
                reference.open().copy_to(target, "data")
                target["data"][key] = value
        except BaseException:
            path.unlink(missing_ok=True)
            raise
        return H5ArrayReference(str(path), reference.shape, reference.dtype)

    def _release_resident(self, size, category):
        with self._resident_lock:
            if not self._closed:
                self._resident[category] -= size
                if self.budget:
                    self.budget.release_memory(size, category)

    def store_array(self, array):
        if isinstance(array, H5Array):
            return array.reference
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
        if isinstance(value, H5ArrayReference):
            self.metrics.add("h5_array_views", 1)
            return value.open()
        if isinstance(value, (ArrayReference, ValueReference)):
            array = value.open()
            self.metrics.add("array_mappings" if isinstance(value, ArrayReference) else "value_loads", 1)
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
