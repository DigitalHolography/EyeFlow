"""Local scratch storage and immutable large-array references (never .tmp)."""

from __future__ import annotations

import os
import shutil
import tempfile
import weakref
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from threading import RLock
from uuid import uuid4

import numpy as np


def is_local_path(path):
    path = Path(path).resolve()
    if os.name != "nt":
        # Linux network mounts, when present, must not host scratch storage.
        try:
            import psutil

            mounts = sorted(
                psutil.disk_partitions(all=True), key=lambda p: len(p.mountpoint), reverse=True
            )
            for mount in mounts:
                if path.is_relative_to(Path(mount.mountpoint)):
                    return mount.fstype.lower() not in {"nfs", "nfs4", "cifs", "smbfs", "sshfs"}
        except OSError:
            pass
        return True
    import ctypes

    # DRIVE_FIXED excludes UNC/mapped network drives and removable media.
    return ctypes.windll.kernel32.GetDriveTypeW(str(Path(path.anchor))) == 3


@lru_cache(maxsize=32)
def _is_ssd(anchor):
    """Query seek penalty on Windows; unknown media keep the default location."""
    if os.name != "nt":
        return None
    import ctypes
    from ctypes import wintypes

    kernel = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel.CreateFileW.argtypes = [
        wintypes.LPCWSTR,
        wintypes.DWORD,
        wintypes.DWORD,
        wintypes.LPVOID,
        wintypes.DWORD,
        wintypes.DWORD,
        wintypes.HANDLE,
    ]
    kernel.CreateFileW.restype = wintypes.HANDLE
    kernel.DeviceIoControl.argtypes = [
        wintypes.HANDLE,
        wintypes.DWORD,
        wintypes.LPVOID,
        wintypes.DWORD,
        wintypes.LPVOID,
        wintypes.DWORD,
        ctypes.POINTER(wintypes.DWORD),
        wintypes.LPVOID,
    ]
    kernel.CloseHandle.argtypes = [wintypes.HANDLE]
    handle = kernel.CreateFileW("\\\\.\\" + anchor.rstrip("\\/"), 0, 3, None, 3, 0, None)
    if handle == ctypes.c_void_p(-1).value:
        return None
    try:
        # STORAGE_PROPERTY_QUERY(StorageDeviceSeekPenaltyProperty, PropertyStandardQuery)
        query = (wintypes.DWORD * 3)(7, 0, 0)
        output = ctypes.create_string_buffer(32)
        returned = wintypes.DWORD()
        if (
            kernel.DeviceIoControl(
                handle,
                0x2D1400,
                query,
                ctypes.sizeof(query),
                output,
                len(output),
                ctypes.byref(returned),
                None,
            )
            and returned.value >= 9
        ):
            return output.raw[8] == 0
    finally:
        kernel.CloseHandle(handle)
    return None


@lru_cache(maxsize=8)
def _choose_scratch_root(configured):
    candidates = (
        [Path(configured)]
        if configured
        else [
            Path(os.environ.get("LOCALAPPDATA", tempfile.gettempdir())) / "EyeFlow" / "scratch",
            Path(tempfile.gettempdir()) / "eyeflow-scratch",
        ]
    )
    if not configured and os.name == "nt":
        import psutil

        for partition in psutil.disk_partitions():
            root = Path(partition.mountpoint)
            if is_local_path(root):
                candidates.append(root / "eyeflow-scratch" / Path.home().name)
        candidates.sort(key=lambda path: not (is_local_path(path) and _is_ssd(path.anchor) is True))
    for path in candidates:
        if is_local_path(path):
            try:
                path.mkdir(parents=True, exist_ok=True)
                return path
            except OSError:
                continue
    raise OSError("EyeFlow needs a local scratch directory; set EYEFLOW_SCRATCH_DIR.")


def local_scratch_root():
    return _choose_scratch_root(os.environ.get("EYEFLOW_SCRATCH_DIR"))


def temporary_directory(prefix="eyeflow-", **kwargs):
    worker = os.environ.get("EYEFLOW_WORKER_DIRECTORY")
    if worker:
        return RetainedDirectory(
            tempfile.mkdtemp(prefix=prefix, dir=worker), str(Path(worker).resolve())
        )
    return tempfile.TemporaryDirectory(prefix=prefix, dir=local_scratch_root(), **kwargs)


@dataclass
class RetainedDirectory:
    """Worker scratch is released by its coordinator, never by process exit."""

    name: str
    root: str

    def cleanup(self):
        target, root = Path(self.name).resolve(), Path(self.root).resolve()
        if target == root or not target.is_relative_to(root):
            raise ValueError("Worker scratch cleanup escaped its owning directory.")
        shutil.rmtree(target)

    def __enter__(self):
        return self.name

    def __exit__(self, *exc):
        self.cleanup()


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
    def __init__(self, *, directory=None, budget=None):
        self.temporary = temporary_directory("eyeflow-products-") if directory is None else None
        self.directory = Path(self.temporary.name if self.temporary else directory)
        self.directory.mkdir(parents=True, exist_ok=True)
        self.budget = budget
        self._resident = {}
        self._resident_lock = RLock()
        self._closed = False
        self._maps = []

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def path(self, suffix):
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
        if shutil.disk_usage(self.directory).free < array.nbytes:
            raise OSError("Insufficient local scratch disk space for an array.")
        path = self.path(".npy")
        output = np.lib.format.open_memmap(path, mode="w+", dtype=array.dtype, shape=array.shape)
        try:
            if array.ndim and array.shape[0]:
                stride = max(
                    1, 16 * 2**20 // max(1, (array.dtype.itemsize * int(np.prod(array.shape[1:]))))
                )
                for first in range(0, array.shape[0], stride):
                    output[first : first + stride] = array[first : first + stride]
            else:
                output[...] = array
            output.flush()
        finally:
            output._mmap.close()
        return ArrayReference(str(path), tuple(array.shape), array.dtype)

    def read(self, value):
        if isinstance(value, (ArrayReference, ValueReference)):
            array = value.open()
            if isinstance(array, np.memmap):
                self._maps.append(weakref.ref(array))
            return array
        return value

    def close(self):
        for ref in self._maps:
            array = ref()
            if array is not None and not array._mmap.closed:
                array._mmap.close()
        self._maps.clear()
        with self._resident_lock:
            self._closed = True
            if self.budget:
                for category, size in self._resident.items():
                    self.budget.release_memory(size, category)
            self._resident.clear()
        if self.temporary:
            self.temporary.cleanup()
            self.temporary = None
