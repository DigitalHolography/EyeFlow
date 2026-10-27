"""Run-scoped storage for large arrays that downstream pipelines must replay."""

from __future__ import annotations

import ctypes
import os
import shutil
from math import prod
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from utils.logger import Logger

_FALLBACK_RAM_BUDGET_BYTES = 1024**3
_RAM_BUDGET_FRACTION = 0.25


class RunArrayStorage:
    """Keep modest arrays in RAM and spill larger arrays beside the run output."""

    def __init__(self, output_dir: Path, *, memory_budget_bytes: int | None = None) -> None:
        self.output_dir = Path(output_dir)
        if memory_budget_bytes is None:
            available = _available_memory_bytes()
            memory_budget_bytes = (
                int(available * _RAM_BUDGET_FRACTION)
                if available is not None else _FALLBACK_RAM_BUDGET_BYTES
            )
        self.memory_budget_bytes = int(memory_budget_bytes)
        if self.memory_budget_bytes < 0:
            raise ValueError("memory_budget_bytes must be non-negative.")
        self._resident_bytes = 0
        self._directory: TemporaryDirectory | None = None
        self._maps: list[np.memmap] = []
        self._closed = False

    def __enter__(self) -> RunArrayStorage:
        return self

    def __exit__(self, _type, _value, _traceback) -> None:
        self.close()

    @property
    def scratch_directory(self) -> Path | None:
        return Path(self._directory.name) if self._directory is not None else None

    def empty_float32(self, shape: tuple[int, ...]) -> np.ndarray:
        """Allocate one writable float32 array under the run's shared RAM budget."""
        if self._closed:
            raise RuntimeError("Run scratch storage is closed.")
        size_bytes = prod(shape) * np.dtype(np.float32).itemsize
        if self._resident_bytes + size_bytes <= self.memory_budget_bytes:
            self._resident_bytes += size_bytes
            return np.empty(shape, dtype=np.float32)

        self.output_dir.mkdir(parents=True, exist_ok=True)
        try:
            free_bytes = shutil.disk_usage(self.output_dir).free
        except OSError:
            free_bytes = None
        if free_bytes is not None and free_bytes < size_bytes:
            Logger.log_warning("Run scratch disk space is insufficient; using RAM.")
            self._resident_bytes += size_bytes
            return np.empty(shape, dtype=np.float32)
        if self._directory is None:
            self._directory = TemporaryDirectory(prefix=".ef_scratch_", dir=self.output_dir)
        path = Path(self._directory.name) / f"array_{len(self._maps)}.npy"
        try:
            array = np.lib.format.open_memmap(path, mode="w+", dtype=np.float32, shape=shape)
        except OSError as exc:
            Logger.log_warning(f"Could not create run scratch storage: {exc}; using RAM.")
            self._resident_bytes += size_bytes
            return np.empty(shape, dtype=np.float32)
        self._maps.append(array)
        Logger.log(f"Retaining {size_bytes / 1024**2:.1f} MiB array in run scratch storage.")
        return array

    def close(self) -> None:
        if self._closed:
            return
        self._closed = True
        try:
            for array in self._maps:
                mapping = getattr(array, "_mmap", None)
                if mapping is not None and not mapping.closed:
                    mapping.close()
            self._maps.clear()
        finally:
            if self._directory is not None:
                self._directory.cleanup()
                self._directory = None


def _available_memory_bytes() -> int | None:
    if os.name == "nt":
        class MemoryStatus(ctypes.Structure):
            _fields_ = [
                ("length", ctypes.c_ulong), ("load", ctypes.c_ulong),
                ("total_phys", ctypes.c_ulonglong), ("available_phys", ctypes.c_ulonglong),
                ("total_page", ctypes.c_ulonglong), ("available_page", ctypes.c_ulonglong),
                ("total_virtual", ctypes.c_ulonglong), ("available_virtual", ctypes.c_ulonglong),
                ("available_extended", ctypes.c_ulonglong),
            ]

        status = MemoryStatus()
        status.length = ctypes.sizeof(status)
        return int(status.available_phys) if ctypes.windll.kernel32.GlobalMemoryStatusEx(
            ctypes.byref(status)
        ) else None

    available = None
    try:
        for line in Path("/proc/meminfo").read_text().splitlines():
            if line.startswith("MemAvailable:"):
                available = int(line.split()[1]) * 1024
                break
    except (OSError, ValueError, IndexError):
        pass
    if available is None:
        try:
            available = os.sysconf("SC_AVPHYS_PAGES") * os.sysconf("SC_PAGE_SIZE")
        except (AttributeError, OSError, ValueError):
            return None
        if available <= 0:
            return None

    try:
        limit = Path("/sys/fs/cgroup/memory.max").read_text().strip()
        if limit != "max":
            used = int(Path("/sys/fs/cgroup/memory.current").read_text().strip())
            available = min(available, max(0, int(limit) - used))
    except (OSError, ValueError):
        pass
    return int(available)
