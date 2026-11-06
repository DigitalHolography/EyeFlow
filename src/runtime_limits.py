"""Runtime limits for native numeric libraries and parallel workers."""

from __future__ import annotations

import math
import os
import sys
from contextlib import contextmanager
from contextvars import ContextVar
from threading import Lock

_THREADS = ContextVar("eyeflow_compute_threads", default=None)
_NATIVE_LOCK = Lock()
_NATIVE_THREADS = None

def configure_numeric_threads(count: int = 1) -> None:
    """Bootstrap defaults; the run resource policy owns active allowances."""
    count = max(1, int(count))
    for name in (
        "OPENBLAS_NUM_THREADS",
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[name] = str(count)
    os.environ["ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS"] = str(count)


def configure_native_threads(count: int = 1) -> None:
    """Set process-wide native policy at coordinator/worker startup."""
    global _NATIVE_THREADS
    count = max(1, int(count))
    with _NATIVE_LOCK:
        if _NATIVE_THREADS == count:
            return
        for name in (
            "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
            "NUMEXPR_NUM_THREADS", "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS",
        ):
            os.environ[name] = str(count)
        from threadpoolctl import threadpool_limits

        threadpool_limits(limits=count)
        if "cv2" in sys.modules:
            sys.modules["cv2"].setNumThreads(count)
        if "SimpleITK" in sys.modules:
            sys.modules["SimpleITK"].ProcessObject.SetGlobalDefaultNumberOfThreads(count)
        _NATIVE_THREADS = count


def default_thread_allowance() -> int:
    """Use three quarters of the affinity-available logical CPUs."""
    import psutil

    try:
        cpus = len(psutil.Process().cpu_affinity())
    except (AttributeError, psutil.Error):
        cpus = os.cpu_count() or 1
    return max(1, math.floor(0.75 * cpus))


def thread_allowance() -> int:
    return _THREADS.get() or default_thread_allowance()


def set_thread_allowance(count: int):
    """Set a spawned worker's allowance; return a token for scoped callers."""
    return _THREADS.set(max(1, int(count)))


@contextmanager
def use_thread_allowance(count: int):
    token = set_thread_allowance(count)
    try:
        yield
    finally:
        _THREADS.reset(token)


def max_parallel_jobs() -> int:
    """Return the maximum Python worker threads to use around BLAS-backed code."""
    allowance = thread_allowance()
    configured = _positive_int_env("EYEFLOW_MAX_PARALLEL_JOBS")
    return min(allowance, configured) if configured else allowance


def cap_parallel_jobs(requested_jobs: int) -> int:
    """Cap a requested worker count to avoid OpenBLAS thread metadata warnings."""
    if requested_jobs < 1:
        return max_parallel_jobs()
    return max(1, min(requested_jobs, max_parallel_jobs()))


def _positive_int_env(name: str) -> int | None:
    try:
        value = int(os.environ.get(name, ""))
    except ValueError:
        return None
    return value if value > 0 else None
