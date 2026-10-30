"""Runtime limits for native numeric libraries and parallel workers."""

from __future__ import annotations

import os

def configure_numeric_threads() -> None:
    """Bootstrap defaults; the run resource policy owns active allowances."""
    for name in (
        "OPENBLAS_NUM_THREADS",
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[name] = "1"
    os.environ["ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS"] = "1"


def max_parallel_jobs() -> int:
    """Return the maximum Python worker threads to use around BLAS-backed code."""
    from pipeline_engine.resources import thread_allowance
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
