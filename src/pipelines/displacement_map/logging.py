"""Timed app-log messages for displacement stages."""

from contextlib import contextmanager
from time import perf_counter

from utils.logger import Logger


@contextmanager
def logged_stage(label: str, *, frames: int | None = None):
    detail = f" ({frames} frames)" if frames is not None else ""
    Logger.log(f"Starting {label}{detail}...")
    started = perf_counter()
    yield
    Logger.log(f"Completed {label} in {perf_counter() - started:.1f}s.")
