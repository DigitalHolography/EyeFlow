"""Application-wide admission accounting and native compute allowances."""

from __future__ import annotations

import math
import multiprocessing
import os
from contextlib import contextmanager
from contextvars import ContextVar
from threading import Condition, Event, Thread
from time import monotonic

import psutil

from utils.logger import Logger

_CURRENT = ContextVar("eyeflow_resources", default=None)
_THREADS = ContextVar("eyeflow_compute_threads", default=None)
_RESERVATION = ContextVar("eyeflow_compute_reservation", default=None)


class MemoryAccount:
    categories = ("transient", "retained", "input", "transfer")

    def __init__(self, limit, counters, owned=None):
        self.memory_limit, self.counters = limit, counters
        self.owned = owned

    @property
    def used(self):
        with self.counters.get_lock():
            return dict(zip(self.categories, self.counters[:]))

    def reserve_memory(self, size, category="retained"):
        with self.counters.get_lock():
            limit = self.memory_limit
            if category in {"retained", "input"}:
                # Retained products must leave room for their next consumer.
                limit -= min(64 * 2**20, self.memory_limit // 4)
            if size < 0 or sum(self.counters[:]) + size > limit:
                return False
            if size > psutil.virtual_memory().available * 0.5:
                return False
            self.counters[self.categories.index(category)] += size
            if self.owned is not None:
                self.owned[self.categories.index(category)] += size
            return True

    def release_memory(self, size, category="retained"):
        with self.counters.get_lock():
            index = self.categories.index(category)
            if self.owned is not None:
                size = min(size, self.owned[index])
                self.owned[index] -= size
            self.counters[index] = max(0, self.counters[index] - size)

    def reclaim(self):
        """Reclaim a dead worker's reservations after its process has joined."""
        if self.owned is not None:
            with self.counters.get_lock():
                for index, amount in enumerate(self.owned):
                    self.counters[index] -= amount
                    self.owned[index] = 0


def usable_cpus():
    try:
        return max(1, len(psutil.Process().cpu_affinity()))
    except (AttributeError, psutil.Error):
        return max(1, os.cpu_count() or 1)


def thread_allowance():
    return _THREADS.get() or max(1, math.floor(0.75 * usable_cpus()))


def current_resources():
    return _CURRENT.get()


def configure_native_threads(count=1):
    """Set defaults before imports, and resize already loaded native pools."""
    count = max(1, int(count))
    for name in (
        "OPENBLAS_NUM_THREADS",
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
        "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS",
    ):
        os.environ[name] = str(count)
    from threadpoolctl import threadpool_limits

    limiter = threadpool_limits(limits=count)
    import sys

    if "cv2" in sys.modules:
        sys.modules["cv2"].setNumThreads(count)
    if "SimpleITK" in sys.modules:
        sys.modules["SimpleITK"].ProcessObject.SetGlobalDefaultNumberOfThreads(count)
    return limiter


class ResourceBudget(MemoryAccount):
    """One budget for retained memory, transient calculations and renderers.

    Reservations are admission estimates, not an OS memory hard limit. Mapped
    arrays live on local disk; their resident pages are observed in tree RSS.
    """

    def __init__(self, *, cpu_threads=None, memory_bytes=None):
        self.cpus = usable_cpus()
        # Sample before starting compute, avoiding subtraction of our own load.
        load = psutil.cpu_percent(interval=0.1) / 100
        self.cpu_limit = max(
            1, int(cpu_threads if cpu_threads is not None else math.floor(0.75 * self.cpus))
        )
        self._dynamic_cpu = max(1, min(self.cpu_limit, math.floor(0.75 * self.cpus * (1 - load))))
        self.memory_limit = int(
            memory_bytes if memory_bytes is not None else psutil.virtual_memory().available * 0.5
        )
        MemoryAccount.__init__(
            self,
            self.memory_limit,
            multiprocessing.get_context("spawn").Array("q", len(self.categories)),
        )
        self.cpu_used = 0
        self.peak_rss = 0
        self._condition = Condition()
        self._sample_time = 0.0
        self._process = psutil.Process()
        self._process.cpu_percent()
        self._processes = {self._process.pid: self._process}
        self._external_load = load
        self._activations = 0
        self._monitor_stop = Event()
        self._monitor = None

    @contextmanager
    def activate(self):
        with self._condition:
            self._activations += 1
            if self._activations == 1:
                self._monitor_stop.clear()
                self._monitor = Thread(
                    target=self._monitor_loop, name="EyeFlowResourceMonitor", daemon=True
                )
                self._monitor.start()
        token = _CURRENT.set(self)
        try:
            yield self
        finally:
            _CURRENT.reset(token)
            with self._condition:
                self._activations -= 1
                last = self._activations == 0
                if last:
                    self._monitor_stop.set()
            if last:
                self._monitor.join()
                Logger.log(f"[RESOURCES] Observed peak RAM: {self.peak_rss / 2**30:.2f} GiB.")

    def _monitor_loop(self):
        while not self._monitor_stop.wait(1):
            self.refresh()

    def refresh(self):
        with self._condition:
            self._refresh()

    def _refresh(self):
        now = monotonic()
        if now - self._sample_time < 1:
            return
        self._sample_time = now
        processes = [self._process, *self._process.children(recursive=True)]
        rss, own_cpu = 0, 0.0
        for process in processes:
            try:
                previous = self._processes.get(process.pid)
                if previous is None:
                    process.cpu_percent()
                    self._processes[process.pid] = process
                else:
                    process = previous
                rss += process.memory_info().rss
                own_cpu += process.cpu_percent() / 100
            except psutil.Error:
                pass
        self.peak_rss = max(self.peak_rss, rss)
        self.observed_rss = rss
        external = max(0.0, psutil.cpu_percent() / 100 - own_cpu / self.cpus)
        self._external_load = 0.8 * self._external_load + 0.2 * external
        self._dynamic_cpu = max(
            1, min(self.cpu_limit, math.floor(0.75 * self.cpus * (1 - self._external_load)))
        )

    def memory_account(self):
        return MemoryAccount(
            self.memory_limit,
            self.counters,
            multiprocessing.get_context("spawn").Array("q", len(self.categories), lock=False),
        )

    def try_acquire(self, threads, memory=0):
        self.refresh()
        with self._condition:
            limit = getattr(self, "_dynamic_cpu", self.cpu_limit)
            if self.cpu_used + threads > limit:
                return False
            if not self.reserve_memory(memory, "transient"):
                return False
            self.cpu_used += threads
            return True

    def release(self, threads, memory=0):
        with self._condition:
            self.cpu_used -= threads
            self.release_memory(memory, "transient")
            self._condition.notify_all()

    @contextmanager
    def compute(self, threads=1, memory=0):
        if _RESERVATION.get() is not None:
            yield
            return
        if memory > self.memory_limit:
            raise MemoryError("Calculation reservation exceeds the run RAM budget.")
        started = monotonic()
        requested = threads
        while True:
            self.refresh()
            free = getattr(self, "_dynamic_cpu", self.cpu_limit) - self.cpu_used
            threads = max(1, min(requested, free))
            if self.try_acquire(threads, memory):
                break
            with self._condition:
                self._condition.wait(0.1)
        delay = monotonic() - started
        token = _THREADS.set(threads)
        allocation = [self, threads, memory]
        reservation = _RESERVATION.set(allocation)
        limiter = None
        try:
            if delay > 0.1:
                Logger.log(f"[RESOURCES] Compute admission delayed {delay:.2f}s.")
            # Python frame pools use the allowance; native pools stay serial.
            limiter = configure_native_threads(1)
            yield
        finally:
            try:
                if limiter is not None:
                    limiter.restore_original_limits()
            finally:
                _RESERVATION.reset(reservation)
                _THREADS.reset(token)
                self.release(allocation[1], memory)


@contextmanager
def suspend_compute():
    """Artifact barriers cannot hold the CPU slots their exporters need."""
    reservation = _RESERVATION.get()
    if reservation is None:
        yield
        return
    budget, threads, memory = reservation
    budget.release(threads, memory)
    token = _RESERVATION.set(None)
    try:
        yield
    finally:
        while True:
            budget.refresh()
            free = getattr(budget, "_dynamic_cpu", budget.cpu_limit) - budget.cpu_used
            threads = max(1, min(threads, free))
            if budget.try_acquire(threads, memory):
                break
            with budget._condition:
                budget._condition.wait(0.1)
        reservation[1] = threads
        _THREADS.set(threads)
        _RESERVATION.reset(token)
