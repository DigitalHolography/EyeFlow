"""Batch-wide CPU admission, shared memory reservations and RSS observation."""

from __future__ import annotations

import multiprocessing
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass
from threading import Condition, Event, Thread
from time import monotonic

import psutil

from runtime_limits import (
    configure_native_threads,
    default_thread_allowance,
    set_thread_allowance,
    use_thread_allowance,
)
from utils.logger import Logger

_CURRENT = ContextVar("eyeflow_resources", default=None)
_RESERVATION = ContextVar("eyeflow_compute_reservation", default=None)


class MemoryAccount:
    """Shared byte counters, optionally tracking one worker's owned bytes."""

    categories = ("transient", "retained", "input", "transfer")

    def __init__(
        self,
        limit,
        counters,
        owned=None,
        *,
        observed_rss=None,
        baseline_rss=0,
    ):
        self.memory_limit, self.counters = limit, counters
        self.owned = owned
        self._observed_rss = observed_rss
        self._baseline_rss = int(baseline_rss)

    @property
    def used(self):
        with self.counters.get_lock():
            return dict(zip(self.categories, self.counters[:]))

    def maximum_reservation(self, category="transient"):
        """Return capacity after persistent categories, excluding active peers."""
        index = self.categories.index(category)
        with self.counters.get_lock():
            other_categories = sum(self.counters[:]) - self.counters[index]
        return max(0, self.memory_limit - other_categories)

    def reserve_memory(self, size, category="retained"):
        index = self.categories.index(category)
        if size == 0:
            return True
        with self.counters.get_lock():
            limit = self.memory_limit
            if category in {"retained", "input"}:
                # Retained products must leave room for their next consumer.
                limit -= min(64 * 2**20, self.memory_limit // 4)
            accounted = sum(self.counters[:])
            observed = self._observed_increment()
            if size < 0 or max(accounted, observed) + size > limit:
                return False
            self.counters[index] += size
            if self.owned is not None:
                self.owned[index] += size
            return True

    def _observed_increment(self):
        if self._observed_rss is None:
            return 0
        with self._observed_rss.get_lock():
            rss = int(self._observed_rss.value)
        return max(0, rss - self._baseline_rss)

    def release_memory(self, size, category="retained"):
        if size < 0:
            raise ValueError("Cannot release a negative memory reservation.")
        index = self.categories.index(category)
        with self.counters.get_lock():
            if self.owned is not None:
                size = min(size, self.owned[index])
                self.owned[index] -= size
            self.counters[index] = max(0, self.counters[index] - size)

    def reclaim(self):
        """Reclaim a worker's remaining reservations after its process exits."""
        if self.owned is not None:
            with self.counters.get_lock():
                for index, amount in enumerate(self.owned):
                    self.counters[index] -= amount
                    self.owned[index] = 0


def pipeline_memory_bytes(pipeline, budget) -> int:
    """Convert declared MiB to an attainable transient reservation.

    Persistent retained/input/transfer data reduce the maximum request. Active
    transient reservations are deliberately excluded here: normal admission
    waits for those peers instead of silently shrinking a heavy calculation.
    """
    requested = int(float(pipeline.memory_mb) * 1024**2)
    capacity = budget.maximum_reservation("transient")
    headroom = min(64 * 2**20, capacity // 4)
    if capacity >= requested + headroom:
        return requested
    # A budget exactly equal to the declaration still contains small retained
    # metadata and process overhead. Keep bounded consumer headroom instead of
    # creating an admission request that can never succeed.
    return max(0, capacity - headroom)


@dataclass(frozen=True)
class ResourceRequest:
    threads: int
    memory: int


@dataclass
class _ComputeReservation:
    budget: ResourceBudget
    threads: int
    memory: int


def current_resources():
    return _CURRENT.get()


class ResourceBudget(MemoryAccount):
    """Deterministic admission estimates, shared across an entire batch.

    CPU slots bound concurrent work; memory counters account for managed
    allocations. Neither is an OS hard limit on native code or mapped pages.
    """

    def __init__(self, *, cpu_threads=None, memory_bytes=None):
        self.cpu_limit = max(
            1, int(cpu_threads if cpu_threads is not None else default_thread_allowance()),
        )
        limit = int(memory_bytes if memory_bytes is not None
                    else psutil.virtual_memory().available * 0.5)
        if limit < 1:
            raise ValueError("The RAM budget must be positive.")
        context = multiprocessing.get_context("spawn")
        self._process = psutil.Process()
        baseline_rss = self._process.memory_info().rss
        super().__init__(
            limit, context.Array("q", len(self.categories)),
            observed_rss=context.Value("q", baseline_rss), baseline_rss=baseline_rss,
        )
        self.cpu_used = 0
        self.peak_rss = 0
        self.observed_rss = 0
        self._condition = Condition()
        self._activations = 0
        self._monitor_stop = Event()
        self._monitor = None

    @property
    def available_threads(self):
        with self._condition:
            return self.cpu_limit - self.cpu_used

    def default_request(self, *, threads=None, memory=None):
        """Working estimate used by calculations, processes and renderers."""
        return ResourceRequest(
            min(4, self.cpu_limit) if threads is None else threads,
            min(64 * 2**20, self.memory_limit // 4) if memory is None else memory,
        )

    @contextmanager
    def activate(self):
        with self._condition:
            if self._activations == 0:
                # Native limits are process policy, never temporarily restored
                # by overlapping calculation/finalization threads.
                configure_native_threads(1)
                self.refresh()
                self._monitor_stop.clear()
                self._monitor = Thread(
                    target=self._monitor_loop, name="EyeFlowResourceMonitor", daemon=True,
                )
                self._monitor.start()
            self._activations += 1
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
                self.refresh()
                Logger.log(f"[RESOURCES] Observed peak RAM: {self.peak_rss / 2**30:.2f} GiB.")

    def _monitor_loop(self):
        while not self._monitor_stop.wait(1):
            self.refresh()

    def refresh(self):
        """Observe process-tree memory for admission and peak reporting."""
        try:
            processes = [self._process, *self._process.children(recursive=True)]
        except psutil.Error:
            processes = [self._process]
        rss = 0
        for process in processes:
            try:
                rss += process.memory_info().rss
            except psutil.Error:
                pass
        with self._condition:
            self.observed_rss = rss
            self.peak_rss = max(self.peak_rss, rss)
            with self._observed_rss.get_lock():
                self._observed_rss.value = rss

    def memory_account(self):
        return MemoryAccount(
            self.memory_limit, self.counters,
            multiprocessing.get_context("spawn").Array("q", len(self.categories), lock=False),
            observed_rss=self._observed_rss,
            baseline_rss=self._baseline_rss,
        )

    def _validate_request(self, threads, memory):
        if threads < 1 or threads > self.cpu_limit:
            raise ValueError("CPU reservation must be within the batch CPU limit.")
        if memory < 0 or memory > self.memory_limit:
            raise MemoryError("Calculation reservation exceeds the run RAM budget.")

    def try_acquire(self, threads, memory=0):
        self._validate_request(threads, memory)
        with self._condition:
            if self.cpu_used + threads > self.cpu_limit:
                return False
            if not self.reserve_memory(memory, "transient"):
                return False
            self.cpu_used += threads
            return True

    def acquire(self, threads=1, memory=0, *, flexible=False):
        """Wait for admission and return the admitted number of CPU slots.

        Flexible calculations use the currently free slots up to their request.
        The timeout also observes memory released by cross-process accounts.
        """
        self._validate_request(threads, memory)
        with self._condition:
            while True:
                admitted = min(threads, self.cpu_limit - self.cpu_used) if flexible else threads
                if admitted > 0 and self.try_acquire(admitted, memory):
                    return admitted
                self._condition.wait(0.1)

    def release(self, threads, memory=0):
        with self._condition:
            self.cpu_used -= threads
            self.release_memory(memory, "transient")
            self._condition.notify_all()

    @contextmanager
    def compute(self, threads=None, memory=None):
        if _RESERVATION.get() is not None:
            yield
            return
        request = self.default_request(threads=threads, memory=memory)
        started = monotonic()
        admitted = self.acquire(request.threads, request.memory, flexible=True)
        allocation = _ComputeReservation(self, admitted, request.memory)
        token = _RESERVATION.set(allocation)
        try:
            with use_thread_allowance(admitted):
                delay = monotonic() - started
                if delay > 0.1:
                    Logger.log(f"[RESOURCES] Compute admission delayed {delay:.2f}s.")
                yield
        finally:
            _RESERVATION.reset(token)
            self.release(allocation.threads, allocation.memory)


@contextmanager
def suspend_compute():
    """Release an invocation's slots while waiting for its artifact exporters."""
    reservation = _RESERVATION.get()
    if reservation is None:
        yield
        return
    budget = reservation.budget
    budget.release(reservation.threads, reservation.memory)
    token = _RESERVATION.set(None)
    try:
        yield
    finally:
        reservation.threads = budget.acquire(
            reservation.threads, reservation.memory, flexible=True,
        )
        set_thread_allowance(reservation.threads)
        _RESERVATION.reset(token)
