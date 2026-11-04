"""Small per-storage counters; no scientific values or scheduling policy."""

from collections import defaultdict
from threading import Lock


class ExecutionMetrics:
    def __init__(self):
        self._values = defaultdict(float)
        self._lock = Lock()

    def add(self, name, value):
        with self._lock:
            self._values[name] += value

    def snapshot(self):
        with self._lock:
            return dict(self._values)
