"""Process-local, bounded ownership of authorized immutable replay mappings."""

from collections import OrderedDict
from time import perf_counter
from weakref import WeakValueDictionary


class ReplayMappings:
    """Cache only extents confirmed by a successful coordinator response.

    Eviction drops our reference, never closes a live consumer's ndarray base.
    Weak entries let those consumers' mappings remain reusable without pinning
    an unbounded set of handles. Explicit client teardown ends their lifetime.
    """

    def __init__(self, metric, limit=32):
        self.metric, self.limit = metric, limit
        self._live = WeakValueDictionary()
        self._retained = OrderedDict()
        self.peak_retained = self.peak_live = 0

    def find(self, source, first=None, stop=None):
        started = perf_counter()
        eligible = [(identity, array) for identity, array in list(self._live.items())
                    if identity[0] == source and (identity[3] if first is None else
                                                 identity[1] <= first and stop <= identity[2])]
        if eligible:
            identity, array = max(eligible, key=lambda item: (item[0][3], item[0][2]-item[0][1]))
            self._retain(identity, array)
            self.metric("input_replay_lookup_seconds", perf_counter() - started)
            self.metric("input_replay_cache_hits", 1)
            self.metric("input_replay_zero_rpc_accesses", 1)
            return array, identity[1]
        self.metric("input_replay_lookup_seconds", perf_counter() - started)
        self.metric("input_replay_cache_misses", 1)
        return None

    def authorize(self, source, first, stop, complete, reference):
        identity = source, first, stop, complete
        array = self._live.get(identity)
        if array is None:
            started = perf_counter()
            array = reference.open()
            array.flags.writeable = False
            self.metric("input_replay_mapping_seconds", perf_counter() - started)
            self.metric("input_replay_mappings_opened", 1)
            self._live[identity] = array
        # Larger grants supersede smaller cached mappings. A consumer holding a
        # smaller view keeps its mapping alive through NumPy's base chain.
        for previous in list(self._retained):
            key, lo, hi, _ = previous
            if key == source and first <= lo and hi <= stop:
                self._retained.pop(previous)
        self._retain(identity, array)
        return array

    def _retain(self, identity, array):
        self._retained[identity] = array
        self._retained.move_to_end(identity)
        while len(self._retained) > self.limit:
            self._retained.popitem(last=False)
        self.peak_retained = max(self.peak_retained, len(self._retained))
        self.peak_live = max(self.peak_live, len(self._live))

    @property
    def retained_count(self):
        return len(self._retained)

    @property
    def live_count(self):
        return len(self._live)

    def close(self):
        for array in list(self._live.values()):
            if not array._mmap.closed:
                array._mmap.close()
        self._retained.clear()
        self._live.clear()
