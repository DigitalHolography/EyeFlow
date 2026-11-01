"""Immutable transfers must reuse bytes without caching mutable values."""

import numpy as np

from pipeline_engine.storage import ProductStorage
from pipeline_engine.transfer import dump_result, load_result


def test_immutable_arrays_reuse_bytes_across_snapshots_and_preserve_aliases(tmp_path):
    value = np.full(32768, 3.0, dtype=np.float32)
    value.flags.writeable = False
    with ProductStorage(directory=tmp_path) as storage:
        first = dump_result({"a": value, "b": value}, storage)
        second = dump_result({"a": value}, storage)
        assert len(list(tmp_path.glob("*.npy"))) == 1
        loaded = load_result(first, storage, {})
        assert loaded["a"] is loaded["b"]
        assert not loaded["a"].flags.writeable
        np.testing.assert_array_equal(load_result(second, storage, {})["a"], value)
        # A restored immutable reference must not be copied on its next handoff.
        dump_result(loaded, storage)
        assert len(list(tmp_path.glob("*.npy"))) == 1


def test_readonly_view_of_mutable_buffer_refreshes_between_snapshots(tmp_path):
    backing = bytearray(np.full(32768, 1.0, dtype=np.float32).tobytes())
    writable = np.frombuffer(backing, dtype=np.float32)
    value = np.frombuffer(backing, dtype=np.float32)
    value.flags.writeable = False
    with ProductStorage(directory=tmp_path) as storage:
        first = dump_result({"a": value, "b": value}, storage)
        writable[0] = 9
        second = dump_result({"a": value, "b": value}, storage)
        assert len(list(tmp_path.glob("*.npy"))) == 2
        before, after = load_result(first, storage, {}), load_result(second, storage, {})
        assert before["a"] is before["b"]
        assert after["a"] is after["b"]
        assert before["a"][0] == 1
        assert after["a"][0] == 9
