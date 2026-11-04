"""Immutable transfers must reuse bytes without caching mutable values."""

import numpy as np

from pipeline_engine.storage import ProductStorage
from pipeline_engine.transfer import dump_result, load_result
from data_manager import DataManager


def test_immutable_arrays_reuse_bytes_across_snapshots_and_preserve_aliases(tmp_path):
    value = np.full(32768, 3.0, dtype=np.float32)
    value.flags.writeable = False
    with ProductStorage(directory=tmp_path) as storage:
        # Reuse is based on established backing ownership, not the producer's
        # reversible WRITEABLE flag.
        data = DataManager(storage)
        value = data.read(data.accept(value))
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


def test_owning_readonly_array_can_reenable_writes_and_must_not_be_cached(tmp_path):
    value = np.ones(32768, dtype=np.float32)
    value.flags.writeable = False
    with ProductStorage(directory=tmp_path) as storage:
        first = dump_result(value, storage)
        value.flags.writeable = True
        value[:] = 9
        value.flags.writeable = False
        second = dump_result(value, storage)
        np.testing.assert_array_equal(load_result(first, storage, {}), 1)
        np.testing.assert_array_equal(load_result(second, storage, {}), 9)
        assert len(list(tmp_path.glob("*.npy"))) == 2


def test_readonly_external_mapping_can_have_a_separate_writable_alias(tmp_path):
    source = tmp_path / "source.npy"
    np.save(source, np.ones(32768, dtype=np.float32))
    with ProductStorage() as storage:
        readonly = np.load(source, mmap_mode="r")
        writable = np.load(source, mmap_mode="r+")
        try:
            first = dump_result(readonly, storage)
            writable[:] = 9
            writable.flush()
            second = dump_result(readonly, storage)
            np.testing.assert_array_equal(load_result(first, storage, {}), 1)
            np.testing.assert_array_equal(load_result(second, storage, {}), 9)
            assert len(list(storage.directory.glob("*.npy"))) == 2
        finally:
            readonly._mmap.close()
            writable._mmap.close()
