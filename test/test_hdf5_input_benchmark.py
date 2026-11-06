"""Protect numerical and buffer-ownership contracts of the isolated experiment."""
import importlib.util
from pathlib import Path
import sys

import h5py
import numpy as np
import pytest

_PATH = Path(__file__).resolve().parents[1] / "scripts" / "benchmark_hdf5_input_io.py"
_SPEC = importlib.util.spec_from_file_location("benchmark_hdf5_input_io", _PATH)
benchmark = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = benchmark
_SPEC.loader.exec_module(benchmark)


@pytest.fixture
def source(tmp_path):
    path = tmp_path / "source.h5"
    values = np.random.default_rng(7).normal(size=(9, 4, 16, 17)).astype("f4")
    values[1, 2, 3, 4] = np.nan
    values[6, 1, 0, 0] = np.inf
    with h5py.File(path, "w") as output:
        output.create_dataset("spectrograms", data=values, chunks=(4, 2, 8, 17), compression="gzip")
    return path, values


def test_owned_prefetch_blocks_survive_advance_and_early_close(source):
    path, expected = source
    with h5py.File(path, "r") as file:
        measured = benchmark.MeasuredSource(file["spectrograms"])
        # One frame per block, deliberately smaller than a time chunk.
        iterator = benchmark.prefetch_batches(measured, expected[0].nbytes)
        _, first = next(iterator)
        _, second = next(iterator)
        iterator.close()
        np.testing.assert_array_equal(first, expected[:1])
        np.testing.assert_array_equal(second, expected[1:2])
        assert not np.shares_memory(first, second)


@pytest.mark.parametrize("mode", ["npy", "npy-direct", "ram", "ram-direct", "thread-full", "slab-sync", "slab-thread"])
def test_strategies_match_real_reducer_with_nonfinite_samples(source, tmp_path, mode):
    path, values = source
    config = dict(mode=mode, source=str(path), dataset="spectrograms", buffer_mib=1,
                  max_ram_mib=1, workers=1, timeout=10)
    directory = tmp_path / mode
    directory.mkdir()
    result = benchmark.trial(config, directory)
    gaps = np.flatnonzero(np.indices(values.shape[-2:]).sum(axis=0).ravel() % 3 == 0)
    from pipelines.low_freq_spectra.spatial_reduction import reduce_spatial_block
    expected = reduce_spatial_block(values.copy(), gaps)
    with h5py.File(directory / "reduced.h5") as output:
        for index, value in enumerate(expected):
            np.testing.assert_array_equal(output[str(index)][()], np.maximum(value, 0))
    assert result["logical_bytes"] == values.nbytes
    # The producer still contains the exact original samples after in-place partition.
    with h5py.File(path) as original:
        np.testing.assert_array_equal(original["spectrograms"][()], values)


def test_full_ram_admission_rejects_before_read(source, tmp_path):
    path, _ = source
    with pytest.raises(MemoryError, match="exceeds"):
        benchmark.trial(dict(mode="ram", source=str(path), dataset="spectrograms",
                             buffer_mib=1, max_ram_mib=0, workers=1), tmp_path)


def test_spawn_shared_transfer_and_reader_exception(source, tmp_path, monkeypatch):
    monkeypatch.syspath_prepend(str(_PATH.parent))
    import multiprocessing as mp
    from multiprocessing.shared_memory import SharedMemory
    path, values = source
    config = dict(mode="process-full", source=str(path), dataset="spectrograms",
                  buffer_mib=1, max_ram_mib=1, workers=1, timeout=20)
    result = benchmark.run_trial(config, tmp_path)
    assert result["exact_output_equivalence"]
    assert result["logical_bytes"] == values.nbytes
    memory = SharedMemory(create=True, size=values.nbytes)
    parent, child = mp.get_context("spawn").Pipe(duplex=False)
    process = mp.get_context("spawn").Process(
        target=benchmark.shared_reader,
        args=(str(path), "missing", memory.name, values.shape, values.dtype.str, child, {}))
    try:
        process.start()
        child.close()
        assert parent.poll(20)
        assert parent.recv()[0] == "started"
        assert parent.poll(20)
        status, error = parent.recv()
        assert status == "error" and "missing" in error
        process.join(20)
        assert process.exitcode == 0
    finally:
        if process.is_alive():
            process.terminate()
            process.join()
        parent.close()
        child.close()
        memory.close()
        memory.unlink()


def test_failed_direct_read_removes_incomplete_replay(tmp_path):
    class BrokenSource:
        shape = (2, 3)
        dtype = np.dtype("f4")
        nbytes = 24

        def read_direct(self, destination):
            destination[0] = 42
            raise OSError("incomplete read")

    from pipeline_engine.data.storage import ProductStorage
    with ProductStorage(directory=tmp_path) as storage:
        with pytest.raises(OSError, match="incomplete"):
            benchmark.stage_direct(storage, BrokenSource())
    assert not list(tmp_path.glob("*.npy"))


def test_workflow_comparison_checks_units_and_zero_signs(tmp_path):
    spec = importlib.util.spec_from_file_location(
        "benchmark_hdf5_workflow", _PATH.with_name("benchmark_hdf5_workflow.py"))
    workflow = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(workflow)
    first, second = (tmp_path / name for name in ("first.h5", "second.h5"))
    for path in (first, second):
        with h5py.File(path, "w") as output:
            output.attrs["source_file"] = str(path)
            data = output.create_dataset("velocity", data=[1., np.nan, -0.])
            data.attrs["unit"] = "mm/s"
    assert workflow.compare_hdf5(first, second)["equivalent"]
    with h5py.File(second, "r+") as output:
        output["velocity"].attrs["unit"] = "1"
    with pytest.raises(AssertionError):
        workflow.compare_hdf5(first, second)
    with h5py.File(second, "r+") as output:
        output["velocity"].attrs["unit"] = "mm/s"
        output["velocity"][2] = 0.
    with pytest.raises(AssertionError, match="Numeric bits differ"):
        workflow.compare_hdf5(first, second)
