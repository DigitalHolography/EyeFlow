"""Native spatial axes, bounded source access and precision-safe reductions."""

import h5py
import numpy as np
import pytest
from unittest.mock import patch

from calculations.math.statistics import partition_median_inplace
from input_output.schema import EyeFlowOutputPaths
from pipelines.low_freq_spectra.spatial_reduction import SpatialReducer, reduce_spatial_block, spectral_batches
from pipelines.low_freq_spectra.spectrograms import preprocess_spectrograms
from test_gap_spectrogram_preprocessing import FrameSource


@pytest.mark.parametrize("pixels", [1, 2, 7, 8])
@pytest.mark.parametrize("dtype", [np.float32, np.float64, np.int32])
def test_partition_median_matches_double_precision_reference(pixels, dtype):
    rng = np.random.default_rng(53)
    values = (rng.normal(size=(3, 5, pixels))*200).astype(dtype)
    expected = np.median(values.astype(np.float64), axis=-1)
    actual = partition_median_inplace(values.copy())
    np.testing.assert_array_equal(actual, expected)


def test_large_float32_means_and_even_medians_remain_finite():
    limit = np.finfo(np.float32).max
    block = np.array([[[[limit, limit]], [[limit, -limit]]]], dtype=np.float32)
    with np.errstate(over="raise", invalid="raise"):
        mean, median, gap = reduce_spatial_block(block, np.array([0, 1]))
    np.testing.assert_array_equal(mean, [[float(limit), 0]])
    np.testing.assert_array_equal(median, mean)
    assert gap is median  # Full gap coverage reuses the global reduction.


def test_even_float64_midpoint_handles_overflow_and_preserves_subnormals():
    maximum = np.finfo(np.float64).max
    minimum = np.nextafter(0., 1.)
    work = np.array([[maximum, maximum], [minimum, minimum], [minimum, 2*minimum], [-maximum, maximum]])
    with np.errstate(over="raise", invalid="raise"):
        actual = partition_median_inplace(work)
    np.testing.assert_array_equal(actual, [maximum, minimum, 2*minimum, 0])


def test_batch_budget_last_slab_and_inplace_work_do_not_modify_source():
    original = np.arange(5*3*2*4, dtype=np.float32).reshape(5, 3, 2, 4)
    before = original.copy()
    source = FrameSource(original)
    frame_bytes = original[0].nbytes
    shapes = []
    for start, block in spectral_batches(source, read_buffer_bytes=2*frame_bytes):
        shapes.append((start, block.shape[0]))
        flat = block.reshape(*block.shape[:2], -1)
        assert np.shares_memory(flat, block)
        assert block.dtype == np.float32 and block.flags.c_contiguous
        block[:] = -1
    assert shapes == [(0, 2), (2, 2), (4, 1)]
    assert source.requests == [slice(0, 2), slice(2, 4), slice(4, 5)]
    assert source.reads == list(range(5))
    np.testing.assert_array_equal(original, before)


@pytest.mark.parametrize("chunks", [None, (2, 1, 2, 4)])
def test_hdf5_direct_reads_match_reference_and_align_time_chunks(tmp_path, chunks):
    rng = np.random.default_rng(24)
    values = rng.uniform(-1, 5, (7, 3, 2, 4)).astype(np.float32)
    zones = np.array([[True, False, True, False], [False, True, False, True]])
    reference = values.astype(np.float64).reshape(7, 3, -1)
    schema = EyeFlowOutputPaths.active()
    with h5py.File(tmp_path / "hd.h5", "w") as h5:
        source = h5.create_dataset("spectrograms", data=values, chunks=chunks)
        budget = values[0].nbytes*3
        blocks = [(start, block.shape[0]) for start, block in spectral_batches(source, read_buffer_bytes=budget)]
        assert blocks == ([(0, 3), (3, 3), (6, 1)] if chunks is None else [(0, 2), (2, 2), (4, 2), (6, 1)])
        result = preprocess_spectrograms(source, zones, read_buffer_bytes=budget)
        np.testing.assert_array_equal(source[()], values)
    for path, expected in (
        (schema.spectrograms_spatial_average, np.mean(reference, axis=-1)),
        (schema.spectrograms_spatial_median, np.median(reference, axis=-1)),
        (schema.spectrograms_masked, np.median(reference[..., zones.ravel()], axis=-1)),
    ):
        np.testing.assert_array_equal(result[path][0], np.maximum(expected, 0))
        assert result[path][1]["spectrogram_read_buffer_bytes"] == budget


def test_strided_input_readonly_arrays_and_nan_support_are_preserved():
    rng = np.random.default_rng(27)
    values = rng.uniform(1, 5, (4, 3, 3, 6)).astype(np.float32)[..., ::2]
    values[0, 0, 0, 1] = np.nan  # Outside the mask, so only global products fail.
    values[1, 1, 1, 1] = np.inf  # Inside the mask, so gap power also fails.
    before = values.copy()
    values.setflags(write=False)
    zones = np.eye(3, dtype=bool)
    result = preprocess_spectrograms(values, zones)
    schema = EyeFlowOutputPaths.active()
    for path, samples in (
        (schema.spectrograms_spatial_average, values.reshape(4, 3, -1)),
        (schema.spectrograms_spatial_median, values.reshape(4, 3, -1)),
        (schema.spectrograms_masked, values[:, :, zones]),
    ):
        finite = np.isfinite(samples).all(axis=-1)
        expected = np.full((4, 3), np.nan)
        reducer = np.mean if path == schema.spectrograms_spatial_average else np.median
        expected[finite] = reducer(samples[finite].astype(np.float64), axis=-1)
        np.testing.assert_array_equal(result[path][0], expected)
    np.testing.assert_array_equal(values, before)


@pytest.mark.parametrize("coverage", ["none", "all", "partial"])
@pytest.mark.parametrize("include_global", [True, False])
def test_threaded_rows_match_serial_for_odd_row_counts_invalid_values_and_large_power(coverage, include_global):
    rng = np.random.default_rng(725)
    values = rng.uniform(-2, 5, (3, 5, 3, 7)).astype(np.float32)
    values[0, 0, 0, 0] = np.nan
    values[0, 1, 0, 1] = np.inf
    values[0, 2, 0, 2] = -np.inf
    values[2, 4] = np.finfo(np.float32).max
    indices = None if coverage == "none" else np.arange(21)
    if coverage == "partial":
        indices = np.arange(1, 21, 2)
    expected = reduce_spatial_block(values.copy(), indices, include_global=include_global)
    with SpatialReducer(workers=4) as reducer:
        actual = reducer(values.copy(), indices, include_global=include_global)
    for a, b in zip(actual, expected):
        if b is None:
            assert a is None
        else:
            np.testing.assert_array_equal(a, b)


def test_reused_thread_pool_preserves_native_pixels_across_hdf5_batch_reuse(tmp_path):
    rng = np.random.default_rng(73)
    values = rng.uniform(1, 5, (7, 5, 3, 4)).astype(np.float32)
    values[1, 2, 0, 0] = np.inf
    values[4, 3, 0, 1] = np.nan
    mask = np.eye(3, 4, dtype=bool)
    with h5py.File(tmp_path/"source.h5", "w") as h5:
        source = h5.create_dataset("spectrograms", data=values, chunks=(2, 1, 3, 4))
        serial = preprocess_spectrograms(source, mask, workers=1, read_buffer_bytes=values[0].nbytes*2)
        threaded = preprocess_spectrograms(source, mask, workers=3, read_buffer_bytes=values[0].nbytes*2)
        np.testing.assert_array_equal(source[()], values)
    assert set(threaded) == set(serial)
    for path in serial:
        np.testing.assert_array_equal(threaded[path][0], serial[path][0])


def test_thread_failures_are_propagated_after_all_workers_finish():
    from pipelines.low_freq_spectra import spatial_reduction as spatial
    import threading

    entered, finished = threading.Event(), threading.Event()

    def fail_or_finish(block, *args, **kwargs):
        if block[0, 0, 0, 0] == 0:
            if not entered.wait(2):
                raise AssertionError("The other worker did not start")
            raise RuntimeError("partition failed")
        entered.set()
        result = reduce_spatial_block(block, *args, **kwargs)
        finished.set()
        return result

    values = np.arange(4*3*2*2, dtype=np.float32).reshape(4, 3, 2, 2)
    with patch.object(spatial, "reduce_spatial_block", side_effect=fail_or_finish):
        with pytest.raises(RuntimeError, match="partition failed"):
            with SpatialReducer(workers=2) as reducer:
                reducer(values)
    assert finished.is_set()


@pytest.mark.parametrize("workers", [0, -1, True, 1.5])
def test_invalid_worker_counts_are_rejected(workers):
    with pytest.raises(ValueError, match="workers"):
        SpatialReducer(workers=workers)
