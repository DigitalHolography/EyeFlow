"""Bounded sequential HD reads and reductions in the native spectral layout."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from numbers import Integral
import os

import numpy as np

from calculations.math.statistics import partition_median_inplace

DEFAULT_READ_BUFFER_BYTES = 64 * 1024**2
MAX_AUTO_SPATIAL_WORKERS = 8
MIN_PARALLEL_BLOCK_BYTES = 8 * 1024**2


class SpatialReducer:
    """Reuse a bounded thread pool for independent time/frequency medians.

    Workers partition disjoint rows of the owned input buffer. HDF5 access and
    output writes stay on the caller thread. Every task finishes before the
    next source read can reuse the buffer. Automatic selection caps workers at
    eight and keeps small blocks serial; an explicit count allows benchmarking.
    """

    def __init__(self, *, workers=None):
        if workers is not None and (
            isinstance(workers, bool) or not isinstance(workers, Integral) or workers < 1
        ):
            raise ValueError("workers must be a positive integer or None.")
        available = getattr(os, "process_cpu_count", os.cpu_count)() or 1
        self.workers = min(MAX_AUTO_SPATIAL_WORKERS, available) if workers is None else int(workers)
        self._automatic = workers is None
        self._executor = None

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        if self._executor is not None:
            self._executor.shutdown(wait=True)
            self._executor = None

    def __call__(self, block, gap_indices=None, *, include_global=True):
        row_count = block.shape[0]*block.shape[1]
        workers = min(self.workers, row_count)
        if workers == 1 or (self._automatic and block.nbytes < MIN_PARALLEL_BLOCK_BYTES):
            return reduce_spatial_block(block, gap_indices, include_global=include_global)
        if self._executor is None:
            self._executor = ThreadPoolExecutor(
                max_workers=self.workers, thread_name_prefix="spectral-median",
            )
        # Flatten only time/frequency rows; preserve native spatial pixel order.
        rows = block.reshape(row_count, 1, *block.shape[-2:])
        bounds = np.linspace(0, row_count, workers+1, dtype=int)
        futures = [self._executor.submit(
            reduce_spatial_block, rows[start:stop], gap_indices, include_global=include_global,
        ) for start, stop in zip(bounds[:-1], bounds[1:])]
        parts = [future.result() for future in futures]
        results = [None if parts[0][i] is None else np.concatenate(
            [part[i] for part in parts], axis=0,
        ).reshape(block.shape[:2]) for i in range(3)]
        if include_global and gap_indices is not None and len(gap_indices) == np.prod(block.shape[-2:]):
            results[2] = results[1]
        return tuple(results)


def spectral_batches(source, *, read_buffer_bytes=DEFAULT_READ_BUFFER_BYTES):
    """Yield owned, reusable C-contiguous (time, freq, y, x) working blocks.

    HDF5 reads directly into the buffer in sequential time slabs. This reduces
    h5py calls on remote files without loading the entire video into memory.
    Align time slabs to storage chunks when they fit the requested read budget.
    At least one frame is needed, even if it exceeds the budget. Array sources
    are copied into the same owned buffer so partitioning never modifies input.
    Consume each yielded block before requesting the next one.
    """
    if not isinstance(read_buffer_bytes, (int, np.integer)) or read_buffer_bytes <= 0:
        raise ValueError("read_buffer_bytes must be a positive integer.")
    frame_bytes = int(np.prod(source.shape[1:])) * np.dtype(source.dtype).itemsize
    count = min(source.shape[0], max(1, read_buffer_bytes // frame_bytes))
    chunks = getattr(source, "chunks", None)
    if chunks is not None and chunks[0] <= count:
        count = max(chunks[0], count // chunks[0] * chunks[0])
    buffer = np.empty((count, *source.shape[1:]), dtype=source.dtype, order="C")
    direct_read = getattr(source, "read_direct", None)
    for start in range(0, source.shape[0], count):
        stop = min(start + count, source.shape[0])
        block = buffer[:stop-start]
        if callable(direct_read):
            direct_read(buffer, source_sel=np.s_[start:stop], dest_sel=np.s_[:stop-start])
        else:
            np.copyto(block, source[start:stop])
        yield start, block


def reduce_spatial_block(block, gap_indices=None, *, include_global=True):
    """Reduce a mutable owned HD block; only partitioning changes its values.

    The native contiguous layout makes flattening spatial axes a view. Global
    and gap products share the read buffer and validity check. Gap samples are
    gathered before the global in-place partition changes pixel locations.
    """
    full = block.reshape(*block.shape[:2], -1)
    mean = median = gap_median = None
    if include_global:
        valid = np.all(np.isfinite(full), axis=-1)
        # Accumulate safely without allocating a float64 copy of every pixel.
        with np.errstate(invalid="ignore", over="ignore"):
            mean = np.mean(full, axis=-1, dtype=np.float64)
        mean[~valid] = np.nan
    if gap_indices is not None and (not include_global or len(gap_indices) != full.shape[-1]):
        # np.take gives a C-contiguous owned array suitable for in-place median.
        samples = np.take(full, gap_indices, axis=-1)
        if include_global and np.all(valid):
            gap_valid = valid
        else:
            gap_valid = np.all(np.isfinite(samples), axis=-1)
        with np.errstate(invalid="ignore"):
            gap_median = partition_median_inplace(samples)
        gap_median[~gap_valid] = np.nan
    if include_global:
        with np.errstate(invalid="ignore"):
            median = partition_median_inplace(full)
        median[~valid] = np.nan
        if gap_indices is not None and len(gap_indices) == full.shape[-1]:
            gap_median = median
    return mean, median, gap_median
