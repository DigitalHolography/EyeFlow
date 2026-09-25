"""Sparse payloads with logical (ring, branch, ...) axes.

Integer ring/branch indexing does not materialize missing segments. NumPy
conversion explicitly creates a dense array for legacy consumers and exports.
Small geometry/fit metadata remains ordinary NumPy arrays.
"""

from __future__ import annotations

import operator

import numpy as np


class SegmentArray:
    def __init__(self, shape, fill_value=np.nan, dtype=np.float32):
        self.shape = tuple(operator.index(n) for n in shape)
        if len(self.shape) < 2 or any(n < 0 for n in self.shape):
            raise ValueError("segment arrays require nonnegative (ring, branch, ...) axes")
        self.dtype = np.dtype(dtype)
        self.fill_value = self.dtype.type(fill_value)
        self._segments: dict[tuple[int, int], np.ndarray] = {}

    @property
    def ndim(self):
        return len(self.shape)

    @property
    def size(self):
        return int(np.prod(self.shape))

    @property
    def nbytes(self):
        """Actual payload bytes, excluding dictionary overhead."""
        return sum(value.nbytes for value in self._segments.values())

    @property
    def segment_indexes(self):
        return tuple(self._segments)

    def _split_index(self, index):
        if not isinstance(index, tuple) or len(index) < 2:
            return None
        if not all(isinstance(i, (int, np.integer)) for i in index[:2]):
            return None
        key = []
        for i, length in zip(index[:2], self.shape[:2]):
            i = int(i)
            if i < 0:
                i += length
            if not 0 <= i < length:
                raise IndexError("segment index out of bounds")
            key.append(i)
        return tuple(key), index[2:]

    def __getitem__(self, index):
        split = self._split_index(index)
        if split is None:
            return np.asarray(self)[index]
        key, tail = split
        segment = self._segments.get(key)
        if segment is None:
            segment = np.broadcast_to(np.asarray(self.fill_value), self.shape[2:])
        return segment[tail]

    def __setitem__(self, index, value):
        split = self._split_index(index)
        if split is None:
            raise IndexError("assign a single integer (ring, branch) segment at a time")
        key, tail = split
        if key not in self._segments:
            self._segments[key] = np.full(self.shape[2:], self.fill_value, self.dtype)
        self._segments[key][tail] = value

    def __array__(self, dtype=None, copy=None):
        if copy is False:
            raise ValueError("dense conversion of SegmentArray requires a copy")
        output = np.full(self.shape, self.fill_value, dtype=dtype or self.dtype)
        for key, segment in self._segments.items():
            output[key] = segment
        return output

    def copy(self):
        result = SegmentArray(self.shape, self.fill_value, self.dtype)
        result._segments = {key: value.copy() for key, value in self._segments.items()}
        return result
