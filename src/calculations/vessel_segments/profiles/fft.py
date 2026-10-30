"""FFT calculations for vessel-aligned segment profiles."""

from __future__ import annotations

from time import perf_counter

import numpy as np

from calculations.math.cycle_boundaries import (
    normalize_cycle_boundaries,
)
from calculations.compute_backend import optional_cupy_backend
from calculations.math import interpft_axis0, nanmean_float32, next_power_of_two
from utils.logger import Logger

DEFAULT_PROFILE_MASK_DILATION_PIXELS = 10
_FFT_PROFILE_X_BATCH = 32


class SegmentFftAccumulator:
    """Accumulate FFT profiles while each rotated segment is resident."""

    def __init__(
        self,
        *,
        frame_count: int,
        ring_count: int,
        branch_count: int,
        canvas_side: int,
        cycle_boundary_indexes,
        index_base: int,
    ) -> None:
        self.boundaries = normalize_cycle_boundaries(
            cycle_boundary_indexes,
            frame_count,
            index_base=index_base,
        )
        self.time_count = next_power_of_two(
            int(np.max(np.diff(self.boundaries)))
        )
        beat_count = self.boundaries.size - 1
        output_shape = (
            canvas_side,
            self.time_count,
            beat_count,
            branch_count,
            ring_count,
        )
        self.unmasked = np.full(output_shape, np.nan, dtype=np.float32)
        self.masked = np.full(output_shape, np.nan, dtype=np.float32)
        self.elapsed_seconds = 0.0
        self._pending_key = None
        self._pending_parts = []

    def observe(
        self,
        ring_index: int,
        branch_index: int,
        frame_slice,
        rotated_chunk=None,
        profile_mask: np.ndarray | None = None,
    ) -> None:
        """Consume one ordered chunk while retaining at most one beat."""

        started = perf_counter()
        if profile_mask is None:
            profile_mask = rotated_chunk
            rotated_chunk = frame_slice
            frame_slice = slice(0, int(rotated_chunk.shape[0]))
        chunk_start = int(frame_slice.start or 0)
        chunk_stop = int(frame_slice.stop)
        mask = np.asarray(profile_mask, dtype=bool)
        if rotated_chunk.ndim != 3 or mask.shape != rotated_chunk.shape[1:]:
            raise ValueError(
                "FFT profile input must contain a (frame, y, x) chunk and "
                "matching (y, x) mask."
            )
        for beat_index in range(self.boundaries.size - 1):
            beat_start = int(self.boundaries[beat_index])
            beat_stop = int(self.boundaries[beat_index + 1]) + 1
            overlap_start = max(chunk_start, beat_start)
            overlap_stop = min(chunk_stop, beat_stop)
            if overlap_start >= overlap_stop:
                continue
            part = rotated_chunk[
                overlap_start - chunk_start : overlap_stop - chunk_start
            ].copy()
            key = (ring_index, branch_index, beat_index)
            if self._pending_key not in (None, key):
                raise RuntimeError("FFT chunks arrived out of segment/frame order.")
            self._pending_key = key
            self._pending_parts.append(part)
            buffered = sum(int(value.shape[0]) for value in self._pending_parts)
            expected = beat_stop - beat_start
            if buffered == expected:
                backend = optional_cupy_backend()
                if (
                    backend is not None
                    and isinstance(self._pending_parts[0], backend.cupy.ndarray)
                ):
                    beat = backend.cupy.concatenate(self._pending_parts, axis=0)
                else:
                    beat = np.concatenate(self._pending_parts, axis=0)
                self._write_beat(
                    ring_index,
                    branch_index,
                    beat_index,
                    beat,
                    mask,
                )
                self._pending_key = None
                self._pending_parts = []
            elif buffered > expected:
                raise RuntimeError("FFT chunk buffering exceeded the current beat.")
        self.elapsed_seconds += perf_counter() - started

    def _write_beat(
        self,
        ring_index: int,
        branch_index: int,
        beat_index: int,
        beat,
        mask: np.ndarray,
    ) -> None:
        backend = optional_cupy_backend()
        if backend is not None and isinstance(beat, backend.cupy.ndarray):
            try:
                self._write_beat_gpu(
                    ring_index,
                    branch_index,
                    beat_index,
                    beat,
                    mask,
                    backend.cupy,
                )
                backend.cupy.cuda.get_current_stream().synchronize()
                return
            except Exception as exc:
                Logger.log_debug(
                    "CuPy segment-profile FFT failed; using CPU fallback: "
                    f"{type(exc).__name__}: {exc}"
                )
                beat = backend.cupy.asnumpy(beat)
        self._write_beat_cpu(
            ring_index,
            branch_index,
            beat_index,
            np.asarray(beat, dtype=np.float32),
            mask,
        )

    def _write_beat_cpu(
        self,
        ring_index: int,
        branch_index: int,
        beat_index: int,
        beat: np.ndarray,
        mask: np.ndarray,
    ) -> None:
        for x_start in range(0, beat.shape[2], _FFT_PROFILE_X_BATCH):
            x_stop = min(x_start + _FFT_PROFILE_X_BATCH, beat.shape[2])
            interpolated = interpft_axis0(
                beat[:, :, x_start:x_stop],
                self.time_count + 1,
            )[:-1]
            magnitude = np.abs(
                np.fft.fft(interpolated, axis=0)
            ).astype(np.float32, copy=False)
            output_slice = (
                slice(x_start, x_stop),
                slice(None),
                beat_index,
                branch_index,
                ring_index,
            )
            self.unmasked[output_slice] = nanmean_float32(
                magnitude,
                axis=1,
            ).T
            self.masked[output_slice] = nanmean_float32(
                np.where(
                    mask[None, :, x_start:x_stop],
                    magnitude,
                    np.float32(np.nan),
                ),
                axis=1,
            ).T

    def _write_beat_gpu(
        self,
        ring_index: int,
        branch_index: int,
        beat_index: int,
        beat,
        mask: np.ndarray,
        cupy,
    ) -> None:
        pixel_count = int(beat.shape[1] * beat.shape[2])
        flattened = beat.reshape((beat.shape[0], pixel_count))
        active_pixels = cupy.any(cupy.isfinite(flattened), axis=0)
        if int(cupy.count_nonzero(active_pixels)) == 0:
            return
        interpolated = interpft_axis0(
            flattened[:, active_pixels],
            self.time_count + 1,
            array_module=cupy,
        )[:-1]
        magnitude_active = cupy.abs(
            cupy.fft.fft(interpolated, axis=0)
        ).astype(cupy.float32, copy=False)
        magnitude = cupy.full(
            (self.time_count, pixel_count),
            cupy.nan,
            dtype=cupy.float32,
        )
        magnitude[:, active_pixels] = magnitude_active
        magnitude = magnitude.reshape(
            (self.time_count, beat.shape[1], beat.shape[2])
        )
        gpu_mask = cupy.asarray(mask, dtype=cupy.bool_)
        unmasked = _gpu_nanmean_axis1(magnitude, cupy)
        masked = _gpu_nanmean_axis1(magnitude, cupy, mask=gpu_mask)
        output_slice = (
            slice(None),
            slice(None),
            beat_index,
            branch_index,
            ring_index,
        )
        self.unmasked[output_slice] = cupy.asnumpy(unmasked).T
        self.masked[output_slice] = cupy.asnumpy(masked).T


def _gpu_nanmean_axis1(values, cupy, *, mask=None):
    """Return a NaN-aware axis-one mean using a CuPy-like array module."""

    finite = cupy.isfinite(values)
    if mask is not None:
        finite &= mask[None, ...]
    counts = cupy.sum(finite, axis=1)
    totals = cupy.sum(
        cupy.where(finite, values, cupy.float32(0.0)),
        axis=1,
        dtype=cupy.float32,
    )
    safe_counts = counts.copy()
    safe_counts[safe_counts == 0] = 1
    result = cupy.divide(totals, safe_counts)
    result[counts == 0] = cupy.nan
    return result.astype(cupy.float32, copy=False)


__all__ = [
    "DEFAULT_PROFILE_MASK_DILATION_PIXELS",
    "SegmentFftAccumulator",
]
