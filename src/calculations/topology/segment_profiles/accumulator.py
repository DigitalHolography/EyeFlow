"""Streaming accumulation for segment-profile measurements."""

from __future__ import annotations

import numpy as np

from calculations.compute_backend import optional_cupy_backend
from calculations.math import nanmean_float32

from ..profiles import longitudinal_profiles, transverse_profiles
from ..workflow import PreparedTopology
from .models import (
    CompactSegmentMaps,
    MaskedArrays,
    SegmentProfileResult,
)


class SegmentProfileAccumulator:
    """Allocate, populate, and finalize one vessel's profile measurements."""

    def __init__(
        self,
        topology: PreparedTopology,
        *,
        frame_count: int,
        retain_segment_maps: bool,
    ) -> None:
        self.topology = topology
        self.frame_count = int(frame_count)
        segment_shape = topology.segment_shape
        canvas_side = topology.profile_side_pixels
        signal_shape = (*segment_shape, self.frame_count)
        profile_shape = (*signal_shape, canvas_side)
        image_shape = (*segment_shape, canvas_side, canvas_side)

        def filled(shape):
            return np.full(shape, np.nan, dtype=np.float32)

        self.segment_signal = filled(signal_shape)
        self.transverse = MaskedArrays(filled(profile_shape), filled(profile_shape))
        self.longitudinal = MaskedArrays(filled(profile_shape), filled(profile_shape))
        self.mean_images = MaskedArrays(filled(image_shape), filled(image_shape))
        self._mean_sum = np.zeros(image_shape, dtype=np.float64)
        self._mean_count = np.zeros(image_shape, dtype=np.int64)
        self._masked_sum = np.zeros(image_shape, dtype=np.float64)
        self._masked_count = np.zeros(image_shape, dtype=np.int64)

        indexes = topology.valid_indexes().reshape((-1, 2))
        self.maps = (
            CompactSegmentMaps(
                values=filled(
                    (indexes.shape[0], self.frame_count, canvas_side, canvas_side)
                ),
                indexes=indexes,
            )
            if retain_segment_maps
            else None
        )
        self._map_rows = np.full(segment_shape, -1, dtype=np.int32)
        if self.maps is not None and indexes.size:
            self._map_rows[indexes[:, 0], indexes[:, 1]] = np.arange(
                indexes.shape[0], dtype=np.int32
            )

    def add_chunk(self, segment, profile_mask: np.ndarray) -> None:
        index = (segment.ring_index, segment.branch_index)
        rotated = segment.rotated
        rotated_masked = segment.rotated_masked
        if rotated_masked is None:
            raise ValueError(
                "segment profile measurement requires the masked segment companion."
            )
        rotated_mask = self.topology.rotated_masks[index]
        _accumulate_means(
            rotated,
            rotated_masked,
            rotated_mask,
            self._mean_sum[index],
            self._mean_count[index],
            self._masked_sum[index],
            self._masked_count[index],
        )

        transverse_unmasked = transverse_profiles(rotated)
        transverse_masked = transverse_profiles(rotated_masked)
        transverse_masked_for_output = (
            transverse_profiles(rotated, profile_mask)
            if not np.array_equal(profile_mask, rotated_mask)
            else transverse_masked
        )
        frame_slice = segment.frame_slice
        self.segment_signal[index][frame_slice] = _to_numpy(
            _profile_mean(transverse_masked)
        )
        self.transverse.unmasked[index][frame_slice] = _to_numpy(
            transverse_unmasked
        )
        self.transverse.masked[index][frame_slice] = _to_numpy(
            transverse_masked_for_output
        )
        self.longitudinal.unmasked[index][frame_slice] = _to_numpy(
            longitudinal_profiles(rotated)
        )
        self.longitudinal.masked[index][frame_slice] = _to_numpy(
            longitudinal_profiles(rotated_masked)
        )
        if self.maps is not None:
            map_row = int(self._map_rows[index])
            if map_row < 0:
                raise ValueError("Missing compact segment-map row for valid segment.")
            self.maps.values[map_row, frame_slice] = _to_numpy(rotated)

    def finish(self, *, sample_spacing_mm: float) -> SegmentProfileResult:
        np.divide(
            self._mean_sum,
            self._mean_count,
            out=self.mean_images.unmasked,
            where=self._mean_count > 0,
        )
        np.divide(
            self._masked_sum,
            self._masked_count,
            out=self.mean_images.masked,
            where=self._masked_count > 0,
        )
        return SegmentProfileResult(
            topology=self.topology,
            segment_signal=self.segment_signal,
            transverse=self.transverse,
            longitudinal=self.longitudinal,
            mean_images=self.mean_images,
            sample_spacing_mm=sample_spacing_mm,
            maps=self.maps,
        )


def _accumulate_means(
    rotated,
    rotated_masked,
    rotated_mask: np.ndarray,
    total: np.ndarray,
    count: np.ndarray,
    masked_total: np.ndarray,
    masked_count: np.ndarray,
) -> None:
    backend = optional_cupy_backend()
    if backend is not None and isinstance(rotated, backend.cupy.ndarray):
        cupy = backend.cupy
        finite = cupy.isfinite(rotated)
        chunk_total = cupy.asnumpy(
            cupy.sum(
                cupy.where(finite, rotated, cupy.float32(0.0)),
                axis=0,
                dtype=cupy.float64,
            )
        )
        chunk_count = cupy.asnumpy(cupy.sum(finite, axis=0, dtype=cupy.int64))
        masked_finite = cupy.isfinite(rotated_masked) & cupy.asarray(
            rotated_mask, dtype=cupy.bool_
        )[None]
        chunk_masked_total = cupy.asnumpy(
            cupy.sum(
                cupy.where(masked_finite, rotated_masked, cupy.float32(0.0)),
                axis=0,
                dtype=cupy.float64,
            )
        )
        chunk_masked_count = cupy.asnumpy(
            cupy.sum(masked_finite, axis=0, dtype=cupy.int64)
        )
    else:
        values = np.asarray(rotated, dtype=np.float32)
        finite = np.isfinite(values)
        chunk_total = np.sum(values, axis=0, dtype=np.float64, where=finite)
        chunk_count = np.sum(finite, axis=0, dtype=np.int64)
        masked_values = np.asarray(rotated_masked, dtype=np.float32)
        masked_finite = np.isfinite(masked_values) & rotated_mask[None]
        chunk_masked_total = np.sum(
            masked_values,
            axis=0,
            dtype=np.float64,
            where=masked_finite,
        )
        chunk_masked_count = np.sum(masked_finite, axis=0, dtype=np.int64)
    total += chunk_total
    count += chunk_count
    masked_total += chunk_masked_total
    masked_count += chunk_masked_count


def _profile_mean(profiles):
    backend = optional_cupy_backend()
    if backend is not None and isinstance(profiles, backend.cupy.ndarray):
        finite = backend.cupy.isfinite(profiles)
        count = backend.cupy.sum(finite, axis=-1, dtype=backend.cupy.int32)
        total = backend.cupy.sum(
            backend.cupy.where(finite, profiles, backend.cupy.float32(0.0)),
            axis=-1,
            dtype=backend.cupy.float32,
        )
        return backend.cupy.where(
            count > 0,
            total / backend.cupy.maximum(count, 1),
            backend.cupy.float32(np.nan),
        ).astype(backend.cupy.float32, copy=False)
    return nanmean_float32(np.asarray(profiles, dtype=np.float32), axis=-1)


def _to_numpy(values) -> np.ndarray:
    backend = optional_cupy_backend()
    if backend is not None and isinstance(values, backend.cupy.ndarray):
        return backend.cupy.asnumpy(values)
    return np.asarray(values)


__all__ = ["SegmentProfileAccumulator"]
