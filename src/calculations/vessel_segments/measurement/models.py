"""Typed configuration and results for segment-profile measurement."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from ..sampling.models import SegmentSamplingPlan


@dataclass(frozen=True)
class SegmentMeasurementSettings:
    """Resource and spatial-scale settings for generic profile measurement."""

    pixel_size_mm: float
    submask_size_percentile_kept: float = 0.95
    working_memory_mb: float = 512.0

    def __post_init__(self) -> None:
        if not np.isfinite(self.pixel_size_mm) or self.pixel_size_mm <= 0:
            raise ValueError("pixel_size_mm must be finite and positive.")
        if not np.isfinite(self.working_memory_mb) or self.working_memory_mb <= 0:
            raise ValueError("working_memory_mb must be finite and positive.")
        if not 0 < self.submask_size_percentile_kept <= 1:
            raise ValueError("submask_size_percentile_kept must be in (0, 1].")

    @classmethod
    def from_value(cls, value) -> SegmentMeasurementSettings:
        if isinstance(value, cls):
            return value
        return cls(
            pixel_size_mm=float(value.pixel_size_mm),
            submask_size_percentile_kept=float(
                value.submask_size_percentile_kept
            ),
            working_memory_mb=float(value.working_memory_mb),
        )


@dataclass(frozen=True, slots=True)
class MaskedArrays:
    """Matched unmasked and masked forms of one measured array."""

    unmasked: np.ndarray
    masked: np.ndarray

    def __post_init__(self) -> None:
        unmasked = np.asarray(self.unmasked)
        masked = np.asarray(self.masked)
        if unmasked.shape != masked.shape:
            raise ValueError(
                "masked and unmasked arrays must have the same shape: "
                f"{masked.shape!r} != {unmasked.shape!r}."
            )
        object.__setattr__(self, "unmasked", unmasked)
        object.__setattr__(self, "masked", masked)

    def select(self, *, masked: bool) -> np.ndarray:
        return self.masked if masked else self.unmasked


@dataclass(frozen=True, slots=True)
class CompactSegmentMaps:
    """Maps retained only for valid segments, plus their grid indexes."""

    values: np.ndarray
    indexes: np.ndarray

    def __post_init__(self) -> None:
        values = np.asarray(self.values)
        indexes = np.asarray(self.indexes)
        if indexes.ndim != 2 or indexes.shape[1] != 2:
            raise ValueError("segment map indexes must have shape (segment, 2).")
        if not np.issubdtype(indexes.dtype, np.integer):
            raise ValueError("segment map indexes must contain integers.")
        if values.ndim != 4:
            raise ValueError(
                "compact segment maps must have shape (segment, frame, y, x)."
            )
        if values.shape[0] != indexes.shape[0]:
            raise ValueError("segment map values and indexes must have equal row counts.")
        if indexes.shape[0] != np.unique(indexes, axis=0).shape[0]:
            raise ValueError("segment map indexes must be unique.")
        object.__setattr__(self, "values", values)
        object.__setattr__(self, "indexes", indexes.astype(np.int32, copy=False))

    def row_for(self, ring_index: int, branch_index: int) -> int | None:
        matches = np.flatnonzero(
            (self.indexes[:, 0] == int(ring_index))
            & (self.indexes[:, 1] == int(branch_index))
        )
        return None if matches.size == 0 else int(matches[0])

    def to_dense(self, segment_shape: tuple[int, int]) -> np.ndarray:
        """Expand compact maps to ``(ring, branch, frame, y, x)``."""

        shape = tuple(int(size) for size in segment_shape)
        if len(shape) != 2 or any(size < 0 for size in shape):
            raise ValueError("segment_shape must contain two non-negative sizes.")
        if self.indexes.size and (
            np.any(self.indexes < 0)
            or np.any(self.indexes[:, 0] >= shape[0])
            or np.any(self.indexes[:, 1] >= shape[1])
        ):
            raise ValueError("segment map indexes fall outside segment_shape.")
        dense = np.full((*shape, *self.values.shape[1:]), np.nan, dtype=self.values.dtype)
        if self.indexes.size:
            dense[self.indexes[:, 0], self.indexes[:, 1]] = self.values
        return dense


@dataclass(frozen=True, slots=True, kw_only=True)
class SegmentMeasurements:
    """Measurements for every segment in one spatial sampling plan.

    The existing ``topology`` field holds that plan, not native geometry;
    its ``native`` member provides the source-coordinate SegmentTopology.
    """

    topology: SegmentSamplingPlan
    segment_signal: np.ndarray
    transverse: MaskedArrays
    longitudinal: MaskedArrays
    mean_images: MaskedArrays
    sample_spacing_mm: float
    maps: CompactSegmentMaps | None = None

    def __post_init__(self) -> None:
        segment_shape = self.topology.segment_shape
        profile_side = self.topology.profile_side_pixels
        if self.segment_signal.ndim != 3 or self.segment_signal.shape[:2] != segment_shape:
            raise ValueError(
                "segment_signal must have shape (ring, branch, frame) matching topology."
            )
        expected_prefix = self.segment_signal.shape
        for name, pair in (
            ("transverse", self.transverse),
            ("longitudinal", self.longitudinal),
        ):
            if pair.unmasked.ndim != 4 or pair.unmasked.shape[:3] != expected_prefix:
                raise ValueError(
                    f"{name} profiles must have shape (ring, branch, frame, sample)."
                )
            if pair.unmasked.shape[-1] != profile_side:
                raise ValueError(f"{name} profiles must match the prepared profile size.")
        if self.mean_images.unmasked.ndim != 4:
            raise ValueError("mean images must have shape (ring, branch, y, x).")
        if self.mean_images.unmasked.shape[:2] != segment_shape:
            raise ValueError("mean images must match the topology segment grid.")
        if self.mean_images.unmasked.shape[-2:] != (profile_side, profile_side):
            raise ValueError("mean images must match the prepared profile size.")
        if not np.isfinite(self.sample_spacing_mm) or self.sample_spacing_mm < 0:
            raise ValueError("sample_spacing_mm must be finite and non-negative.")
        if self.maps is not None:
            if self.maps.values.shape[1:] != (
                self.frame_count,
                profile_side,
                profile_side,
            ):
                raise ValueError("retained maps must match the frame and profile sizes.")
            if not np.array_equal(self.maps.indexes, self.topology.valid_indexes()):
                raise ValueError("retained map indexes must match valid topology segments.")

    def profile(
        self,
        direction: Literal["transverse", "longitudinal"],
        *,
        masked: bool,
    ) -> np.ndarray:
        if direction == "transverse":
            return self.transverse.select(masked=masked)
        if direction == "longitudinal":
            return self.longitudinal.select(masked=masked)
        raise ValueError("direction must be 'transverse' or 'longitudinal'.")

    def require_maps(self) -> CompactSegmentMaps:
        if self.maps is None:
            raise RuntimeError(
                "Per-segment maps were not retained; request them during profile analysis."
            )
        return self.maps

    @property
    def frame_count(self) -> int:
        return int(self.segment_signal.shape[-1])

    @property
    def segment_shape(self) -> tuple[int, int]:
        return self.topology.segment_shape


__all__ = [
    "CompactSegmentMaps",
    "MaskedArrays",
    "SegmentMeasurements",
    "SegmentMeasurementSettings",
]
