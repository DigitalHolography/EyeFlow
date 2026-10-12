"""Reusable spatial sampling plans and sampled segment chunk contracts."""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from typing import NamedTuple

import numpy as np

from calculations.topology.segments import SegmentTopology


@dataclass(frozen=True)
class SegmentSamplingPlan:
    """Native topology, spatial orientations, and masks ready for sampling.

    Centerline orientations can be refined using a reference signal map.
    The plan contains spatial geometry; temporal execution options belong to
    the streaming functions.
    """

    native: SegmentTopology
    rotation_degrees: np.ndarray
    interpolated_masks: np.ndarray
    rotated_masks: np.ndarray

    @property
    def topology(self) -> SegmentTopology:
        """Compatibility alias for the source-coordinate topology."""

        return self.native

    @property
    def segment_shape(self) -> tuple[int, int]:
        return self.native.segment_shape

    @property
    def profile_side_pixels(self) -> int:
        return int(self.rotated_masks.shape[-1])

    @property
    def valid_segments(self) -> np.ndarray:
        """Segments with both native geometry and a resolved rotation."""

        return self.native.valid_segments & np.isfinite(self.rotation_degrees)

    def valid_indexes(self) -> np.ndarray:
        """Return valid ``(ring, branch)`` indexes in stable row-major order."""

        return np.argwhere(self.valid_segments).astype(np.int32, copy=False)

    def is_aligned_with(self, other: SegmentSamplingPlan) -> bool:
        """Return whether two prepared objects describe the same segment grid."""

        if self is other:
            return True
        if not isinstance(other, SegmentSamplingPlan):
            return False
        if (
            self.segment_shape != other.segment_shape
            or self.profile_side_pixels != other.profile_side_pixels
            or self.native.spatial_shape != other.native.spatial_shape
            or self.native.window_side_pixels != other.native.window_side_pixels
        ):
            return False
        if not np.array_equal(self.native.branch_ids, other.native.branch_ids):
            return False
        if not np.array_equal(self.native.labels, other.native.labels):
            return False
        if not np.array_equal(
            self.native.window_bounds_xyxy,
            other.native.window_bounds_xyxy,
        ):
            return False
        if not np.array_equal(self.rotated_masks, other.rotated_masks):
            return False
        return bool(
            np.allclose(
                self.native.segment_centers_xy,
                other.native.segment_centers_xy,
                equal_nan=True,
            )
            and np.allclose(
                self.rotation_degrees,
                other.rotation_degrees,
                equal_nan=True,
            )
        )


class SampledSegment(NamedTuple):
    """One uniformly resized and upright segment from a retinal map."""

    ring_index: int
    branch_index: int
    rotated: object


SampledSegments = Iterator[SampledSegment]


class SampledSegmentChunk(NamedTuple):
    """One bounded temporal slice of a uniformly transformed segment."""

    ring_index: int
    branch_index: int
    frame_slice: slice
    rotated: object
    rotated_masked: object | None = None


SampledSegmentChunks = Iterator[SampledSegmentChunk]
