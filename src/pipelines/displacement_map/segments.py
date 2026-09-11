"""Prepare topology-aligned displacement segments and profile measurements."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

import numpy as np

from calculations.topology import (
    PreparedTopology,
    dilate_segment_masks,
    longitudinal_profiles,
    prepare_segments,
    transverse_profiles,
)


PROFILE_MASK_DILATION_ITERATIONS = 10


@dataclass(frozen=True)
class DisplacementSegmentResult:
    """Displacement measurements indexed by ``(radius, branch, frame, ...)``."""

    maps: np.ndarray | None
    transverse_profiles_unmasked: np.ndarray
    transverse_profiles_masked: np.ndarray
    longitudinal_profiles_unmasked: np.ndarray
    longitudinal_profiles_masked: np.ndarray
    x_sum_profile: np.ndarray
    y_sum_profile: np.ndarray
    radial_movement_amplitude: np.ndarray
    radial_asymmetry_index: np.ndarray


def analyze_displacement_segments(
    displacement_maps: Mapping[str, object],
    topologies: Mapping[str, PreparedTopology],
    *,
    retain_maps: bool,
) -> dict[str, DisplacementSegmentResult]:
    """Prepare and measure one dense displacement map per named vessel."""

    results: dict[str, DisplacementSegmentResult] = {}
    for vessel, displacement_map in displacement_maps.items():
        topology = topologies.get(vessel)
        if topology is None:
            raise KeyError(f"Missing prepared topology for vessel {vessel!r}.")
        ring_count, branch_count = topology.rotation_degrees.shape
        frame_count = int(displacement_map.shape[0])
        canvas_side = int(topology.rotated_masks.shape[-1])
        signal_shape = (ring_count, branch_count, frame_count)
        profile_shape = (*signal_shape, canvas_side)
        map_shape = (*signal_shape, canvas_side, canvas_side, 2)
        profile_masks = dilate_segment_masks(
            topology.rotated_masks,
            iterations=PROFILE_MASK_DILATION_ITERATIONS,
        )
        maps = np.full(map_shape, np.nan, dtype=np.float32) if retain_maps else None
        transverse_unmasked = np.full(profile_shape, np.nan, dtype=np.float32)
        transverse_masked = np.full(profile_shape, np.nan, dtype=np.float32)
        longitudinal_unmasked = np.full(profile_shape, np.nan, dtype=np.float32)
        longitudinal_masked = np.full(profile_shape, np.nan, dtype=np.float32)
        x_sum = np.full(signal_shape, np.nan, dtype=np.float32)
        y_sum = np.full(signal_shape, np.nan, dtype=np.float32)
        radial_amplitude = np.full(signal_shape, np.nan, dtype=np.float32)
        radial_asymmetry = np.full(signal_shape, np.nan, dtype=np.float32)
        for prepared in prepare_segments(
            displacement_map,
            topology,
            spatial_axes=(1, 2),
        ):
            index = (prepared.ring_index, prepared.branch_index)
            vectors = _local_displacement_vectors(
                prepared.rotated,
                float(topology.rotation_degrees[index]),
            )
            magnitude = np.hypot(vectors[..., 0], vectors[..., 1]).astype(
                np.float32,
                copy=False,
            )
            if maps is not None:
                maps[index] = vectors
            transverse_unmasked[index] = transverse_profiles(magnitude)
            transverse_masked[index] = transverse_profiles(
                magnitude,
                profile_masks[index],
            )
            longitudinal_unmasked[index] = longitudinal_profiles(magnitude)
            longitudinal_masked[index] = longitudinal_profiles(
                magnitude,
                profile_masks[index],
            )
            x_sum[index] = _finite_sum(vectors[..., 0], axis=(-2, -1))
            y_sum[index] = _finite_sum(vectors[..., 1], axis=(-2, -1))
            amplitude, asymmetry = _radial_metrics(
                vectors,
                topology.rotated_masks[index],
            )
            radial_amplitude[index] = amplitude
            radial_asymmetry[index] = asymmetry
            del prepared, vectors, magnitude
        results[vessel] = DisplacementSegmentResult(
            maps=maps,
            transverse_profiles_unmasked=transverse_unmasked,
            transverse_profiles_masked=transverse_masked,
            longitudinal_profiles_unmasked=longitudinal_unmasked,
            longitudinal_profiles_masked=longitudinal_masked,
            x_sum_profile=x_sum,
            y_sum_profile=y_sum,
            radial_movement_amplitude=radial_amplitude,
            radial_asymmetry_index=radial_asymmetry,
        )
    return results


def _local_displacement_vectors(
    rotated_components: np.ndarray,
    rotation_degrees: float,
) -> np.ndarray:
    values = np.asarray(rotated_components, dtype=np.float32)
    if values.ndim != 4 or values.shape[1] != 2:
        raise ValueError(
            "prepared displacement segments must have shape "
            "(frame, component, y, x)."
        )
    angle = np.deg2rad(np.float32(rotation_degrees))
    cosine = np.cos(angle)
    sine = np.sin(angle)
    dx = values[:, 0]
    dy = values[:, 1]
    local_x = cosine * dx + sine * dy
    local_y = -sine * dx + cosine * dy
    return np.stack((local_x, local_y), axis=-1).astype(np.float32, copy=False)


def _radial_metrics(
    vectors: np.ndarray,
    vessel_mask: np.ndarray,
    *,
    epsilon: float = 1e-6,
) -> tuple[np.ndarray, np.ndarray]:
    radial_strength = np.abs(vectors[..., 0])
    left_region, right_region = _wall_regions(vessel_mask)
    left = _mean_in_region(radial_strength, left_region)
    right = _mean_in_region(radial_strength, right_region)
    amplitude = np.float32(0.5) * (left + right)
    denominator = left + right + np.float32(epsilon)
    asymmetry = np.divide(
        left - right,
        denominator,
        out=np.full_like(left, np.nan),
        where=np.isfinite(denominator),
    )
    return amplitude, asymmetry


def _wall_regions(vessel_mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mask = np.asarray(vessel_mask, dtype=bool)
    height, width = mask.shape
    left_region = np.zeros_like(mask)
    right_region = np.zeros_like(mask)
    x = np.arange(width)
    for row in np.flatnonzero(np.any(mask, axis=1)):
        mask_x = np.flatnonzero(mask[row])
        left_wall = int(mask_x[0])
        right_wall = int(mask_x[-1])
        centerline = np.float32(0.5 * (left_wall + right_wall))
        band_width = max(int(np.ceil((right_wall - left_wall + 1) / 2.0)), 1)
        left_region[row] = (
            (x >= max(left_wall - band_width, 0))
            & (x <= left_wall)
            & (x < centerline)
        )
        right_region[row] = (
            (x >= right_wall)
            & (x < min(right_wall + band_width + 1, width))
            & (x > centerline)
        )
    return left_region, right_region


def _mean_in_region(values: np.ndarray, region: np.ndarray) -> np.ndarray:
    finite = np.isfinite(values) & region[None, ...]
    count = np.sum(finite, axis=(1, 2), dtype=np.int32)
    total = np.sum(
        values,
        axis=(1, 2),
        dtype=np.float32,
        where=finite,
    )
    return np.divide(
        total,
        count,
        out=np.full(values.shape[0], np.nan, dtype=np.float32),
        where=count > 0,
    )


def _finite_sum(values: np.ndarray, *, axis: tuple[int, ...]) -> np.ndarray:
    finite = np.isfinite(values)
    count = np.sum(finite, axis=axis, dtype=np.int32)
    total = np.sum(values, axis=axis, dtype=np.float32, where=finite)
    total[count == 0] = np.nan
    return total


__all__ = ["DisplacementSegmentResult", "analyze_displacement_segments"]
