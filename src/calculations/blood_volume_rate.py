"""Pure blood-volume-rate calculations for retinal vessel segments."""

from __future__ import annotations

import numpy as np

from calculations.math import nanmedian
from calculations.topology import annulus_widths_pixels, segment_mask_areas_pixels
from calculations.topology.geometry import AnnulusGeometry

TOTAL_MASKED_EDGES_WINDOW_SIZE = 9
TOTAL_MASKED_EDGES_WINDOW_STRIDE = 1


def circular_lumen_profile_flow(
    profiles,
    left_edge,
    right_edge,
    *,
    profile_pixel_size_mm: float,
) -> np.ndarray:
    """Integrate a transverse velocity profile over a circular lumen.

    Velocity is piecewise linear between profile samples.  At each transverse
    coordinate it is multiplied by the circular chord implied by the two edge
    positions, then integrated analytically and converted to ``mm^3/s``.
    """

    values = np.asarray(profiles, dtype=np.float64)
    left = np.asarray(left_edge, dtype=np.float64)
    right = np.asarray(right_edge, dtype=np.float64)
    if values.ndim != 5:
        raise ValueError("velocity profile must have dimensions (x, t, b, k, r).")
    if left.shape != values.shape[1:] or right.shape != values.shape[1:]:
        raise ValueError(
            "lumen edges must match velocity profile dimensions (t, b, k, r)."
        )
    if not np.isfinite(profile_pixel_size_mm) or profile_pixel_size_mm <= 0:
        raise ValueError("profile_pixel_size_mm must be finite and positive.")

    center = 0.5 * (left + right)
    radius = 0.5 * (right - left)
    valid_edges = np.isfinite(left) & np.isfinite(right) & (radius > 0.0)
    safe_radius = np.where(valid_edges, radius, 1.0)
    integral = np.zeros(left.shape, dtype=np.float64)
    has_finite_interval = np.zeros(left.shape, dtype=bool)

    for index in range(max(values.shape[0] - 1, 0)):
        start = np.maximum(left, float(index))
        stop = np.minimum(right, float(index + 1))
        start_value = values[index]
        stop_value = values[index + 1]
        active = (
            valid_edges
            & (stop > start)
            & np.isfinite(start_value)
            & np.isfinite(stop_value)
        )
        slope = stop_value - start_value
        velocity_at_center = start_value + slope * (center - float(index))
        start_u = start - center
        stop_u = stop - center
        contribution = _linear_velocity_chord_antiderivative(
            stop_u,
            safe_radius,
            velocity_at_center,
            slope,
        ) - _linear_velocity_chord_antiderivative(
            start_u,
            safe_radius,
            velocity_at_center,
            slope,
        )
        integral += np.where(active, contribution, 0.0)
        has_finite_interval |= active

    integral *= float(profile_pixel_size_mm) ** 2
    integral[~has_finite_interval] = np.nan
    return integral.astype(np.float32)


def _linear_velocity_chord_antiderivative(
    coordinate,
    radius,
    center_velocity,
    slope,
) -> np.ndarray:
    clipped = np.clip(np.asarray(coordinate, dtype=np.float64), -radius, radius)
    radial_square = np.maximum(radius**2 - clipped**2, 0.0)
    root = np.sqrt(radial_square)
    normalized = np.clip(clipped / radius, -1.0, 1.0)
    constant_term = clipped * root + radius**2 * np.arcsin(normalized)
    linear_term = -(2.0 / 3.0) * radial_square * root
    return center_velocity * constant_term + slope * linear_term


def mask_derived_lumen_geometry(
    topologies,
    *,
    pixel_size_mm: float,
) -> tuple[tuple[np.ndarray, ...], tuple[np.ndarray, ...], np.ndarray]:
    """Return equivalent diameters, native mask areas, and annulus widths."""

    if not np.isfinite(pixel_size_mm) or pixel_size_mm <= 0:
        raise ValueError("pixel_size_mm must be finite and positive.")
    topologies = tuple(topologies)
    if not topologies:
        raise ValueError("at least one segment topology is required.")

    ring_settings: AnnulusGeometry | None = None
    spatial_shape: tuple[int, int] | None = None
    mask_areas: list[np.ndarray] = []
    for value in topologies:
        topology = getattr(value, "topology", value)
        nested = getattr(topology, "topology", None)
        topology = nested if nested is not None else topology
        settings = topology.ring_settings
        if not isinstance(settings, AnnulusGeometry):
            raise ValueError("segment topology must retain its ring settings.")
        if ring_settings is None:
            ring_settings = settings
        elif settings != ring_settings:
            raise ValueError("artery and vein segment ring settings must match.")
        labels = np.asarray(topology.labels, dtype=np.int32)
        if spatial_shape is None:
            spatial_shape = labels.shape
        elif labels.shape != spatial_shape:
            raise ValueError("artery and vein segment image shapes must match.")
        mask_areas.append(segment_mask_areas_pixels(topology))

    if ring_settings is None or spatial_shape is None:
        raise RuntimeError("segment geometry preparation produced no topology.")
    ring_count = int(mask_areas[0].shape[1])
    if any(area.shape[1] != ring_count for area in mask_areas):
        raise ValueError("artery and vein segment ring counts must match.")
    radial_widths = annulus_widths_pixels(spatial_shape, ring_settings, ring_count)

    diameters: list[np.ndarray] = []
    for area in mask_areas:
        diameter = np.full(area.shape, np.nan, dtype=np.float32)
        valid_widths = radial_widths > 0
        if np.any(valid_widths):
            calculated = (
                area[:, valid_widths].astype(np.float32)
                * np.float32(pixel_size_mm)
                / radial_widths[None, valid_widths]
            )
            diameter[:, valid_widths] = np.where(
                area[:, valid_widths] > 0,
                calculated,
                np.float32(np.nan),
            )
        diameters.append(diameter)
    return tuple(diameters), tuple(mask_areas), radial_widths


def masked_edges_flow(velocity, diameter_mm) -> np.ndarray:
    """Multiply safe per-beat velocity by equivalent circular lumen area."""

    velocity_tbkr = np.asarray(velocity, dtype=np.float32)
    diameter = np.asarray(diameter_mm, dtype=np.float32)
    if velocity_tbkr.ndim != 4:
        raise ValueError(
            "safe per-beat segment velocity must have dimensions "
            "(time, beat, branch, radius)."
        )
    if velocity_tbkr.shape[2:] != diameter.shape:
        raise ValueError(
            "safe per-beat velocity branch/radius dimensions must match geometry."
        )
    area_mm2 = np.float32(np.pi / 4.0) * diameter**2
    return (velocity_tbkr * area_mm2[None, None, :, :]).astype(
        np.float32,
        copy=False,
    )


def total_masked_edges_flow(masked_edges) -> np.ndarray:
    """Apply the established temporal, branch, and radius reductions."""

    values = np.asarray(masked_edges, dtype=np.float32)
    if values.ndim != 4:
        raise ValueError(
            "masked-edge blood-volume rate must have dimensions "
            "(time, beat, branch, radius)."
        )
    if values.shape[0] > 0:
        half_window = TOTAL_MASKED_EDGES_WINDOW_SIZE // 2
        periodic = np.pad(
            values,
            ((half_window, half_window), (0, 0), (0, 0), (0, 0)),
            mode="wrap",
        )
        windows = np.lib.stride_tricks.sliding_window_view(
            periodic,
            window_shape=TOTAL_MASKED_EDGES_WINDOW_SIZE,
            axis=0,
        )[::TOTAL_MASKED_EDGES_WINDOW_STRIDE]
        values = np.mean(windows, axis=-1, dtype=np.float32)

    finite = np.isfinite(values)
    rate_tbr = np.sum(
        np.where(finite, values, np.float32(0.0)),
        axis=2,
        dtype=np.float32,
    )
    rate_tbr[~np.any(finite, axis=2)] = np.nan
    return nanmedian(rate_tbr, axis=2).astype(np.float32, copy=False)


__all__ = [
    "TOTAL_MASKED_EDGES_WINDOW_SIZE",
    "TOTAL_MASKED_EDGES_WINDOW_STRIDE",
    "circular_lumen_profile_flow",
    "mask_derived_lumen_geometry",
    "masked_edges_flow",
    "total_masked_edges_flow",
]
