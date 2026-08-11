"""Annular segment geometry helpers ported from Tools/diskMask.m."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class SegmentRingSettings:
    inner_radius_frac: float
    outer_radius_frac: float
    ring_width_frac: float
    ring_count: int
    segment_length_frac: float | None = None


def ring_masks(
    image_shape: tuple[int, int],
    optic_disc_center,
    settings: SegmentRingSettings,
) -> np.ndarray:
    masks = [
        annulus_mask(
            image_shape,
            optic_disc_center,
            *_ring_bounds(settings, ring_index),
        )
        for ring_index in range(settings.ring_count)
    ]
    return np.asarray(masks, dtype=bool)


def section_masks(
    image_shape: tuple[int, int],
    optic_disc_center,
    settings: SegmentRingSettings,
) -> np.ndarray:
    length = settings.segment_length_frac
    if length is None:
        length = settings.ring_width_frac
    masks = [
        annulus_mask(
            image_shape,
            optic_disc_center,
            *_ring_bounds(settings, ring_index, length),
        )
        for ring_index in range(settings.ring_count)
    ]
    return np.asarray(masks, dtype=bool)


def annulus_mask(
    image_shape: tuple[int, int],
    optic_disc_center,
    inner_radius_frac: float,
    outer_radius_frac: float,
) -> np.ndarray:
    """Return a circular annulus centered on the optic disc.

    Radii are fractions of the distance from the optic-disc center to the
    furthest image corner. This makes an outer radius of ``1.0`` reach the
    image boundary while keeping every annulus circular in pixel coordinates.
    """
    ny, nx = image_shape
    cy, cx = optic_disc_center_yx(optic_disc_center, ny, nx)
    corner_radius = max_corner_radius(ny, nx, cy, cx)
    scale = np.float32(1.0 / max(corner_radius, 1.0))
    y_distance = (
        np.arange(ny, dtype=np.float32)[:, None] - np.float32(cy)
    ) * scale
    x_distance = (
        np.arange(nx, dtype=np.float32)[None, :] - np.float32(cx)
    ) * scale
    radius_sq = x_distance**2 + y_distance**2
    return (radius_sq > inner_radius_frac**2) & (radius_sq <= outer_radius_frac**2)


def optic_disc_center_yx(optic_disc_center, ny: int, nx: int) -> tuple[float, float]:
    if optic_disc_center is None:
        return ny / 2.0, nx / 2.0
    center = np.asarray(optic_disc_center, dtype=np.float32).reshape(-1)
    if center.size < 2 or not np.all(np.isfinite(center[:2])):
        return ny / 2.0, nx / 2.0
    return float(center[1]), float(center[0])


def _ring_inner(settings: SegmentRingSettings, ring_index: int) -> float:
    return settings.inner_radius_frac + ring_index * settings.ring_width_frac


def _ring_bounds(
    settings: SegmentRingSettings,
    ring_index: int,
    length: float | None = None,
) -> tuple[float, float]:
    inner = _ring_inner(settings, ring_index)
    if length is None:
        length = settings.ring_width_frac
    outer = min(settings.outer_radius_frac, inner + length)
    return inner, outer


def max_corner_radius(ny: int, nx: int, cy: float, cx: float) -> float:
    y_dist = max(abs(cy), abs((ny - 1) - cy))
    x_dist = max(abs(cx), abs((nx - 1) - cx))
    return float(np.hypot(y_dist, x_dist))
