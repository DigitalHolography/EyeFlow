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
            _ring_inner(settings, ring_index),
            _ring_inner(settings, ring_index) + settings.ring_width_frac,
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
            _ring_inner(settings, ring_index),
            _ring_inner(settings, ring_index) + length,
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
    ny, nx = image_shape
    cy, cx = optic_disc_center_yx(optic_disc_center, ny, nx)
    y_grid = np.linspace(0.0, 1.0, ny, dtype=np.float32)[:, None]
    x_grid = np.linspace(0.0, 1.0, nx, dtype=np.float32)[None, :]
    center_y = np.float32(cy / max(ny, 1))
    center_x = np.float32(cx / max(nx, 1))
    radius_sq = (x_grid - center_x) ** 2 + (y_grid - center_y) ** 2
    return (radius_sq > inner_radius_frac**2) & (radius_sq <= outer_radius_frac**2)


def largest_centered_circle_radius_frac(reference_image, optic_disc_center) -> np.float32:
    shape = np.shape(reference_image)
    if len(shape) < 2:
        raise ValueError("reference_image must have at least two dimensions.")
    ny, nx = int(shape[-2]), int(shape[-1])
    cy, cx = optic_disc_center_yx(optic_disc_center, ny, nx)
    center_y = np.float32(cy / max(ny, 1))
    center_x = np.float32(cx / max(nx, 1))
    corner_y = np.asarray([0.0, 0.0, 1.0, 1.0], dtype=np.float32)
    corner_x = np.asarray([0.0, 1.0, 0.0, 1.0], dtype=np.float32)
    distances = np.sqrt(
        (corner_y - center_y) ** 2 + (corner_x - center_x) ** 2,
        dtype=np.float32,
    )
    return np.float32(max(np.max(distances), np.finfo(np.float32).eps))


def optic_disc_center_yx(optic_disc_center, ny: int, nx: int) -> tuple[float, float]:
    if optic_disc_center is None:
        return ny / 2.0, nx / 2.0
    center = np.asarray(optic_disc_center, dtype=np.float32).reshape(-1)
    if center.size < 2 or not np.all(np.isfinite(center[:2])):
        return ny / 2.0, nx / 2.0
    return float(center[1]), float(center[0])


def _ring_inner(settings: SegmentRingSettings, ring_index: int) -> float:
    return settings.inner_radius_frac + ring_index * settings.ring_width_frac
