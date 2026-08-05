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
    optic_disc_width: float | None = None
    optic_disc_height: float | None = None


def ring_masks(
    image_shape: tuple[int, int],
    optic_disc_center,
    settings: SegmentRingSettings,
) -> np.ndarray:
    masks = [
        segment_annulus_mask(
            image_shape,
            optic_disc_center,
            settings,
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
        segment_annulus_mask(
            image_shape,
            optic_disc_center,
            settings,
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
    *,
    optic_disc_width=None,
    optic_disc_height=None,
    optic_disc_boundary_radius_frac: float | None = None,
) -> np.ndarray:
    """Return an annulus in either frame- or optic-disc-relative coordinates.

    The legacy coordinate system expresses radii as fractions of the image
    axes.  When valid optic-disc dimensions and a boundary radius are supplied,
    the coordinate system is rescaled independently along X and Y so that the
    detected optic-disc ellipse lies exactly on ``optic_disc_boundary_radius_frac``.
    """
    ny, nx = image_shape
    cy, cx = optic_disc_center_yx(optic_disc_center, ny, nx)
    disc_geometry = _optic_disc_geometry(
        optic_disc_width,
        optic_disc_height,
        optic_disc_boundary_radius_frac,
    )
    if disc_geometry is None:
        y_grid = np.linspace(0.0, 1.0, ny, dtype=np.float32)[:, None]
        x_grid = np.linspace(0.0, 1.0, nx, dtype=np.float32)[None, :]
        y_distance = y_grid - np.float32(cy / max(ny, 1))
        x_distance = x_grid - np.float32(cx / max(nx, 1))
    else:
        width, height, boundary = disc_geometry
        y_distance = (
            np.arange(ny, dtype=np.float32)[:, None] - np.float32(cy)
        ) * np.float32(2.0 * boundary / height)
        x_distance = (
            np.arange(nx, dtype=np.float32)[None, :] - np.float32(cx)
        ) * np.float32(2.0 * boundary / width)
    radius_sq = x_distance**2 + y_distance**2
    return (radius_sq > inner_radius_frac**2) & (radius_sq <= outer_radius_frac**2)


def segment_annulus_mask(
    image_shape: tuple[int, int],
    optic_disc_center,
    settings: SegmentRingSettings,
    inner_radius_frac: float,
    outer_radius_frac: float,
) -> np.ndarray:
    """Build one segment annulus using the geometry carried by its settings."""
    return annulus_mask(
        image_shape,
        optic_disc_center,
        inner_radius_frac,
        outer_radius_frac,
        optic_disc_width=getattr(settings, "optic_disc_width", None),
        optic_disc_height=getattr(settings, "optic_disc_height", None),
        optic_disc_boundary_radius_frac=settings.inner_radius_frac,
    )


def optic_disc_center_yx(optic_disc_center, ny: int, nx: int) -> tuple[float, float]:
    if optic_disc_center is None:
        return ny / 2.0, nx / 2.0
    center = np.asarray(optic_disc_center, dtype=np.float32).reshape(-1)
    if center.size < 2 or not np.all(np.isfinite(center[:2])):
        return ny / 2.0, nx / 2.0
    return float(center[1]), float(center[0])


def _ring_inner(settings: SegmentRingSettings, ring_index: int) -> float:
    return settings.inner_radius_frac + ring_index * settings.ring_width_frac


def _optic_disc_geometry(width, height, boundary):
    values = []
    for value in (width, height, boundary):
        if value is None:
            return None
        array = np.asarray(value, dtype=np.float32).reshape(-1)
        if array.size == 0 or not np.isfinite(array[0]) or array[0] <= 0:
            return None
        values.append(float(array[0]))
    return tuple(values)
