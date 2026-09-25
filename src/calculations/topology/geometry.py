"""Geometry for optic-disc-centered retinal regions."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class AnnulusGeometry:
    """Radial regions used to identify and sample vessel segments."""

    inner_radius_frac: float
    outer_radius_frac: float
    ring_width_frac: float
    ring_count: int
    segment_length_frac: float | None = None


def ring_masks(
    image_shape: tuple[int, int],
    optic_disc_center,
    settings: AnnulusGeometry,
) -> np.ndarray:
    """Return the configured non-overlapping annuli."""

    return np.asarray(
        [
            annulus_mask(
                image_shape,
                optic_disc_center,
                *_ring_bounds(settings, ring_index),
            )
            for ring_index in range(settings.ring_count)
        ],
        dtype=bool,
    )


def section_masks(
    image_shape: tuple[int, int],
    optic_disc_center,
    settings: AnnulusGeometry,
) -> np.ndarray:
    """Return the annuli in which branch-centered maps are sampled."""

    length = settings.segment_length_frac
    if length is None:
        length = settings.ring_width_frac
    return np.asarray(
        [
            annulus_mask(
                image_shape,
                optic_disc_center,
                *_ring_bounds(settings, ring_index, length),
            )
            for ring_index in range(settings.ring_count)
        ],
        dtype=bool,
    )


def annulus_mask(
    image_shape: tuple[int, int],
    optic_disc_center,
    inner_radius_frac: float,
    outer_radius_frac: float,
) -> np.ndarray:
    """Return a circular annulus centered on the optic disc.

    Radii are fractions of the image half-diagonal, so moving the optic disc
    translates the annulus without changing its size.
    """

    ny, nx = _validated_image_shape(image_shape)
    cx, cy = _validated_center_xy(optic_disc_center)
    scale = np.float32(1.0 / max(image_half_diagonal(ny, nx), 1.0))
    y_distance = (
        np.arange(ny, dtype=np.float32)[:, None] - np.float32(cy)
    ) * scale
    x_distance = (
        np.arange(nx, dtype=np.float32)[None, :] - np.float32(cx)
    ) * scale
    radius_sq = x_distance**2 + y_distance**2
    return (radius_sq > inner_radius_frac**2) & (radius_sq <= outer_radius_frac**2)


def _ring_bounds(
    settings: AnnulusGeometry,
    ring_index: int,
    length: float | None = None,
) -> tuple[float, float]:
    inner = settings.inner_radius_frac + ring_index * settings.ring_width_frac
    if length is None:
        length = settings.ring_width_frac
    return inner, min(settings.outer_radius_frac, inner + length)


def section_bounds(
    settings: AnnulusGeometry,
    ring_index: int,
) -> tuple[float, float]:
    """Return the radial bounds used for one measured vessel segment."""

    return _ring_bounds(settings, ring_index, settings.segment_length_frac)


def image_half_diagonal(ny: int, nx: int) -> float:
    """Return the center-independent image half-diagonal in pixels."""

    return float(np.hypot((ny - 1) / 2.0, (nx - 1) / 2.0))


def _validated_image_shape(image_shape) -> tuple[int, int]:
    shape = tuple(int(size) for size in image_shape)
    if len(shape) != 2 or any(size < 1 for size in shape):
        raise ValueError(f"image_shape must contain two positive sizes, got {shape}.")
    return shape


def _validated_center_xy(center) -> tuple[float, float]:
    values = np.asarray(center, dtype=np.float64).reshape(-1)
    if values.size != 2 or not np.all(np.isfinite(values)):
        raise ValueError("optic-disc center must contain two finite (x, y) values.")
    return float(values[0]), float(values[1])
