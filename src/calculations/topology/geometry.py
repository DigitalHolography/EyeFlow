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
    ring_stride_frac: float | None = None

    @property
    def stride_frac(self) -> float:
        """Radial displacement, independently of measurement support."""
        return self.ring_width_frac if self.ring_stride_frac is None else self.ring_stride_frac

    def metadata(self, image_shape) -> dict:
        """Resolved geometry in both native pixels and half-diagonal fractions."""
        scale = max(image_half_diagonal(*image_shape), 1.0)
        bounds = np.asarray([section_bounds(self, i) for i in range(self.ring_count)],
                            dtype=np.float64).reshape(-1, 2)
        return {
            "annulus_size_pixels": self.ring_width_frac * scale,
            "annulus_stride_pixels": self.stride_frac * scale,
            "annulus_count": self.ring_count,
            "annulus_bounds_pixels": bounds * scale,
            "annulus_bounds_fraction": bounds,
            "annulus_radius_scale_pixels": scale,
            "annulus_radial_coordinate": "fraction_of_image_half_diagonal",
        }


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
    ).reshape((settings.ring_count, *image_shape))


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
    inner = settings.inner_radius_frac + ring_index * settings.stride_frac
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


def annulus_count(extent_pixels: float, size_pixels: float, stride_pixels: float) -> int:
    """Count starts before the outer limit, preserving useful clipped sections.

    Suppress only an extra overlapping tail thinner than min(1 native pixel,
    5% of the size, stride). This keeps deliberately subpixel sampling useful.
    With 50% overlap, even zero-based indices are legacy annuli: retain their
    positive clipped tails even below that threshold. Non-overlapping/gapped
    configurations also retain their final partial section. A scale-aware
    float64 tolerance excludes numerically zero tails at an exact division.
    """
    tolerance = 32 * np.finfo(float).eps * max(extent_pixels, size_pixels, stride_pixels, 1.0)
    if extent_pixels <= tolerance:
        return 0
    count = max(1, int(np.ceil(extent_pixels / stride_pixels)))
    while count and extent_pixels - (count - 1) * stride_pixels <= tolerance:
        count -= 1
    if count and stride_pixels < size_pixels:
        remaining = extent_pixels - (count - 1) * stride_pixels
        legacy_position = (
            abs(2 * stride_pixels - size_pixels) <= tolerance and (count - 1) % 2 == 0
        )
        minimum = min(1.0, 0.05 * size_pixels, stride_pixels)
        if not legacy_position and remaining < minimum - tolerance:
            count -= 1
    return count


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
