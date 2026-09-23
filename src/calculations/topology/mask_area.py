"""Fractional native-pixel areas for exact optic-disc-centered annuli."""

from __future__ import annotations

from collections.abc import Iterator

import numpy as np

from .geometry import (
    AnnulusGeometry,
    image_half_diagonal,
    section_bounds,
)


def exact_annulus_pixel_coverages(
    image_shape: tuple[int, int],
    optic_disc_center,
    settings: AnnulusGeometry,
    ring_count: int,
) -> Iterator[tuple[int, np.ndarray, float]]:
    """Yield ``(ring, fractional_pixel_coverage, radial_width_pixels)``.

    Vessel masks remain native binary measurements. Only the analytically
    known annulus boundary receives fractional pixel coverage.
    """

    radius_scale = image_half_diagonal(*image_shape)
    center = np.asarray(optic_disc_center, dtype=np.float64).reshape(-1)
    if center.size != 2 or not np.all(np.isfinite(center)):
        raise ValueError("optic-disc center must contain two finite (x, y) values.")
    cx, cy = (float(value) for value in center)
    previous_outer_frac: float | None = None
    previous_outer_coverage: np.ndarray | None = None
    for ring_index in range(ring_count):
        inner_frac, outer_frac = section_bounds(settings, ring_index)
        radial_width = max(outer_frac - inner_frac, 0.0) * radius_scale
        inner_coverage = (
            previous_outer_coverage
            if previous_outer_coverage is not None
            and previous_outer_frac is not None
            and np.isclose(previous_outer_frac, inner_frac, rtol=0.0, atol=1e-12)
            else circle_pixel_coverage(
                image_shape,
                cy=cy,
                cx=cx,
                radius_pixels=inner_frac * radius_scale,
            )
        )
        outer_coverage = circle_pixel_coverage(
            image_shape,
            cy=cy,
            cx=cx,
            radius_pixels=outer_frac * radius_scale,
        )
        yield (
            ring_index,
            np.clip(outer_coverage - inner_coverage, 0.0, 1.0),
            float(radial_width),
        )
        previous_outer_frac = outer_frac
        previous_outer_coverage = outer_coverage


def annulus_widths_pixels(
    image_shape: tuple[int, int],
    settings: AnnulusGeometry,
    ring_count: int,
) -> np.ndarray:
    """Return exact radial widths, including a clipped final annulus."""

    radius_scale = image_half_diagonal(*image_shape)
    widths = np.zeros(ring_count, dtype=np.float32)
    for ring_index in range(ring_count):
        inner, outer = section_bounds(settings, ring_index)
        widths[ring_index] = np.float32(max(outer - inner, 0.0) * radius_scale)
    return widths


def circle_pixel_coverage(
    image_shape: tuple[int, int],
    *,
    cy: float,
    cx: float,
    radius_pixels: float,
) -> np.ndarray:
    """Return the fraction of every native pixel square inside a circle."""

    ny, nx = image_shape
    coverage = np.zeros((ny, nx), dtype=np.float64)
    radius = float(radius_pixels)
    if not np.isfinite(radius) or radius <= 0:
        return coverage

    abs_y = np.abs(np.arange(ny, dtype=np.float64) - float(cy))
    abs_x = np.abs(np.arange(nx, dtype=np.float64) - float(cx))
    min_y = np.maximum(abs_y - 0.5, 0.0)
    min_x = np.maximum(abs_x - 0.5, 0.0)
    max_y = abs_y + 0.5
    max_x = abs_x + 0.5
    radius_sq = radius**2
    nearest_sq = min_y[:, None] ** 2 + min_x[None, :] ** 2
    farthest_sq = max_y[:, None] ** 2 + max_x[None, :] ** 2
    fully_inside = farthest_sq <= radius_sq
    fully_outside = nearest_sq >= radius_sq
    coverage[fully_inside] = 1.0
    for y, x in np.argwhere(~(fully_inside | fully_outside)):
        coverage[y, x] = _pixel_circle_overlap_area(
            x_center=float(x) - float(cx),
            y_center=float(y) - float(cy),
            radius=radius,
        )
    return coverage


def _pixel_circle_overlap_area(
    *,
    x_center: float,
    y_center: float,
    radius: float,
) -> float:
    """Return analytic circle/pixel-square intersection area."""

    x0, x1 = x_center - 0.5, x_center + 0.5
    y0, y1 = y_center - 0.5, y_center + 0.5
    vertices = (
        (x0, y0),
        (x1, y0),
        (x1, y1),
        (x0, y1),
    )
    signed_area = sum(
        _circle_edge_area(start, stop, radius)
        for start, stop in zip(vertices, (*vertices[1:], vertices[0]), strict=True)
    )
    return min(max(abs(float(signed_area)), 0.0), 1.0)


def _circle_edge_area(
    start: tuple[float, float],
    stop: tuple[float, float],
    radius: float,
) -> float:
    """Signed circle-intersection area contributed by one polygon edge."""

    start_x, start_y = start
    direction_x = stop[0] - start_x
    direction_y = stop[1] - start_y
    a = direction_x**2 + direction_y**2
    b = 2.0 * (start_x * direction_x + start_y * direction_y)
    c = start_x**2 + start_y**2 - radius**2
    parameters = [0.0, 1.0]
    discriminant = b**2 - 4.0 * a * c
    if discriminant > 0.0 and a > 0.0:
        root = np.sqrt(discriminant)
        for value in ((-b - root) / (2.0 * a), (-b + root) / (2.0 * a)):
            if 0.0 < value < 1.0:
                parameters.append(float(value))
    parameters.sort()

    area = 0.0
    for left, right in zip(parameters[:-1], parameters[1:], strict=True):
        first_x = start_x + left * direction_x
        first_y = start_y + left * direction_y
        second_x = start_x + right * direction_x
        second_y = start_y + right * direction_y
        middle = 0.5 * (left + right)
        midpoint_x = start_x + middle * direction_x
        midpoint_y = start_y + middle * direction_y
        cross = first_x * second_y - first_y * second_x
        if midpoint_x**2 + midpoint_y**2 <= radius**2:
            area += 0.5 * cross
        else:
            dot = first_x * second_x + first_y * second_y
            area += 0.5 * radius**2 * np.arctan2(cross, dot)
    return area
