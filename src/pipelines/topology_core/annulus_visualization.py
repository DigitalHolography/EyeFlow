"""Shared monochrome annulus outlines for diagnostic images and HDF5 views."""

from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi

from calculations.topology.geometry import AnnulusGeometry, image_half_diagonal, section_bounds

_FOUR_CONNECTED = ndi.generate_binary_structure(2, 1)


def annulus_gray_levels(settings: AnnulusGeometry) -> np.ndarray:
    """Cycle white through gray over the number of overlapping stride phases.

    Size/stride = 2 gives white, gray, white, gray. Adjacent or gapped annuli
    have one phase (white); denser overlap adds evenly spaced grayscale levels.
    Only phases represented by the calculated annuli need a palette entry.
    """
    if settings.ring_count == 0:
        return np.asarray([255], dtype=np.uint8)
    ratio = settings.ring_width_frac / settings.stride_frac
    tolerance = 32 * np.finfo(float).eps * max(ratio, 1.0)
    phases = min(settings.ring_count, max(1, int(np.ceil(ratio - tolerance))))
    return np.rint(np.linspace(255, 128, phases)).astype(np.uint8)


def circle_outline(image_shape, center_xy, radius_pixels: float) -> np.ndarray:
    """Draw the existing one-pixel, native-coordinate circular boundary."""
    center_x, center_y = (float(value) for value in center_xy)
    y, x = np.ogrid[: image_shape[0], : image_shape[1]]
    circle = (y - center_y) ** 2 + (x - center_x) ** 2 <= float(radius_pixels) ** 2
    return circle & ~ndi.binary_erosion(circle, structure=_FOUR_CONNECTED)


def annulus_outline_image(
    image_shape: tuple[int, int],
    center_xy,
    settings: AnnulusGeometry,
    *,
    r0_radius_pixels: float | None = None,
) -> np.ndarray:
    """Return uint8 outlines with zero background and stride-dependent shades.

    Both boundaries of each measured annulus share its shade. Shared/clipped
    boundaries and subpixel raster collisions keep the earliest annulus's
    shade. At 50% overlap, shared full-width boundaries have the same shade.
    This is a display overlay; vessel masks and annulus membership are unchanged.
    """
    scale = max(image_half_diagonal(*image_shape), 1.0)
    if r0_radius_pixels is None:
        r0_radius_pixels = settings.inner_radius_frac * scale
    levels = annulus_gray_levels(settings)
    image = np.zeros(image_shape, dtype=np.uint8)
    image[circle_outline(image_shape, center_xy, r0_radius_pixels)] = levels[0]
    for index in range(settings.ring_count):
        for radius in section_bounds(settings, index):
            outline = circle_outline(image_shape, center_xy, radius * scale)
            image[outline & (image == 0)] = levels[index % levels.size]
    return image
