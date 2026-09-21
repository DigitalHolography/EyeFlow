"""Native-pixel widths for optic-disc-centered annuli."""

from __future__ import annotations

import numpy as np

from .segment_geometry import (
    SegmentRingSettings,
    image_half_diagonal,
    section_bounds,
)


def annulus_widths_pixels(
    image_shape: tuple[int, int],
    settings: SegmentRingSettings,
    ring_count: int,
) -> np.ndarray:
    """Return exact radial widths, including a clipped final annulus."""

    radius_scale = image_half_diagonal(*image_shape)
    widths = np.zeros(ring_count, dtype=np.float32)
    for ring_index in range(ring_count):
        inner, outer = section_bounds(settings, ring_index)
        widths[ring_index] = np.float32(max(outer - inner, 0.0) * radius_scale)
    return widths
