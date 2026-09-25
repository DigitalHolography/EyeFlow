"""Physical scale helpers shared by retinal topology consumers."""

from __future__ import annotations

import numpy as np

REFERENCE_OPTIC_DISC_DIAMETER_MM = 1.91
DEFAULT_PIXEL_SIZE_MM = 0.0191
SPATIAL_INTERPOLATION_FACTOR = 1.0


def retinal_pixel_size_mm(optic_disc) -> float:
    """Derive native retinal pixel size from available optic-disc geometry."""

    values = np.asarray(
        [optic_disc.width, optic_disc.height],
        dtype=np.float32,
    ).reshape(-1)
    if values.size == 2 and np.all(np.isfinite(values)):
        return float(REFERENCE_OPTIC_DISC_DIAMETER_MM / np.mean(values))
    return DEFAULT_PIXEL_SIZE_MM / (2.0**SPATIAL_INTERPOLATION_FACTOR)


__all__ = ["retinal_pixel_size_mm"]
