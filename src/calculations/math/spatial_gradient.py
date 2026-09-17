"""Temporal median and framewise Sobel filters for moment0 stacks."""

from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi

TEMPORAL_MEDIAN_WINDOW = 9


def temporal_median_window(stack, *, window=TEMPORAL_MEDIAN_WINDOW, array_module=np):
    """Filter each pixel in time, replicating the first and last frames."""
    if window < 1 or window % 2 == 0:
        raise ValueError("temporal median window must be a positive odd integer.")
    xp = array_module
    values = xp.asarray(stack, dtype=xp.float32)
    filtered = xp.empty_like(values)
    offsets = xp.arange(window) - window // 2
    for index in range(len(values)):
        indexes = xp.clip(index + offsets, 0, len(values) - 1)
        filtered[index] = xp.median(values[indexes], axis=0)
    return filtered


def spatial_gradient(frame, *, array_module=np, ndimage=ndi):
    """Return the magnitude of the horizontal and vertical 3x3 Sobel filters."""
    xp = array_module
    image = xp.asarray(frame, dtype=xp.float32)
    if image.ndim != 2:
        raise ValueError(f"A moment0 frame must be 2-D, got shape {image.shape}.")
    # ImageJ's Find Edges command combines the Sobel responses as hypot(Gx, Gy).
    # Nearest-edge extension matches its border handling.
    finite_image = xp.nan_to_num(image, nan=0.0, posinf=0.0, neginf=0.0)
    horizontal = ndimage.sobel(finite_image, axis=1, mode="nearest")
    vertical = ndimage.sobel(finite_image, axis=0, mode="nearest")
    return xp.hypot(horizontal, vertical)


def sobel_spatial_gradient(stack, *, array_module=np, ndimage=ndi):
    """Apply Sobel independently to each frame, preserving missing pixels."""
    xp = array_module
    values = xp.asarray(stack, dtype=xp.float32)
    gradients = xp.empty_like(values)
    for index, frame in enumerate(values):
        gradients[index] = spatial_gradient(frame, array_module=xp, ndimage=ndimage)
        gradients[index][~xp.isfinite(frame)] = xp.nan
    return gradients
