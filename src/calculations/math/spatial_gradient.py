"""Temporal averaging, Sobel magnitude and ImageJ spatial filters for stacks."""

from __future__ import annotations

from functools import lru_cache
from math import ceil, exp, log, sqrt

import numpy as np
from scipy import ndimage as ndi

TEMPORAL_MOVING_AVERAGE_WINDOW = 7
GAUSSIAN_BLUR_RADIUS = 4.0
UNSHARP_MASK_RADIUS = 4.0
UNSHARP_MASK_WEIGHT = 0.6
_IMAGEJ_GAUSSIAN_ACCURACY = 0.0002
_IMAGEJ_UNSHARP_ACCURACY = 0.01


def moving_avg_window(stack, *, window=TEMPORAL_MOVING_AVERAGE_WINDOW, array_module=np):
    """Average each pixel in time, truncating centered windows at the edges.

    Windows contain only available frames, with no padding or duplication.
    Any NaN in a pixel's window propagates to that pixel's average.
    """
    if window < 1 or window % 2 == 0:
        raise ValueError("moving average window must be a positive odd integer.")
    xp = array_module
    values = xp.asarray(stack, dtype=xp.float32)
    filtered = xp.empty_like(values)
    half_window = window // 2
    for index in range(len(values)):
        start = max(0, index - half_window)
        stop = min(len(values), index + half_window + 1)
        filtered[index] = xp.mean(values[start:stop], axis=0)
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


@lru_cache(maxsize=2)
def _imagej_gaussian_kernel(sigma, accuracy):
    """Build ImageJ's finite Gaussian kernel for our radii of 6 and 8 pixels.

    GaussianBlur.makeGaussianKernel smooths the kernel tails before normalizing.
    At these radii ImageJ uses direct convolution, without downscaling the image.
    Reference: https://github.com/imagej/ImageJ/blob/master/ij/plugin/filter/GaussianBlur.java
    """
    radius = ceil(sigma * sqrt(-2.0 * log(accuracy))) + 1
    half = np.asarray(
        [exp(-0.5 * index * index / sigma / sigma) for index in range(radius)],
        dtype=np.float32,
    )
    sqrt_slope = float("inf")
    for index in range(radius - 1, radius // 2 - 1, -1):
        slope = sqrt(float(half[index])) / (radius - index)
        if slope >= sqrt_slope:
            break
        sqrt_slope = slope
    for tail_index in range(index + 2, radius):
        half[tail_index] = (radius - tail_index) ** 2 * sqrt_slope**2
    normalization = float(half[0]) + sum(2.0 * float(value) for value in half[1:])
    half = (half.astype(np.float64) / normalization).astype(np.float32)
    return np.concatenate((half[:0:-1], half))


def _zero_valid_border(frame, array_module):
    """Zero the outermost rows and columns without replacing propagated NaNs."""
    xp = array_module
    for edge in (frame[0], frame[-1], frame[:, 0], frame[:, -1]):
        edge[...] = xp.where(xp.isnan(edge), xp.nan, xp.float32(0.0))


def gaussian2d_blur(stack, *, array_module=np, ndimage=ndi):
    """Apply ImageJ's Gaussian Blur independently to each spatial frame.

    Uses sigma=6 and ImageJ's float-image accuracy of 0.0002 with nearest-edge
    extension. NaNs propagate through the spatial kernel; valid output pixels
    on the outermost rows and columns are set to zero.
    """
    xp = array_module
    values = xp.asarray(stack, dtype=xp.float32)
    kernel = xp.asarray(
        _imagej_gaussian_kernel(GAUSSIAN_BLUR_RADIUS, _IMAGEJ_GAUSSIAN_ACCURACY),
        dtype=xp.float32,
    )
    filtered = xp.empty_like(values)
    for index, frame in enumerate(values):
        blurred = ndimage.correlate1d(frame, kernel, axis=1, mode="nearest")
        blurred = ndimage.correlate1d(blurred, kernel, axis=0, mode="nearest")
        _zero_valid_border(blurred, xp)
        filtered[index] = blurred
    return filtered


def unsharpen(stack, *, array_module=np, ndimage=ndi):
    """Apply ImageJ's unsharp mask independently to each spatial frame.

    Uses sigma=8, mask weight=0.6 and nearest-edge Gaussian extension, with
    ImageJ's formula (input - weight * blurred) / (1 - weight). NaNs propagate
    through the spatial kernel, negative results are clipped to zero, and
    finite pixels on the outermost rows and columns are set to zero.
    """
    xp = array_module
    values = xp.asarray(stack, dtype=xp.float32)
    kernel = xp.asarray(
        _imagej_gaussian_kernel(UNSHARP_MASK_RADIUS, _IMAGEJ_UNSHARP_ACCURACY),
        dtype=xp.float32,
    )
    filtered = xp.empty_like(values)
    weight = xp.float32(UNSHARP_MASK_WEIGHT)
    for index, frame in enumerate(values):
        blurred = ndimage.correlate1d(frame, kernel, axis=1, mode="nearest")
        blurred = ndimage.correlate1d(blurred, kernel, axis=0, mode="nearest")
        sharpened = (frame - weight * blurred) / (xp.float32(1.0) - weight)
        xp.maximum(sharpened, xp.float32(0.0), out=sharpened)
        _zero_valid_border(sharpened, xp)
        filtered[index] = sharpened
    return filtered
