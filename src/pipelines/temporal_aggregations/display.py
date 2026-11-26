"""HoloDoppler-style PNG scaling for temporal image statistics."""

import numpy as np


def png_pixels(image):
    """Rescale a grayscale or RGB image to 8-bit using one shared range."""
    finite = np.isfinite(image)
    low = np.min(image, where=finite, initial=np.inf)
    high = np.max(image, where=finite, initial=-np.inf)
    if high <= low:
        return np.zeros(image.shape, dtype=np.uint8)
    scaled = np.subtract(image, low, out=np.zeros_like(image), where=finite) / (high - low)
    return np.floor(scaled * 255 + 0.5).astype(np.uint8)
