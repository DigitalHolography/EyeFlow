"""Coordinate layouts used by EyeFlow output datasets."""

from __future__ import annotations

import numpy as np


def serialize_spatial_image(image: np.ndarray) -> np.ndarray:
    """Store an aligned image as lower-left-origin (x, y) coordinates."""
    return np.flip(np.asarray(image), axis=0).T.copy()


def serialize_label_map(image: np.ndarray) -> np.ndarray:
    return serialize_spatial_image(np.asarray(image, dtype=np.int32))


def serialize_segment_velocity(values: np.ndarray) -> np.ndarray:
    """Store (radius, branch, frame) signals as (frame, branch, radius)."""
    return np.asarray(values).transpose(2, 1, 0)


def serialize_segment_masks(masks: np.ndarray) -> np.ndarray:
    """Store (radius, branch, y, x) masks as (x, y, branch, radius)."""
    return np.asarray(masks).transpose(3, 2, 1, 0)
