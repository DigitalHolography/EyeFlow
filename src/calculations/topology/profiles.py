"""Generic profile reductions for uniformly oriented segment maps."""

from __future__ import annotations

import numpy as np

from calculations.compute_backend import optional_cupy_backend

from calculations.math import nanmean_float32


def transverse_profiles(
    segments,
    segment_masks: np.ndarray | None = None,
) -> np.ndarray:
    """Average along Y, retaining X and the input's NumPy/CuPy device."""

    return _profile_nanmean(_masked_segments(segments, segment_masks), axis=-2)


def longitudinal_profiles(
    segments,
    segment_masks: np.ndarray | None = None,
) -> np.ndarray:
    """Average along X, retaining Y and the input's NumPy/CuPy device."""

    return _profile_nanmean(_masked_segments(segments, segment_masks), axis=-1)


def mean_profiles(profiles, *, axis: int) -> np.ndarray:
    """Return the NaN-aware mean of every profile along one named axis."""

    return nanmean_float32(np.asarray(profiles, dtype=np.float32), axis=axis)


def profile_deviation_power(
    profiles,
    mean: np.ndarray | None = None,
    *,
    axis: int,
) -> np.ndarray:
    """Return squared deviations from a supplied or calculated profile mean."""

    values = np.asarray(profiles, dtype=np.float32)
    normalized_axis = int(axis) % values.ndim
    if mean is None:
        profile_mean = mean_profiles(values, axis=normalized_axis)
    else:
        profile_mean = np.asarray(mean, dtype=np.float32)
    expected_shape = (*values.shape[:normalized_axis], *values.shape[normalized_axis + 1 :])
    if profile_mean.shape != expected_shape:
        raise ValueError(
            f"mean must have shape {expected_shape}, got {profile_mean.shape}."
        )
    centered = values - np.expand_dims(profile_mean, axis=normalized_axis)
    return np.square(centered).astype(np.float32, copy=False)


def _masked_segments(
    segments,
    segment_masks: np.ndarray | None,
) -> np.ndarray:
    xp = _array_module(segments)
    values = xp.asarray(segments, dtype=xp.float32)
    if segment_masks is None:
        return values

    masks = xp.asarray(segment_masks, dtype=bool)
    expected_shape = (*values.shape[:2], *values.shape[-2:])
    if masks.shape != expected_shape:
        raise ValueError(
            "segment_masks must match the segment, annulus, branch, and spatial axes."
        )
    expanded_shape = (
        *masks.shape[:2],
        *((1,) * (values.ndim - 4)),
        *masks.shape[-2:],
    )
    return xp.where(masks.reshape(expanded_shape), values, xp.float32(np.nan))


def _array_module(values):
    backend = optional_cupy_backend()
    if backend is not None and isinstance(values, backend.cupy.ndarray):
        return backend.cupy
    return np


def _profile_nanmean(values, *, axis: int):
    """Reduce on the input device without downloading segment movies."""
    xp = _array_module(values)
    if xp is np:
        return nanmean_float32(values, axis=axis)
    finite = xp.isfinite(values)
    count = xp.sum(finite, axis=axis, dtype=xp.int32)
    total = xp.sum(xp.where(finite, values, xp.float32(0)), axis=axis, dtype=xp.float32)
    return xp.where(
        count > 0, total / xp.maximum(count, 1), xp.float32(np.nan),
    ).astype(xp.float32, copy=False)
