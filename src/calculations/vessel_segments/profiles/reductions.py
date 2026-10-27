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


def _masked_segments(
    segments,
    segment_masks: np.ndarray | None,
) -> np.ndarray:
    xp = _array_module(segments)
    values = xp.asarray(segments, dtype=xp.float32)
    if segment_masks is None:
        return values

    masks = xp.asarray(segment_masks, dtype=bool)
    if (
        masks.ndim < 2
        or masks.ndim > values.ndim
        or masks.shape[-2:] != values.shape[-2:]
    ):
        raise ValueError(
            "segment_masks must match the segment spatial axes."
        )
    leading_mask_shape = masks.shape[:-2]
    if tuple(values.shape[: len(leading_mask_shape)]) != leading_mask_shape:
        raise ValueError(
            "segment_masks leading axes must match the corresponding segment axes."
        )
    expanded_shape = (
        *leading_mask_shape,
        *((1,) * (values.ndim - masks.ndim)),
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
