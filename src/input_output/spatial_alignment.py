"""Shared spatial alignment of source arrays to HoloDoppler coordinates."""

from __future__ import annotations

import numpy as np


def align_mask(
    array,
    shape: tuple[int, int],
    name: str,
) -> tuple[np.ndarray, bool]:
    value = np.asarray(array)
    if value.shape[-2:] == shape:
        return np.asarray(value, dtype=bool), False
    if value.shape[-2:] == shape[::-1]:
        return np.asarray(np.swapaxes(value, -1, -2), dtype=bool), True
    raise ValueError(f"{name} spatial shape {value.shape[-2:]} does not match {shape}.")


def apply_mask_alignment(
    array,
    shape: tuple[int, int],
    name: str,
    *,
    swapped: bool,
) -> np.ndarray:
    return np.asarray(
        apply_spatial_alignment(array, shape, name, swapped=swapped),
        dtype=bool,
    )


def apply_spatial_alignment(
    array,
    shape: tuple[int, int],
    name: str,
    *,
    swapped: bool,
) -> np.ndarray:
    value = np.asarray(array)
    aligned = np.swapaxes(value, -1, -2) if swapped else value
    if aligned.shape[-2:] != shape:
        raise ValueError(
            f"{name} does not share the selected DopplerView spatial orientation; "
            f"aligned shape {aligned.shape[-2:]} does not match {shape}."
        )
    return aligned


def apply_optional_mask_alignment(
    array,
    shape: tuple[int, int],
    name: str,
    *,
    swapped: bool,
) -> np.ndarray | None:
    if array is None:
        return None
    return apply_spatial_alignment(array, shape, name, swapped=swapped)
