"""Pure blood-volume-rate calculations for retinal vessel segments."""

from __future__ import annotations

import numpy as np

from calculations.math import (
    SlidingWindowMethod,
    centered_sliding_window,
    nanmedian,
)
from calculations.topology import annulus_widths_pixels, segment_mask_areas_pixels
from calculations.topology.geometry import AnnulusGeometry

TOTAL_MASKED_EDGES_WINDOW_SIZE = 9
TOTAL_MASKED_EDGES_WINDOW_STRIDE = 1


def mask_derived_lumen_geometry(
    topologies,
    *,
    pixel_size_mm: float,
) -> tuple[tuple[np.ndarray, ...], tuple[np.ndarray, ...], np.ndarray]:
    """Return equivalent diameters, native mask areas, and annulus widths."""

    if not np.isfinite(pixel_size_mm) or pixel_size_mm <= 0:
        raise ValueError("pixel_size_mm must be finite and positive.")
    topologies = tuple(topologies)
    if not topologies:
        raise ValueError("at least one segment topology is required.")

    ring_settings: AnnulusGeometry | None = None
    spatial_shape: tuple[int, int] | None = None
    mask_areas: list[np.ndarray] = []
    for value in topologies:
        topology = getattr(value, "topology", value)
        nested = getattr(topology, "topology", None)
        topology = nested if nested is not None else topology
        settings = topology.ring_settings
        if not isinstance(settings, AnnulusGeometry):
            raise ValueError("segment topology must retain its ring settings.")
        if ring_settings is None:
            ring_settings = settings
        elif settings != ring_settings:
            raise ValueError("artery and vein segment ring settings must match.")
        labels = np.asarray(topology.labels, dtype=np.int32)
        if spatial_shape is None:
            spatial_shape = labels.shape
        elif labels.shape != spatial_shape:
            raise ValueError("artery and vein segment image shapes must match.")
        mask_areas.append(segment_mask_areas_pixels(topology))

    if ring_settings is None or spatial_shape is None:
        raise RuntimeError("segment geometry preparation produced no topology.")
    ring_count = int(mask_areas[0].shape[1])
    if any(area.shape[1] != ring_count for area in mask_areas):
        raise ValueError("artery and vein segment ring counts must match.")
    radial_widths = annulus_widths_pixels(spatial_shape, ring_settings, ring_count)

    diameters: list[np.ndarray] = []
    for area in mask_areas:
        diameter = np.full(area.shape, np.nan, dtype=np.float32)
        valid_widths = radial_widths > 0
        if np.any(valid_widths):
            calculated = (
                area[:, valid_widths].astype(np.float32)
                * np.float32(pixel_size_mm)
                / radial_widths[None, valid_widths]
            )
            diameter[:, valid_widths] = np.where(
                area[:, valid_widths] > 0,
                calculated,
                np.float32(np.nan),
            )
        diameters.append(diameter)
    return tuple(diameters), tuple(mask_areas), radial_widths


def circular_lumen_flow(velocity, diameter_mm) -> np.ndarray:
    """Multiply per-beat velocity by circular lumen area."""

    velocity_tbkr = np.asarray(velocity, dtype=np.float32)
    diameter = np.asarray(diameter_mm, dtype=np.float32)
    if velocity_tbkr.ndim != 4:
        raise ValueError(
            "per-beat segment velocity must have dimensions "
            "(time, beat, branch, radius)."
        )
    if diameter.shape == velocity_tbkr.shape[2:]:
        diameter = diameter[None, None, :, :]
    elif diameter.shape != velocity_tbkr.shape:
        raise ValueError(
            "lumen diameter must have dimensions (branch, radius) or match "
            "per-beat velocity dimensions."
        )
    area_mm2 = np.float32(np.pi / 4.0) * diameter**2
    return (velocity_tbkr * area_mm2).astype(
        np.float32,
        copy=False,
    )


def total_masked_edges_flow(masked_edges) -> np.ndarray:
    """Apply the established temporal, branch, and radius reductions."""

    values = np.asarray(masked_edges, dtype=np.float32)
    if values.ndim != 4:
        raise ValueError(
            "masked-edge blood-volume rate must have dimensions "
            "(time, beat, branch, radius)."
        )
    values = centered_sliding_window(
        values,
        TOTAL_MASKED_EDGES_WINDOW_SIZE,
        SlidingWindowMethod.AVERAGE,
        window_stride=TOTAL_MASKED_EDGES_WINDOW_STRIDE,
        axis=0,
    )

    finite = np.isfinite(values)
    rate_tbr = np.sum(
        np.where(finite, values, np.float32(0.0)),
        axis=2,
        dtype=np.float32,
    )
    rate_tbr[~np.any(finite, axis=2)] = np.nan
    return nanmedian(rate_tbr, axis=2).astype(np.float32, copy=False)


__all__ = [
    "TOTAL_MASKED_EDGES_WINDOW_SIZE",
    "TOTAL_MASKED_EDGES_WINDOW_STRIDE",
    "circular_lumen_flow",
    "mask_derived_lumen_geometry",
    "total_masked_edges_flow",
]
