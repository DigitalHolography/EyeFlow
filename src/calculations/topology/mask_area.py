"""Native-pixel areas and widths for optic-disc-centered annuli."""

from __future__ import annotations

import numpy as np

from .geometry import AnnulusGeometry, image_half_diagonal, section_bounds


def annulus_widths_pixels(
    image_shape: tuple[int, int],
    settings: AnnulusGeometry,
    ring_count: int,
) -> np.ndarray:
    """Return radial widths, including a clipped final annulus."""

    radius_scale = image_half_diagonal(*image_shape)
    widths = np.zeros(ring_count, dtype=np.float32)
    for ring_index in range(ring_count):
        inner, outer = section_bounds(settings, ring_index)
        widths[ring_index] = np.float32(
            max(outer - inner, 0.0) * radius_scale
        )
    return widths


def segment_mask_areas_pixels(topology) -> np.ndarray:
    """Count labeled native pixels for every branch and annular section.

    The returned array uses the public ``(branch, radius)`` order shared by
    per-beat segment outputs.
    """

    labels = np.asarray(topology.labels, dtype=np.int32)
    annuli = np.asarray(topology.annulus_masks, dtype=bool)
    branch_ids = np.asarray(topology.branch_ids, dtype=np.int32).reshape(-1)
    if labels.ndim != 2:
        raise ValueError("segment topology labels must be a 2-D array.")
    if annuli.ndim != 3 or tuple(annuli.shape[1:]) != labels.shape:
        raise ValueError(
            "segment topology annulus masks must match its label image."
        )
    if np.any(branch_ids < 1):
        raise ValueError("segment topology branch IDs must be positive.")

    areas = np.zeros((branch_ids.size, annuli.shape[0]), dtype=np.int32)
    label_count = max(int(labels.max()) + 1, 1)
    for radius_index, annulus in enumerate(annuli):
        counts = np.bincount(
            labels[annulus].ravel(),
            minlength=label_count,
        )
        for branch_index, branch_id in enumerate(branch_ids):
            label_id = int(branch_id)
            if label_id < counts.size:
                areas[branch_index, radius_index] = np.int32(counts[label_id])
    return areas


__all__ = ["annulus_widths_pixels", "segment_mask_areas_pixels"]
