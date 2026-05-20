"""Branch labeling ported from CrossSection/labelVesselBranches.m."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import ndimage as ndi
from skimage.measure import label as label_components
from skimage.morphology import skeletonize
from skimage.segmentation import watershed

from .segment_geometry import SegmentRingSettings, annulus_mask


LOW_RES_SMALL_BRANCH_PIXELS = 10
LOW_RES_PER_CIRCLE_MIN_AREA_PIXELS = 51
STREL_SIZE = 3


@dataclass(frozen=True)
class BranchIdentityStages:
    # Matlab-style intermediate arrays used to inspect branch labeling.
    vessel: np.ndarray
    section: np.ndarray
    skeleton: np.ndarray
    branch_points: np.ndarray
    cleaned_skeleton: np.ndarray
    marker_labels: np.ndarray
    distance_topography: np.ndarray
    imposed_minima_topography: np.ndarray
    watershed_labels: np.ndarray
    annulus_refined_labels: np.ndarray
    per_circle_cleaned_labels: np.ndarray


@dataclass(frozen=True)
class BranchIdentityResult:
    # Final branch labels plus optional diagnostic stages from the same run.
    labels: np.ndarray
    branch_ids: np.ndarray
    stages: BranchIdentityStages


def label_vessel_branches(
    vessel_mask,
    optic_disc_center,
    settings: SegmentRingSettings,
    *,
    small_branch_pixels: int = LOW_RES_SMALL_BRANCH_PIXELS,
    strel_size: int = STREL_SIZE,
) -> BranchIdentityResult:
    vessel = np.asarray(vessel_mask, dtype=bool)
    if vessel.ndim != 2:
        raise ValueError("vessel_mask must be a 2-D array.")

    stages = _branch_identity_stages(
        vessel,
        optic_disc_center,
        settings,
        small_branch_pixels=small_branch_pixels,
        strel_size=strel_size,
    )
    labels = stages.per_circle_cleaned_labels
    return BranchIdentityResult(
        labels.astype(np.int32, copy=False),
        np.arange(1, int(labels.max()) + 1, dtype=np.int32),
        stages,
    )


def _branch_identity_stages(
    vessel: np.ndarray,
    optic_disc_center,
    settings: SegmentRingSettings,
    *,
    small_branch_pixels: int = LOW_RES_SMALL_BRANCH_PIXELS,
    strel_size: int = STREL_SIZE,
) -> BranchIdentityStages:
    section = annulus_mask(
        vessel.shape,
        optic_disc_center,
        settings.inner_radius_frac,
        settings.outer_radius_frac,
    )
    skeleton, branch_points, cleaned_skeleton, marker_labels = _seed_labels(
        vessel,
        int(small_branch_pixels),
        int(strel_size),
    )
    distance_topography = _distance_topography(vessel, section)
    imposed_topography = _impose_marker_minima(distance_topography, marker_labels > 0)
    watershed_labels = _watershed_labels(vessel, imposed_topography, marker_labels)
    annulus_refined = _annulus_refined_labels(watershed_labels, section)
    cleaned_labels = _per_circle_cleaned_labels(
        annulus_refined,
        optic_disc_center,
        settings,
    )
    return BranchIdentityStages(
        vessel,
        section,
        skeleton,
        branch_points,
        cleaned_skeleton,
        marker_labels,
        distance_topography,
        imposed_topography,
        watershed_labels,
        annulus_refined,
        cleaned_labels,
    )


def _seed_labels(
    vessel: np.ndarray,
    small_branch_pixels: int,
    strel_size: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    skeleton = skeletonize(vessel)
    branch_points = _branch_points(skeleton)
    cleaned_skeleton = _clean_skeleton(
        skeleton,
        branch_points,
        small_branch_pixels,
        strel_size,
    )
    return skeleton, branch_points, cleaned_skeleton, _marker_labels(cleaned_skeleton)


def _clean_skeleton(
    skeleton: np.ndarray,
    branch_points: np.ndarray,
    small_branch_pixels: int,
    strel_size: int = STREL_SIZE,
) -> np.ndarray:
    structure = np.ones((max(1, strel_size), max(1, strel_size)), dtype=bool)
    cleaned = skeleton & ~ndi.binary_dilation(branch_points, structure=structure)
    return _remove_small(cleaned, max(0, small_branch_pixels))


def _branch_points(skeleton: np.ndarray) -> np.ndarray:
    padded = np.pad(skeleton.astype(np.uint8), 1)
    neighbors = np.stack(
        (
            padded[:-2, 1:-1],
            padded[:-2, 2:],
            padded[1:-1, 2:],
            padded[2:, 2:],
            padded[2:, 1:-1],
            padded[2:, :-2],
            padded[1:-1, :-2],
            padded[:-2, :-2],
        ),
        axis=0,
    )
    neighbors = neighbors.astype(np.int8, copy=False)
    transitions = np.sum(np.abs(neighbors - np.roll(neighbors, -1, axis=0)), axis=0)
    return skeleton & ((transitions // 2) > 2)


def _marker_labels(skeleton: np.ndarray) -> np.ndarray:
    return label_components(skeleton, connectivity=2).astype(np.int32, copy=False)


def _distance_topography(vessel: np.ndarray, section: np.ndarray) -> np.ndarray:
    topography = -ndi.distance_transform_edt(vessel).astype(np.float32, copy=False)
    topography[~(vessel & section)] = -np.inf
    return topography


def _impose_marker_minima(
    topography: np.ndarray,
    marker_mask: np.ndarray,
) -> np.ndarray:
    imposed = np.full(topography.shape, np.inf, dtype=np.float32)
    imposed[marker_mask] = -np.inf
    return imposed


def _watershed_labels(
    vessel: np.ndarray,
    topography: np.ndarray,
    marker_labels: np.ndarray,
) -> np.ndarray:
    if int(marker_labels.max()) == 0:
        return np.zeros(vessel.shape, dtype=np.int32)
    labels = watershed(
        topography,
        markers=marker_labels,
        watershed_line=True,
    )
    labels = (labels * vessel).astype(np.int32, copy=False)
    return _split_touching_watershed_labels(labels)


def _split_touching_watershed_labels(labels: np.ndarray) -> np.ndarray:
    separated = np.asarray(labels, dtype=np.int32).copy()
    separated[_different_label_contacts(separated)] = 0
    return separated


def _different_label_contacts(labels: np.ndarray) -> np.ndarray:
    contacts = np.zeros(labels.shape, dtype=bool)
    for shifted in _neighbor_labels(labels):
        contacts |= (labels > 0) & (shifted > 0) & (labels != shifted)
    return contacts


def _neighbor_labels(labels: np.ndarray):
    padded = np.pad(labels, 1)
    return (
        padded[:-2, 1:-1],
        padded[:-2, 2:],
        padded[1:-1, 2:],
        padded[2:, 2:],
        padded[2:, 1:-1],
        padded[2:, :-2],
        padded[1:-1, :-2],
        padded[:-2, :-2],
    )


def _annulus_refined_labels(
    labels: np.ndarray,
    section: np.ndarray,
) -> np.ndarray:
    masked_labels = labels * section.astype(np.int32)
    return label_components(masked_labels > 0, connectivity=2).astype(np.int32, copy=False)


def _per_circle_cleaned_labels(
    labels: np.ndarray,
    optic_disc_center,
    settings: SegmentRingSettings,
) -> np.ndarray:
    cleaned = np.zeros(labels.shape, dtype=bool)
    min_area = LOW_RES_PER_CIRCLE_MIN_AREA_PIXELS
    for ring_index in range(settings.ring_count):
        ring_inner = settings.inner_radius_frac + ring_index * settings.ring_width_frac
        ring = annulus_mask(
            labels.shape,
            optic_disc_center,
            ring_inner,
            ring_inner + settings.ring_width_frac,
        )
        cleaned |= _remove_small((labels > 0) & ring, min_area)
    return label_components(cleaned, connectivity=2).astype(np.int32, copy=False)


def _remove_small(mask: np.ndarray, min_area: int) -> np.ndarray:
    if min_area <= 1:
        return mask
    labeled, _ = ndi.label(mask, structure=np.ones((3, 3), dtype=bool))
    sizes = np.bincount(labeled.reshape(-1))
    keep = sizes >= int(min_area)
    keep[0] = False
    return keep[labeled]
