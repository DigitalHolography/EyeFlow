"""Construct native branch/annulus geometry and local vessel masks."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .branch_identity import BranchIdentityResult, label_vessel_branches
from .geometry import AnnulusGeometry, image_half_diagonal, section_masks
from .mask_area import annulus_widths_pixels
from .optic_disc import OpticDisc
from .windows import centered_window_bounds, window_target_slices


@dataclass(frozen=True)
class SegmentTopology:
    """Map-independent locations of vessel segments.

    Arrays indexed by segment use ``(annulus, branch, ...)`` ordering.
    ``centerline`` is the full-frame skeleton used during branch identification.
    ``segment_masks`` are fixed-size local masks rather than full-frame masks.
    """

    optic_disc_center_xy: tuple[float, float]
    branches: BranchIdentityResult
    annulus_masks: np.ndarray
    segment_masks: np.ndarray
    segment_centers_xy: np.ndarray
    window_bounds_xyxy: np.ndarray
    delta_radius: np.ndarray | None = None
    optic_disc_mask: np.ndarray | None = None
    ring_settings: AnnulusGeometry | None = None

    @property
    def spatial_shape(self) -> tuple[int, int]:
        return tuple(int(size) for size in self.branches.labels.shape)

    @property
    def labels(self) -> np.ndarray:
        return self.branches.labels

    @property
    def centerline(self) -> np.ndarray:
        return self.branches.centerline

    @property
    def branch_ids(self) -> np.ndarray:
        return self.branches.branch_ids

    @property
    def branch_identity(self) -> BranchIdentityResult:
        """Compatibility name for the authoritative branch composition."""

        return self.branches

    @property
    def window_side_pixels(self) -> int:
        return int(self.segment_masks.shape[-1])

    @property
    def segment_shape(self) -> tuple[int, int]:
        return tuple(int(size) for size in self.segment_centers_xy.shape[:2])

    @property
    def valid_segments(self) -> np.ndarray:
        return np.all(np.isfinite(self.segment_centers_xy), axis=-1)


def build_segment_topology(
    vessel_mask,
    optic_disc: OpticDisc,
    settings: AnnulusGeometry,
    *,
    window_size_percentile_kept: float = 0.95,
    window_side_pixels: int | None = None,
) -> SegmentTopology:
    """Find branch segments around an optic disc.

    Args:
        vessel_mask: Two-dimensional mask containing one vessel class.
        optic_disc: Authoritative optic-disc geometry in the vessel-mask frame.
        settings: Annulus placement and sampling settings.
        window_size_percentile_kept: Fraction of segment widths and heights
            represented by the automatically selected square window.
        window_side_pixels: Explicit odd window side, overriding the percentile.

    Returns:
        Geometry reusable to extract the same segments from any retinal map.
    """

    vessel = np.asarray(vessel_mask, dtype=bool)
    assert vessel.ndim == 2
    return _build_segment_topology(
        vessel,
        optic_disc,
        settings,
        window_size_percentile_kept=window_size_percentile_kept,
        window_side_pixels=window_side_pixels,
    )


def competing_segment_masks(
    topology: SegmentTopology,
    vessel_mask,
    other_vessel_mask=None,
) -> np.ndarray:
    """Return local masks for vessels other than each segment's own branch."""

    vessels = np.asarray(vessel_mask, dtype=bool)
    if vessels.shape != topology.spatial_shape:
        raise ValueError(
            "vessel_mask must match the topology spatial shape, got "
            f"{vessels.shape}."
        )
    other_vessels = (
        np.zeros_like(vessels)
        if other_vessel_mask is None
        else np.asarray(other_vessel_mask, dtype=bool)
    )
    if other_vessels.shape != topology.spatial_shape:
        raise ValueError(
            "other_vessel_mask must match the topology spatial shape, got "
            f"{other_vessels.shape}."
        )
    ring_count, branch_count = topology.segment_centers_xy.shape[:2]
    side = topology.window_side_pixels
    masks = np.zeros((ring_count, branch_count, side, side), dtype=bool)
    if side == 0:
        return masks

    for ring_index, branch_index in np.argwhere(topology.valid_segments):
        index = (int(ring_index), int(branch_index))
        bounds = topology.window_bounds_xyxy[index]
        center = topology.segment_centers_xy[index]
        target_y, target_x = window_target_slices(bounds, center, side)
        x_start, x_stop, y_start, y_stop = bounds
        branch_id = int(topology.branch_ids[index[1]])
        source_y = slice(int(y_start), int(y_stop))
        source_x = slice(int(x_start), int(x_stop))
        own_branch = topology.labels[source_y, source_x] == branch_id
        masks[index][target_y, target_x] = (
            (vessels[source_y, source_x] & ~own_branch)
            | (other_vessels[source_y, source_x] & ~own_branch)
        )
    return masks


def resize_segment_topology_windows(
    topology: SegmentTopology,
    window_side_pixels: int,
) -> SegmentTopology:
    """Rebuild local windows without repeating branch identification."""

    side = int(window_side_pixels)
    if side < 0 or (side > 0 and side % 2 == 0):
        raise ValueError("window_side_pixels must be zero or a positive odd integer.")
    if side == topology.window_side_pixels:
        return topology

    ring_count, branch_count = topology.segment_centers_xy.shape[:2]
    bounds = np.full((ring_count, branch_count, 4), -1, dtype=np.int32)
    masks = np.zeros((ring_count, branch_count, side, side), dtype=bool)
    if side:
        for ring_index, branch_index in np.argwhere(topology.valid_segments):
            center = tuple(
                int(value)
                for value in topology.segment_centers_xy[ring_index, branch_index]
            )
            segment_bounds = centered_window_bounds(
                topology.spatial_shape,
                center,
                side,
            )
            bounds[ring_index, branch_index] = segment_bounds
            target_y, target_x = window_target_slices(segment_bounds, center, side)
            x_start, x_stop, y_start, y_stop = segment_bounds
            branch_id = int(topology.branch_ids[branch_index])
            full_mask = topology.annulus_masks[ring_index] & (
                topology.labels == branch_id
            )
            masks[ring_index, branch_index, target_y, target_x] = full_mask[
                y_start:y_stop,
                x_start:x_stop,
            ]

    return SegmentTopology(
        optic_disc_center_xy=topology.optic_disc_center_xy,
        optic_disc_mask=topology.optic_disc_mask,
        ring_settings=topology.ring_settings,
        branches=topology.branches,
        annulus_masks=topology.annulus_masks,
        segment_masks=masks,
        segment_centers_xy=topology.segment_centers_xy,
        window_bounds_xyxy=bounds,
        delta_radius=topology.delta_radius,
    )


def _build_segment_topology(
    vessel_mask: np.ndarray,
    optic_disc: OpticDisc,
    settings: AnnulusGeometry,
    *,
    window_size_percentile_kept: float,
    window_side_pixels: int | None,
) -> SegmentTopology:
    branches = label_vessel_branches(
        vessel_mask,
        optic_disc,
        settings,
    )
    fallback_radius = (
        float(settings.inner_radius_frac)
        * max(image_half_diagonal(*vessel_mask.shape), 1.0)
    )
    optic_disc_mask = optic_disc.centered_circle_mask_for(
        vessel_mask.shape,
        fallback_radius_pixels=fallback_radius,
    )
    optic_disc_center_xy = optic_disc.center
    centerline = branches.centerline
    annuli = section_masks(vessel_mask.shape, optic_disc_center_xy, settings)
    annuli &= ~optic_disc_mask[None, ...]
    side = (
        _segment_window_side(
            branches.labels,
            branches.branch_ids,
            annuli,
            window_size_percentile_kept,
        )
        if window_side_pixels is None
        else int(window_side_pixels)
    )
    if side < 0 or (side > 0 and side % 2 == 0):
        raise ValueError("window_side_pixels must be zero or a positive odd integer.")

    ring_count = int(annuli.shape[0])
    branch_count = int(branches.branch_ids.size)
    delta_radius = annulus_widths_pixels(
        vessel_mask.shape,
        settings,
        ring_count,
    )
    centers = np.full((ring_count, branch_count, 2), np.nan, dtype=np.float32)
    bounds = np.full((ring_count, branch_count, 4), -1, dtype=np.int32)
    masks = np.zeros((ring_count, branch_count, side, side), dtype=bool)
    for ring_index, annulus in enumerate(annuli):
        for branch_index, branch_id in enumerate(branches.branch_ids):
            mask = annulus & (branches.labels == int(branch_id))
            segment_centerline = centerline & mask
            center = _segment_center_xy(mask, segment_centerline)
            if center is None:
                continue
            segment_bounds = centered_window_bounds(vessel_mask.shape, center, side)
            centers[ring_index, branch_index] = center
            bounds[ring_index, branch_index] = segment_bounds
            target_y, target_x = window_target_slices(segment_bounds, center, side)
            x_start, x_stop, y_start, y_stop = segment_bounds
            masks[ring_index, branch_index, target_y, target_x] = mask[
                y_start:y_stop,
                x_start:x_stop,
            ]

    return SegmentTopology(
        optic_disc_center_xy=tuple(float(value) for value in optic_disc_center_xy),
        optic_disc_mask=optic_disc_mask.copy(),
        ring_settings=settings,
        branches=branches,
        annulus_masks=annuli,
        segment_masks=masks,
        segment_centers_xy=centers,
        window_bounds_xyxy=bounds,
        delta_radius=delta_radius,
    )


def _segment_center_xy(
    segment_mask: np.ndarray,
    segment_centerline: np.ndarray,
) -> tuple[int, int] | None:
    source = segment_centerline if np.any(segment_centerline) else segment_mask
    if not np.any(source):
        return None
    point_y, point_x = np.nonzero(source)
    center_x = int(np.floor(np.median(point_x) + 0.5))
    center_y = int(np.floor(np.median(point_y) + 0.5))
    return center_x, center_y


def _segment_window_side(
    labels: np.ndarray,
    branch_ids: np.ndarray,
    annuli: np.ndarray,
    percentile_kept: float,
) -> int:
    percentile = float(percentile_kept)
    if not 0.0 < percentile <= 1.0:
        raise ValueError("window_size_percentile_kept must be in (0, 1].")

    widths: list[int] = []
    heights: list[int] = []
    for annulus in annuli:
        for branch_id in branch_ids:
            segment_y, segment_x = np.nonzero(annulus & (labels == int(branch_id)))
            if segment_x.size == 0:
                continue
            widths.append(int(segment_x.max() - segment_x.min() + 1))
            heights.append(int(segment_y.max() - segment_y.min() + 1))
    if not widths:
        return 0

    width = int(np.quantile(widths, percentile, method="higher"))
    height = int(np.quantile(heights, percentile, method="higher"))
    side = max(width, height)
    return side if side % 2 == 1 else side + 1
