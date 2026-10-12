"""Build shared spatial sampling plans and resolve reference-based orientations."""

from __future__ import annotations

from collections.abc import Mapping, MutableMapping
from dataclasses import replace

import numpy as np
from scipy import ndimage as ndi

from calculations.topology.geometry import AnnulusGeometry, image_half_diagonal
from calculations.topology.optic_disc import OpticDisc
from calculations.topology.orientation import determine_segment_rotations
from calculations.topology.segments import (
    build_segment_topology,
    resize_segment_topology_windows,
)
from utils.logger import Logger

from .cache import TopologyCacheKey, topology_cache_key
from .extraction import extract_segment
from .models import SegmentSamplingPlan
from .transforms import interpolate_segment_masks, interpolate_segments, rotate_segment_masks


def prepare_sampling_plan(
    vessel_mask,
    optic_disc: OpticDisc,
    settings: AnnulusGeometry,
    *,
    output_side_pixels: int = 128,
    window_size_percentile_kept: float = 0.95,
    window_side_pixels: int | None = None,
) -> SegmentSamplingPlan:
    """Build segment geometry, orientations, and uniform masks once."""

    topology = build_segment_topology(
        vessel_mask,
        optic_disc,
        settings,
        window_size_percentile_kept=window_size_percentile_kept,
        window_side_pixels=window_side_pixels,
    )
    rotation_degrees = determine_segment_rotations(topology)
    interpolated_masks = interpolate_segment_masks(
        topology.segment_masks,
        output_side_pixels,
    )
    return SegmentSamplingPlan(
        native=topology,
        rotation_degrees=rotation_degrees,
        interpolated_masks=interpolated_masks,
        rotated_masks=rotate_segment_masks(
            interpolated_masks,
            rotation_degrees,
        ),
    )


def prepare_sampling_plans(
    vessel_masks: Mapping[str, object],
    optic_disc: OpticDisc,
    settings: AnnulusGeometry,
    *,
    source_id: str,
    cache: MutableMapping[TopologyCacheKey, object] | None = None,
    output_side_pixels: int = 128,
    window_size_percentile_kept: float = 0.95,
    window_side_pixels: int | None = None,
) -> dict[str, SegmentSamplingPlan]:
    """Prepare consistently sized topology for every named vessel mask."""

    masks = {
        str(name): np.asarray(mask, dtype=bool)
        for name, mask in vessel_masks.items()
    }
    if not masks:
        return {}
    image_shape = next(iter(masks.values())).shape
    if any(mask.shape != image_shape for mask in masks.values()):
        raise ValueError("All vessel masks must have the same spatial shape.")
    disc = _topology_optic_disc_mask(optic_disc, image_shape, settings)

    if window_side_pixels is not None:
        return {
            name: _cached_topology(
                name,
                mask,
                optic_disc,
                settings,
                source_id=source_id,
                cache=cache,
                output_side_pixels=output_side_pixels,
                window_size_percentile_kept=window_size_percentile_kept,
                window_side_pixels=window_side_pixels,
            )
            for name, mask in masks.items()
        }

    initial = {
        name: _cached_topology(
            name,
            mask,
            optic_disc,
            settings,
            source_id=source_id,
            cache=cache,
            output_side_pixels=output_side_pixels,
            window_size_percentile_kept=window_size_percentile_kept,
            window_side_pixels=None,
        )
        for name, mask in masks.items()
    }
    shared_side = _shared_window_side(initial, window_size_percentile_kept)
    prepared_topologies = {
        name: (
            topology
            if topology.native.window_side_pixels == shared_side
            else _resize_prepared_topology(
                topology,
                shared_side,
                output_side_pixels,
            )
        )
        for name, topology in initial.items()
    }
    if cache is not None:
        for name, prepared in prepared_topologies.items():
            key = topology_cache_key(
                source_id, name, masks[name], disc, settings,
                optic_disc_center=optic_disc.center,
                output_side_pixels=output_side_pixels,
                window_size_percentile_kept=window_size_percentile_kept,
                window_side_pixels=None,
            )
            # Cache the final joint-window geometry, not just the initial
            # per-vessel geometry, so every pipeline reuses the same objects.
            cache[key] = prepared
    Logger.log(f"Prepared topology: vessels={tuple(masks)}, shared_window={shared_side}px.")
    return prepared_topologies


def _resize_prepared_topology(
    prepared: SegmentSamplingPlan,
    window_side_pixels: int,
    output_side_pixels: int,
) -> SegmentSamplingPlan:
    topology = resize_segment_topology_windows(
        prepared.native,
        window_side_pixels,
    )
    interpolated_masks = interpolate_segment_masks(
        topology.segment_masks,
        output_side_pixels,
    )
    return SegmentSamplingPlan(
        native=topology,
        rotation_degrees=prepared.rotation_degrees,
        interpolated_masks=interpolated_masks,
        rotated_masks=rotate_segment_masks(
            interpolated_masks,
            prepared.rotation_degrees,
        ),
    )


def resolve_segment_rotations(
    prepared: SegmentSamplingPlan,
    reference_map,
    *,
    working_memory_mb: float = 512.0,
) -> SegmentSamplingPlan:
    """Resolve indeterminate centerline angles once from a mean reference map.

    The projection-score fallback is intentionally limited to segments whose
    topology centerline cannot determine an orientation. Segments remain
    unresolved when both methods fail and are skipped by segment iterators.
    """

    rotations = prepared.rotation_degrees.copy()
    unresolved = np.argwhere(
        prepared.native.valid_segments & ~np.isfinite(rotations)
    )
    if len(unresolved) == 0:
        return prepared
    frame_count = int(reference_map.shape[0])
    if frame_count < 1:
        return prepared
    side = max(prepared.native.window_side_pixels, 1)
    budget = int(float(working_memory_mb) * 1024**2)
    if budget <= 0:
        raise ValueError("working_memory_mb must be finite and positive.")
    frames_per_block = max(1, min(frame_count, budget // max(side * side * 8, 1)))
    for ring, branch in unresolved:
        total = np.zeros((side, side), dtype=np.float64)
        count = np.zeros((side, side), dtype=np.int64)
        for start in range(0, frame_count, frames_per_block):
            extracted = extract_segment(
                _frame_slice(reference_map, start, min(start + frames_per_block, frame_count)),
                prepared.native,
                int(ring),
                int(branch),
            )
            finite = np.isfinite(extracted)
            total += np.sum(np.where(finite, extracted, 0.0), axis=0, dtype=np.float64)
            count += np.sum(finite, axis=0, dtype=np.int64)
        mean_native = np.full((side, side), np.nan, dtype=np.float32)
        np.divide(total, count, out=mean_native, where=count > 0)
        mean_image = interpolate_segments(
            mean_native,
            prepared.interpolated_masks.shape[-1],
        )
        mask = prepared.interpolated_masks[int(ring), int(branch)]
        mean_image[~mask] = np.nan
        angle = _projection_rotation_angle(
            mean_image,
            prepared.native.segment_centers_xy[int(ring), int(branch)],
            prepared.native.optic_disc_center_xy,
        )
        if np.isfinite(angle):
            rotations[int(ring), int(branch)] = np.float32(angle)

    return replace(
        prepared,
        rotation_degrees=rotations,
        rotated_masks=rotate_segment_masks(prepared.interpolated_masks, rotations),
    )


def _projection_rotation_angle(
    mean_image: np.ndarray,
    segment_center_xy,
    optic_disc_center_xy,
) -> float:
    image = np.asarray(mean_image, dtype=np.float32).copy()
    yy, xx = np.indices(image.shape, dtype=np.float32)
    cy = np.float32((image.shape[0] - 1) / 2.0)
    cx = np.float32((image.shape[1] - 1) / 2.0)
    image[(yy - cy) ** 2 + (xx - cx) ** 2 > min(cx, cy) ** 2] = np.nan
    image = np.nan_to_num(image, nan=0.0)
    image[image < 0] = 0
    if not np.any(image > 0):
        return float("nan")
    loc_x, loc_y = (float(value) for value in segment_center_xy)
    disc_x, disc_y = (float(value) for value in optic_disc_center_xy)
    alpha = np.degrees(np.arctan2(loc_y - disc_y, loc_x - disc_x))
    beta = int(np.mod(90 + np.floor(np.mod(alpha, 360.0) + 0.5), 180))
    angles = np.mod(np.arange(beta - 90, beta + 91), 180)
    scores = np.asarray([_projection_score(image, angle) for angle in angles])
    if not np.any(np.isfinite(scores)):
        return float("nan")
    return float(angles[int(np.nanargmax(scores))])


def _projection_score(image: np.ndarray, angle: float) -> float:
    rotated = ndi.rotate(
        image,
        angle,
        reshape=False,
        order=0,
        mode="constant",
        cval=0.0,
    )
    start = max(int(np.floor(rotated.shape[0] / 3)) - 1, 0)
    stop = int(np.ceil(2 * rotated.shape[0] / 3))
    projection = np.sum(rotated[start:stop], axis=0, dtype=np.float32)
    total = np.sum(projection, dtype=np.float32)
    return float("nan") if total <= 0 else float(np.max(projection) / total)


def _frame_slice(data_map, start: int, stop: int):
    slices = [slice(None)] * len(data_map.shape)
    slices[0] = slice(int(start), int(stop))
    return data_map[tuple(slices)]


def _cached_topology(
    vessel_name: str,
    vessel_mask: np.ndarray,
    optic_disc: OpticDisc,
    settings: AnnulusGeometry,
    *,
    source_id: str,
    cache: MutableMapping[TopologyCacheKey, object] | None,
    output_side_pixels: int,
    window_size_percentile_kept: float,
    window_side_pixels: int | None,
) -> SegmentSamplingPlan:
    optic_disc_mask = _topology_optic_disc_mask(
        optic_disc,
        vessel_mask.shape,
        settings,
    )
    key = topology_cache_key(
        source_id,
        vessel_name,
        vessel_mask,
        optic_disc_mask,
        settings,
        optic_disc_center=optic_disc.center,
        output_side_pixels=output_side_pixels,
        window_size_percentile_kept=window_size_percentile_kept,
        window_side_pixels=window_side_pixels,
    )
    if cache is not None:
        found = cache.get(key)
        if found is not None:
            if not isinstance(found, SegmentSamplingPlan):
                raise TypeError("Topology cache values must be SegmentSamplingPlan.")
            Logger.log(f"Topology cache hit: {vessel_name}; reusing prepared topology.")
            return found

    cache_status = "miss" if cache is not None else "disabled"
    Logger.log(
        f"Topology cache {cache_status}: {vessel_name}; preparing topology."
    )
    prepared = prepare_sampling_plan(
        vessel_mask,
        optic_disc,
        settings,
        output_side_pixels=output_side_pixels,
        window_size_percentile_kept=window_size_percentile_kept,
        window_side_pixels=window_side_pixels,
    )
    if cache is not None:
        cache[key] = prepared
    return prepared


def _topology_optic_disc_mask(
    optic_disc: OpticDisc,
    image_shape: tuple[int, int],
    settings: AnnulusGeometry,
) -> np.ndarray:
    fallback_radius = (
        float(settings.inner_radius_frac)
        * max(image_half_diagonal(*image_shape), 1.0)
    )
    return optic_disc.centered_circle_mask_for(
        image_shape,
        fallback_radius_pixels=fallback_radius,
    )


def _shared_window_side(
    topologies: Mapping[str, SegmentSamplingPlan],
    percentile_kept: float,
) -> int:
    percentile = float(percentile_kept)
    if not 0.0 < percentile <= 1.0:
        raise ValueError("window_size_percentile_kept must be in (0, 1].")

    widths: list[int] = []
    heights: list[int] = []
    for prepared in topologies.values():
        topology = prepared.native
        for annulus in topology.annulus_masks:
            for branch_id in topology.branch_ids:
                y, x = np.nonzero(annulus & (topology.labels == int(branch_id)))
                if x.size == 0:
                    continue
                widths.append(int(x.max() - x.min() + 1))
                heights.append(int(y.max() - y.min() + 1))

    if not widths:
        return 0

    width = int(np.quantile(widths, percentile, method="higher"))
    height = int(np.quantile(heights, percentile, method="higher"))
    side = max(width, height)
    return side if side % 2 == 1 else side + 1
