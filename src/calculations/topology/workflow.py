"""Logical execution order for reusable retinal segment topology."""

from __future__ import annotations

from collections.abc import Mapping, MutableMapping
from dataclasses import dataclass
from time import perf_counter

import numpy as np

from utils.logger import Logger

from .cache import TopologyCacheKey, topology_cache_key
from .geometry import SegmentRingSettings
from .segments import SegmentTopology, build_segment_topology, extract_segments
from .transforms import (
    determine_segment_rotations,
    interpolate_segment_masks,
    interpolate_segments,
    rotate_segment_masks,
    rotate_segments,
)


@dataclass(frozen=True)
class PreparedTopology:
    """Map-independent segment geometry ready for map transformations."""

    topology: SegmentTopology
    rotation_degrees: np.ndarray
    interpolated_masks: np.ndarray
    rotated_masks: np.ndarray


@dataclass(frozen=True)
class PreparedSegments:
    """Uniform non-rotated and upright views extracted from one retinal map."""

    interpolated: np.ndarray
    rotated: np.ndarray


def prepare_topology(
    vessel_mask,
    optic_disc_mask,
    settings: SegmentRingSettings,
    *,
    output_side_pixels: int = 128,
    window_size_percentile_kept: float = 0.95,
    window_side_pixels: int | None = None,
) -> PreparedTopology:
    """Build segment geometry, orientations, and uniform masks once."""

    topology = build_segment_topology(
        vessel_mask,
        optic_disc_mask,
        settings,
        window_size_percentile_kept=window_size_percentile_kept,
        window_side_pixels=window_side_pixels,
    )
    rotation_degrees = determine_segment_rotations(topology)
    interpolated_masks = interpolate_segment_masks(
        topology.segment_masks,
        output_side_pixels,
    )
    return PreparedTopology(
        topology=topology,
        rotation_degrees=rotation_degrees,
        interpolated_masks=interpolated_masks,
        rotated_masks=rotate_segment_masks(
            interpolated_masks,
            rotation_degrees,
        ),
    )


def prepare_topologies(
    vessel_masks: Mapping[str, object],
    optic_disc_mask,
    settings: SegmentRingSettings,
    *,
    source_id: str,
    cache: MutableMapping[TopologyCacheKey, object] | None = None,
    optic_disc_center=None,
    output_side_pixels: int = 128,
    window_size_percentile_kept: float = 0.95,
    window_side_pixels: int | None = None,
) -> dict[str, PreparedTopology]:
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
    disc = _resolved_optic_disc_mask(
        optic_disc_mask,
        optic_disc_center,
        image_shape,
    )

    if window_side_pixels is not None:
        return {
            name: _cached_topology(
                name,
                mask,
                disc,
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
            disc,
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
    return {
        name: (
            topology
            if topology.topology.window_side_pixels == shared_side
            else _cached_topology(
                name,
                masks[name],
                disc,
                settings,
                source_id=source_id,
                cache=cache,
                output_side_pixels=output_side_pixels,
                window_size_percentile_kept=window_size_percentile_kept,
                window_side_pixels=shared_side,
            )
        )
        for name, topology in initial.items()
    }


def prepare_segments(
    data_map,
    prepared_topology: PreparedTopology,
    *,
    spatial_axes: tuple[int, int] = (-2, -1),
) -> PreparedSegments:
    """Extract, uniformly interpolate, and rotate one map's segment arrays."""

    total_started = perf_counter()
    Logger.log(
        "Starting topology segment preparation: "
        f"source_shape={tuple(int(size) for size in data_map.shape)}, "
        f"source_type={type(data_map).__name__}."
    )
    started = perf_counter()
    extracted = extract_segments(
        data_map,
        prepared_topology.topology,
        spatial_axes=spatial_axes,
    )
    Logger.log(
        "Completed segment extraction in "
        f"{perf_counter() - started:.2f}s; {_array_summary(extracted)}."
    )
    started = perf_counter()
    interpolated = interpolate_segments(
        extracted,
        prepared_topology.interpolated_masks.shape[-1],
    )
    Logger.log(
        "Completed segment interpolation in "
        f"{perf_counter() - started:.2f}s; {_array_summary(interpolated)}."
    )
    started = perf_counter()
    rotated = rotate_segments(
        interpolated,
        prepared_topology.rotation_degrees,
    )
    Logger.log(
        "Completed segment rotation in "
        f"{perf_counter() - started:.2f}s; {_array_summary(rotated)}."
    )
    prepared = PreparedSegments(
        interpolated=interpolated,
        rotated=rotated,
    )
    Logger.log(
        f"Completed topology segment preparation in {perf_counter() - total_started:.2f}s."
    )
    return prepared


def _cached_topology(
    vessel_name: str,
    vessel_mask: np.ndarray,
    optic_disc_mask: np.ndarray,
    settings: SegmentRingSettings,
    *,
    source_id: str,
    cache: MutableMapping[TopologyCacheKey, object] | None,
    output_side_pixels: int,
    window_size_percentile_kept: float,
    window_side_pixels: int | None,
) -> PreparedTopology:
    key = topology_cache_key(
        source_id,
        vessel_name,
        vessel_mask,
        optic_disc_mask,
        settings,
        output_side_pixels=output_side_pixels,
        window_size_percentile_kept=window_size_percentile_kept,
        window_side_pixels=window_side_pixels,
    )
    if cache is not None:
        found = cache.get(key)
        if found is not None:
            if not isinstance(found, PreparedTopology):
                raise TypeError("Topology cache values must be PreparedTopology.")
            return found

    prepared = prepare_topology(
        vessel_mask,
        optic_disc_mask,
        settings,
        output_side_pixels=output_side_pixels,
        window_size_percentile_kept=window_size_percentile_kept,
        window_side_pixels=window_side_pixels,
    )
    if cache is not None:
        cache[key] = prepared
    return prepared


def _resolved_optic_disc_mask(
    optic_disc_mask,
    optic_disc_center,
    image_shape: tuple[int, int],
) -> np.ndarray:
    if optic_disc_mask is not None:
        disc = np.asarray(optic_disc_mask, dtype=bool)
        if disc.shape != image_shape:
            raise ValueError(
                f"optic_disc_mask must have shape {image_shape}, got {disc.shape}."
            )
        return disc

    disc = np.zeros(image_shape, dtype=bool)
    if optic_disc_center is None:
        y = image_shape[0] // 2
        x = image_shape[1] // 2
    else:
        center = np.asarray(optic_disc_center, dtype=np.float32).reshape(-1)
        if center.size < 2 or not np.all(np.isfinite(center[:2])):
            y = image_shape[0] // 2
            x = image_shape[1] // 2
        else:
            x = int(np.clip(np.floor(center[0] + 0.5), 0, image_shape[1] - 1))
            y = int(np.clip(np.floor(center[1] + 0.5), 0, image_shape[0] - 1))
    disc[y, x] = True
    return disc


def _shared_window_side(
    topologies: Mapping[str, PreparedTopology],
    percentile_kept: float,
) -> int:
    percentile = float(percentile_kept)
    if not 0.0 < percentile <= 1.0:
        raise ValueError("window_size_percentile_kept must be in (0, 1].")

    widths: list[int] = []
    heights: list[int] = []
    for prepared in topologies.values():
        topology = prepared.topology
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


def _array_summary(values: np.ndarray) -> str:
    array = np.asarray(values)
    return (
        f"shape={array.shape}, dtype={array.dtype}, "
        f"allocated={array.nbytes / (1024 ** 3):.2f} GiB"
    )
