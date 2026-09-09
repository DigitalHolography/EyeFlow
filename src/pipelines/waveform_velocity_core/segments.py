"""Prepare and analyze velocity segments for configured vessel classes."""

from __future__ import annotations

from collections.abc import Mapping, MutableMapping
from time import perf_counter

import numpy as np

from calculations.compute_backend import optional_cupy_backend
from calculations.topology import (
    SegmentRingSettings,
    TopologyCacheKey,
    prepare_segments,
    prepare_topologies,
)
from calculations.blood_flow_velocity.cross_section.generate_cross_section_signals import (
    CrossSectionSignalResult,
    CrossSectionSignalSettings,
    _generate_cross_section_signals_from_prepared,
    _validate_velocity_map,
)
from utils.logger import Logger


def analyze_velocity_segments(
    velocity_map,
    vessel_masks: Mapping[str, object],
    optic_disc_center,
    ring_settings: SegmentRingSettings,
    cross_section_settings: CrossSectionSignalSettings,
    *,
    optic_disc_mask=None,
    source_id: str = "",
    topology_cache: MutableMapping[TopologyCacheKey, object] | None = None,
) -> dict[str, CrossSectionSignalResult]:
    """Analyze velocity-map segments for every named vessel mask.

    Vessel names are retained as result keys. All selected vessels contribute
    to one shared segment-window size before their maps are prepared.
    """

    masks = {
        str(name): np.asarray(mask, dtype=bool)
        for name, mask in vessel_masks.items()
    }
    if not masks:
        return {}
    for mask in masks.values():
        _validate_velocity_map(velocity_map, mask)

    backend = optional_cupy_backend()
    Logger.log(
        "Cross-section compute backend: "
        + ("CuPy/CUDA" if backend is not None else "CPU/SciPy")
        + "."
    )
    Logger.log(
        "Segment input: "
        f"map_shape={tuple(int(size) for size in velocity_map.shape)}, "
        f"map_type={type(velocity_map).__name__}, vessels={tuple(masks)}."
    )
    topology_started = perf_counter()
    topologies = prepare_topologies(
        masks,
        optic_disc_mask,
        ring_settings,
        source_id=source_id,
        cache=topology_cache,
        optic_disc_center=optic_disc_center,
        window_size_percentile_kept=(
            cross_section_settings.submask_size_percentile_kept
        ),
    )
    Logger.log(
        f"Completed topology preparation in {perf_counter() - topology_started:.2f}s."
    )

    results: dict[str, CrossSectionSignalResult] = {}
    for name, topology in topologies.items():
        geometry = topology.topology
        Logger.log(
            f"Preparing {name} segments: radii={geometry.annulus_masks.shape[0]}, "
            f"branches={geometry.branch_ids.size}, "
            f"valid_segments={int(np.count_nonzero(geometry.valid_segments))}, "
            f"native_window={geometry.window_side_pixels}px."
        )
        preparation_started = perf_counter()
        segments = prepare_segments(velocity_map, topology)
        Logger.log(
            f"Prepared {name} segment arrays in "
            f"{perf_counter() - preparation_started:.2f}s."
        )
        measurement_started = perf_counter()
        results[name] = _generate_cross_section_signals_from_prepared(
            velocity_map,
            topology,
            segments,
            ring_settings,
            cross_section_settings,
        )
        Logger.log(
            f"Completed {name} segment profile measurements in "
            f"{perf_counter() - measurement_started:.2f}s."
        )
    return results
