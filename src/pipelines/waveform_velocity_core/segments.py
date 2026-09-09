"""Prepare and analyze velocity segments for configured vessel classes."""

from __future__ import annotations

from collections.abc import Mapping, MutableMapping

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
        + ("CuPy/CUDA" if backend is not None else "CPU with parallel segments")
        + "."
    )
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

    results: dict[str, CrossSectionSignalResult] = {}
    for name, topology in topologies.items():
        segments = prepare_segments(velocity_map, topology)
        results[name] = _generate_cross_section_signals_from_prepared(
            velocity_map,
            topology,
            segments,
            ring_settings,
            cross_section_settings,
        )
    return results
