"""Run-scoped cache identity for prepared retinal topology."""

from __future__ import annotations

from collections.abc import MutableMapping
from dataclasses import dataclass
from hashlib import blake2b

import numpy as np

from .geometry import SegmentRingSettings


TOPOLOGY_CACHE_STATE = "topology.prepared"


@dataclass(frozen=True, slots=True)
class TopologyCacheKey:
    """Inputs that completely identify one prepared vessel topology."""

    source_id: str
    vessel_name: str
    vessel_mask_fingerprint: str
    competing_vessel_mask_fingerprint: str
    optic_disc_mask_fingerprint: str
    settings: SegmentRingSettings
    output_side_pixels: int
    window_size_percentile_kept: float
    window_side_pixels: int | None


def topology_cache_key(
    source_id: str,
    vessel_name: str,
    vessel_mask,
    optic_disc_mask,
    settings: SegmentRingSettings,
    *,
    output_side_pixels: int,
    window_size_percentile_kept: float,
    window_side_pixels: int | None,
    competing_vessel_mask=None,
) -> TopologyCacheKey:
    """Create a stable cache key from topology-defining inputs."""

    return TopologyCacheKey(
        source_id=str(source_id),
        vessel_name=str(vessel_name),
        vessel_mask_fingerprint=_mask_fingerprint(vessel_mask),
        competing_vessel_mask_fingerprint=_mask_fingerprint(
            np.zeros_like(np.asarray(vessel_mask, dtype=bool))
            if competing_vessel_mask is None
            else competing_vessel_mask
        ),
        optic_disc_mask_fingerprint=_mask_fingerprint(optic_disc_mask),
        settings=settings,
        output_side_pixels=int(output_side_pixels),
        window_size_percentile_kept=float(window_size_percentile_kept),
        window_side_pixels=(
            None if window_side_pixels is None else int(window_side_pixels)
        ),
    )


def run_topology_cache(
    run_state: MutableMapping[str, object],
) -> dict[TopologyCacheKey, object]:
    """Return the prepared-topology dictionary owned by one pipeline run."""

    found = run_state.get(TOPOLOGY_CACHE_STATE)
    if found is None:
        cache: dict[TopologyCacheKey, object] = {}
        run_state[TOPOLOGY_CACHE_STATE] = cache
        return cache
    if not isinstance(found, dict):
        raise TypeError(
            f"Run state '{TOPOLOGY_CACHE_STATE}' must contain a dictionary."
        )
    return found


def topology_source_id(*source_names: str | None) -> str:
    """Identify the external inputs whose aligned segmentation is being used."""

    return "|".join(str(name or "") for name in source_names)


def _mask_fingerprint(mask) -> str:
    values = np.ascontiguousarray(np.asarray(mask, dtype=bool))
    digest = blake2b(digest_size=16)
    digest.update(np.asarray(values.shape, dtype=np.int64).tobytes())
    digest.update(values.view(np.uint8))
    return digest.hexdigest()
