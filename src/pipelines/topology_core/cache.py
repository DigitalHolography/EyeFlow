"""Run-state ownership for shared spatial sampling-plan caching."""

from __future__ import annotations

from collections.abc import MutableMapping

from calculations.vessel_segments.sampling.cache import TopologyCacheKey

# Keep the existing run-state key while moving its owner out of calculations.
TOPOLOGY_CACHE_STATE = "topology.prepared"


def run_topology_cache(
    run_state: MutableMapping[str, object],
) -> dict[TopologyCacheKey, object]:
    """Return the sampling-plan dictionary owned by one pipeline run."""

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
