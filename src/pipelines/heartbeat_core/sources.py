"""Resolve canonical source data required to detect heartbeat boundaries."""

from __future__ import annotations

from input_output.schema import RetinalSourceData
from pipelines.vessel_inputs import load_vessel_topology_inputs


def load_heartbeat_inputs(ctx) -> RetinalSourceData:
    """Load and spatially align only the arrays used for beat detection."""

    return load_vessel_topology_inputs(ctx)


__all__ = ["load_heartbeat_inputs"]
