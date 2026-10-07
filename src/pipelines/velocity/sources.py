"""Resolve canonical source data required for retinal velocity."""

from __future__ import annotations

from input_output.schema import RetinalSourceData
from pipelines.vessel_inputs import load_vessel_topology_inputs


def load_velocity_inputs(ctx) -> RetinalSourceData:
    """Load and spatially align the retinal-velocity inputs."""

    return load_vessel_topology_inputs(ctx)


__all__ = ["load_velocity_inputs"]
