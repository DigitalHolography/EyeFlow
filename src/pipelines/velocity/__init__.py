"""Hidden core pipeline for velocity estimation and cardiac-cycle timing."""

from pipeline_engine.imports import pipeline

from .models import RetinalVelocity, VelocityWorkflows
from .runner import (
    cardiac_cycle_indexes,
    cardiac_cycles,
    run_velocity,
    velocity,
    velocity_workflows,
)


@pipeline(
    name="velocity",
    description=(
        "Compute base retinal velocity, average maps, and cardiac-cycle timing."
    ),
    requires=["numpy", "h5py", "scipy", "skimage"],
    dag_produces=["velocity", "cardiac_cycles", "velocity_workflows"],
    input_slot="both",
    visibility="hidden",
    produces_execution_variants=True,
)
def run(ctx) -> None:
    run_velocity(ctx)
    # Each workflow writes through its own output namespace.


__all__ = [
    "RetinalVelocity",
    "VelocityWorkflows",
    "cardiac_cycle_indexes",
    "cardiac_cycles",
    "run",
    "velocity",
    "velocity_workflows",
]
