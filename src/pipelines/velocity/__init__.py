"""Hidden core pipeline for velocity estimation and cardiac-cycle timing."""

from pipeline_engine import pipeline

from .models import RetinalVelocity, VelocityWorkflows
from .runner import (
    cardiac_cycle_indexes,
    run_velocity,
    velocity_graph,
    velocity,
    velocity_workflows,
)


@pipeline(
    name="velocity",
    description=(
        "Compute base retinal velocity, average maps, and cardiac-cycle timing."
    ),
    requires=["numpy", "h5py", "scipy", "skimage", "cv2"],
    dag_produces=["velocity", "cardiac_cycles", "velocity_workflows"],
    graph_factory=velocity_graph,
    input_slot="both",
    visibility="hidden",
    produces_execution_variants=True,
)
def run(ctx) -> None:
    return run_velocity(ctx)


__all__ = [
    "RetinalVelocity",
    "VelocityWorkflows",
    "cardiac_cycle_indexes",
    "run",
    "velocity",
    "velocity_workflows",
]
