"""Hidden core pipeline for velocity estimation and cardiac-cycle timing."""

from pipeline_engine.imports import ProcessResult, pipeline

from .models import RetinalVelocity
from .runner import (
    cardiac_cycle_indexes,
    cardiac_cycles,
    run_velocity,
    velocity,
)


@pipeline(
    name="velocity",
    description=(
        "Compute base retinal velocity, average maps, and cardiac-cycle timing."
    ),
    requires=["numpy", "h5py", "scipy", "skimage"],
    dag_produces=["velocity", "cardiac_cycles"],
    input_slot="both",
    visibility="hidden",
)
def run(ctx) -> ProcessResult:
    velocity, metrics = run_velocity(ctx)
    return ProcessResult(
        metrics=metrics,
        attrs={
            "analysis_source": "eyeflow_velocity",
            "cardiac_cycle_detection_source": velocity.cardiac_cycle_source,
        },
    )


__all__ = [
    "RetinalVelocity",
    "cardiac_cycle_indexes",
    "cardiac_cycles",
    "run",
    "velocity",
]
