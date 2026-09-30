"""Hidden core pipeline for retinal velocity and cardiac-cycle timing."""

from pipeline_engine.imports import ProcessResult, pipeline

from .models import RetinalVelocity
from .runner import (
    cardiac_cycle_indexes,
    cardiac_cycles,
    retinal_velocity,
    run_retinal_velocity,
)


@pipeline(
    name="retinal_velocity",
    description=(
        "Compute base retinal velocity, average maps, and cardiac-cycle timing."
    ),
    requires=["numpy", "h5py", "scipy", "skimage"],
    dag_produces=["retinal_velocity", "cardiac_cycles"],
    input_slot="both",
    visibility="hidden",
)
def run(ctx) -> ProcessResult:
    velocity, metrics = run_retinal_velocity(ctx)
    return ProcessResult(
        metrics=metrics,
        attrs={
            "analysis_source": "eyeflow_retinal_velocity",
            "cardiac_cycle_detection_source": velocity.cardiac_cycle_source,
        },
    )


__all__ = [
    "RetinalVelocity",
    "cardiac_cycle_indexes",
    "cardiac_cycles",
    "retinal_velocity",
    "run",
]
