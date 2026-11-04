"""Dense retinal displacement-map pipeline."""

from pipeline_engine.imports import PipelineOption, pipeline

from .runner import run_displacement_map, run_displacement_profiles


@pipeline(
    name="displacement_map",
    description=(
        "Estimate dense retinal displacement from a HoloDoppler moment sequence "
        "inside a DopplerView vessel mask."
    ),
    requires=["numpy", "h5py", "cv2", "SimpleITK", "scipy"],
    dag_produces=["displacement_map"],
    input_slot="both",
    execution_mode="process",
    native_parallelism=True,
    run_once=True,
)
def run(ctx):
    """Entrypoint called by the pipeline runtime."""

    return run_displacement_map(ctx)


@pipeline(
    name="displacement_profiles",
    description=(
        "Measure topology-aligned, beat-resolved displacement profiles and "
        "optional dense per-segment maps."
    ),
    requires=["numpy", "scipy", "cv2"],
    dag_requires=["displacement_map", "prepared_topology", "velocity_workflows"],
    dag_produces=["displacement_profiles"],
    options=[
        PipelineOption(
            "segment_maps",
            "Per-segment displacement maps",
            "Dense per-beat local displacement maps for every vessel segment.",
            default_enabled=False,
        ),
    ],
    input_slot="both",
    execution_mode="process",
    run_once=True,
)
def run_profiles(ctx):
    """Run the later analysis after fields, topology, and timing are ready."""

    return run_displacement_profiles(ctx)


__all__ = ["run", "run_profiles"]
