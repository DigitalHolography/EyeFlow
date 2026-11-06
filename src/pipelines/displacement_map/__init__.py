"""Dense retinal displacement-map pipeline."""

from pipeline_engine import PipelineOption, pipeline
from pipeline_engine.product_exchange import ProductRequirement

from .runner import run_displacement_map, run_displacement_profiles
from .runner import displacement_graph, displacement_profiles_graph


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
    native_parallelism=False,
    run_once=True,
    graph_factory=displacement_graph,
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
    dag_requires=["prepared_topology", "velocity_workflows"],
    product_requires=[ProductRequirement(alias="registration", producer_pipeline="displacement_map",
                                         product="registered_fields", variant="shared")],
    graph_factory=displacement_profiles_graph,
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
