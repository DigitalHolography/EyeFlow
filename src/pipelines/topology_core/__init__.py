"""Hidden shared retinal topology preparation pipeline."""

from pipeline_engine.imports import pipeline

from .runner import run_topology_core


@pipeline(
    name="topology_core",
    description=(
        "Prepare canonical retinal vessel topology and capillary analysis zones."
    ),
    requires=["numpy", "h5py", "scipy", "skimage"],
    dag_produces=["prepared_topology", "capillary_analysis_zones"],
    input_slot="both",
    visibility="hidden",
)
def run(ctx):
    return run_topology_core(ctx)


__all__ = ["run"]
