"""Optional topology-aligned spatial-gradient profile pipeline."""

from pipeline_engine import pipeline

from .runner import run_spatial_gradient_moment0, spatial_gradient_graph


@pipeline(
    name="spatial_gradient_moment0",
    description=(
        "Compute topology-aligned moment0 spatial-gradient profiles, lumen "
        "metrics, and reusable lumen-edge outputs."
    ),
    requires=["numpy", "h5py", "scipy", "skimage", "matplotlib"],
    dag_requires=["cardiac_cycles", "prepared_topology"],
    dag_produces=["spatial_gradient_moment0", "spatial_gradient_edges"],
    graph_factory=spatial_gradient_graph,
    input_slot="both",
)
def run(ctx):
    return run_spatial_gradient_moment0(ctx)


__all__ = ["run"]
