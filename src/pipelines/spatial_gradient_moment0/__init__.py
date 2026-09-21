"""Optional topology-aligned spatial-gradient profile pipeline."""

from pipeline_engine.imports import pipeline

from .runner import run_spatial_gradient_moment0


@pipeline(
    name="spatial_gradient_moment0",
    description=(
        "Compute topology-aligned moment0 spatial-gradient profiles, lumen "
        "metrics, and gradient-derived blood-volume-rate outputs."
    ),
    requires=["numpy", "h5py", "scipy", "skimage", "matplotlib"],
    dag_requires=["waveform_velocity_core"],
    dag_produces=["spatial_gradient_moment0"],
    input_slot="both",
)
def run(ctx):
    return run_spatial_gradient_moment0(ctx)


__all__ = ["run"]
