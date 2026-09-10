"""Weighted velocity-profile analysis pipeline."""

from pipeline_engine.imports import pipeline

from .runner import run_velocity_profile_analysis


@pipeline(
    name="velocity_profile_analysis",
    description="Analyze artery and vein velocity profiles with weighted quadratic fits.",
    requires=["numpy", "h5py", "scipy", "skimage"],
    dag_requires=["waveform_velocity"],
    dag_produces=["velocity_profile_analysis"],
    input_slot="both",
)
def run(ctx):
    return run_velocity_profile_analysis(ctx)


__all__ = ["run"]
