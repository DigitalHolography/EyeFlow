"""Selectable waveform velocity products."""

from pipeline_engine.imports import PipelineOption, pipeline

from .runner import run_waveform_velocity


@pipeline(
    name="waveform_velocity",
    description=(
        "Compute raw and filtered waveform velocity with optional derived products."
    ),
    requires=["numpy", "h5py", "scipy", "skimage"],
    dag_requires=["waveform_velocity_core"],
    dag_produces=["waveform_velocity"],
    options=[
        PipelineOption(
            "velocity_profiles",
            "Velocity profiles",
            "Per-beat transverse cross-section velocity profiles.",
        ),
        PipelineOption(
            "per_beat",
            "Per beat",
            "Raw and band-limited vessel velocity for each beat.",
        ),
        PipelineOption(
            "hemifield",
            "Hemifield",
            "Eight-region velocity and per-beat velocity aggregates.",
        ),
    ],
    input_slot="both",
)
def run(ctx):
    return run_waveform_velocity(ctx)


__all__ = ["run"]
