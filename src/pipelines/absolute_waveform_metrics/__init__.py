"""Gain-sensitive absolute waveform metric products."""

from pipeline_engine.imports import PipelineOption, pipeline

from .runner import run_absolute_waveform_metrics


@pipeline(
    name="absolute_waveform_metrics",
    description=(
        "Compute gain-sensitive absolute velocity, VTI, derivative, energy, "
        "harmonic, raw-vs-band-limited QC, segment, and hemifield metrics."
    ),
    requires=["numpy", "h5py"],
    dag_requires=["waveform_velocity"],
    dag_produces=["absolute_waveform_metrics"],
    options=[
        PipelineOption(
            "per_beat",
            "Per beat",
            "Global absolute waveform metrics for each beat.",
        ),
        PipelineOption(
            "segments",
            "Segments",
            "Per-segment and segment-aggregate absolute metrics.",
            requires=("per_beat",),
        ),
        PipelineOption(
            "hemifield",
            "Hemifield",
            "Eight-region absolute metric aggregates.",
            requires=("per_beat",),
        ),
    ],
    input_slot="both",
)
def run(ctx):
    return run_absolute_waveform_metrics(ctx)


__all__ = ["run_absolute_waveform_metrics"]
