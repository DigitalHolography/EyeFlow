"""Gain-sensitive absolute waveform metric products."""

from pipeline_engine import PipelineOption, pipeline

from .runner import run_absolute_waveform_metrics, absolute_waveform_graph


@pipeline(
    name="absolute_waveform_metrics",
    description=(
        "Compute gain-sensitive absolute velocity, VTI, derivative, energy, "
        "harmonic, raw-vs-band-limited QC, segment, and quadrant metrics."
    ),
    requires=["numpy", "h5py"],
    dag_requires=["velocity_analysis"],
    dag_produces=["absolute_waveform_metrics"],
    options=[
        PipelineOption(
            "segments",
            "Segments",
            "Per-segment and segment-aggregate absolute metrics.",
        ),
        PipelineOption(
            "quadrants",
            "Quadrants",
            "Four-quadrant absolute metric aggregates.",
        ),
    ],
    graph_factory=absolute_waveform_graph,
    input_slot="both",
)
def run(ctx):
    return run_absolute_waveform_metrics(ctx)


__all__ = ["run_absolute_waveform_metrics"]
