"""Selectable waveform-shape metric products."""

from pipeline_engine.imports import PipelineOption, pipeline

from .runner import run_waveform_shape_metrics


@pipeline(
    name="waveform_shape_metrics",
    description=(
        "Compute global, by-segment, and quadrant waveform-shape metrics."
    ),
    requires=["numpy", "h5py", "scipy", "skimage"],
    dag_requires=["waveform_velocity"],
    options=[
        PipelineOption(
            "per_beat",
            "Per beat",
            "Global waveform-shape metrics per beat.",
        ),
        PipelineOption(
            "segments",
            "Segments",
            "By-segment waveform-shape metrics per beat.",
            requires=("per_beat",),
        ),
        PipelineOption(
            "quadrants",
            "Quadrants",
            "Four-quadrant waveform-shape metric aggregates.",
            requires=("per_beat",),
        ),
    ],
    dag_produces=["waveform_shape_metrics"],
    input_slot="both",
)
def run(ctx):
    return run_waveform_shape_metrics(ctx)
