"""Selectable waveform-shape metric products."""

from pipeline_engine.imports import PipelineOption, pipeline

from .runner import run_waveform_shape_metrics


@pipeline(
    name="waveform_shape_metrics",
    description=(
        "Compute global, by-segment, and hemifield waveform-shape metrics."
    ),
    requires=["numpy", "h5py", "scipy", "skimage"],
    dag_requires=["waveform_velocity_core"],
    options=[
        PipelineOption(
            "per_beat",
            "Per beat",
            "Global and segment waveform-shape metrics per beat.",
        ),
        PipelineOption(
            "hemifield",
            "Hemifield",
            "Eight-region waveform-shape metric aggregates.",
            requires=("per_beat",),
        ),
    ],
    dag_produces=["waveform_shape_metrics"],
    input_slot="both",
)
def run(ctx):
    return run_waveform_shape_metrics(ctx)
