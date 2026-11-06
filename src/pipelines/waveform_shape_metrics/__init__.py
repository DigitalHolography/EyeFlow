"""Selectable waveform-shape metric products."""

from pipeline_engine import PipelineOption, pipeline

from .runner import run_waveform_shape_metrics, waveform_shape_graph


@pipeline(
    name="waveform_shape_metrics",
    description=(
        "Compute global, by-segment, and quadrant waveform-shape metrics."
    ),
    requires=["numpy", "h5py", "scipy", "skimage"],
    dag_requires=["velocity_analysis"],
    options=[
        PipelineOption(
            "segments",
            "Segments",
            "By-segment waveform-shape metrics per beat.",
        ),
        PipelineOption(
            "quadrants",
            "Quadrants",
            "Four-quadrant waveform-shape metric aggregates.",
        ),
    ],
    dag_produces=["waveform_shape_metrics"],
    graph_factory=waveform_shape_graph,
    input_slot="both",
)
def run(ctx):
    return run_waveform_shape_metrics(ctx)
