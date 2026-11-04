"""Selectable blood-volume-rate products."""

from pipeline_engine import PipelineOption, pipeline

from .runner import run_blood_volume_rate, blood_volume_rate_graph


@pipeline(
    name="blood_volume_rate",
    description=(
        "Compute signed physical blood-volume rate from gradient- or mask-derived lumen geometry."
    ),
    requires=["numpy", "h5py", "scipy", "skimage", "matplotlib"],
    dag_requires=["velocity_analysis"],
    options=[
        PipelineOption(
            "gradient_edges",
            "Spatial-gradient edges",
            "Dynamic- and static-edge blood-volume rate.",
            default_enabled=False,
            dag_requires=("spatial_gradient_edges",),
        ),
        PipelineOption(
            "masked_edges",
            "Mask-derived geometry",
            "Masked-edge and total masked-edge blood-volume rate.",
        ),
    ],
    dag_produces=["blood_volume_rate"],
    graph_factory=blood_volume_rate_graph,
    input_slot="both",
    default_selected=True,
)
def run(ctx):
    return run_blood_volume_rate(ctx)


__all__ = ["run"]
