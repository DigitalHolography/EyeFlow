"""Selectable blood-volume-rate products."""

from pipeline_engine.imports import PipelineOption, pipeline

from .runner import run_blood_volume_rate


@pipeline(
    name="blood_volume_rate",
    description=(
        "Compute signed physical blood-volume rate from gradient- or mask-derived lumen geometry."
    ),
    requires=["numpy", "h5py", "scipy", "skimage"],
    options=[
        PipelineOption(
            "gradient_edges",
            "Spatial-gradient edges",
            "Dynamic- and static-edge blood-volume rate.",
            dag_requires=("spatial_gradient_edges", "velocity_profiles"),
        ),
        PipelineOption(
            "masked_edges",
            "Mask-derived geometry",
            "Masked-edge and total masked-edge blood-volume rate.",
            dag_requires=("segment_velocity_per_beat", "prepared_topology"),
        ),
    ],
    dag_produces=["blood_volume_rate"],
    input_slot="both",
    default_selected=True,
)
def run(ctx):
    return run_blood_volume_rate(ctx)


__all__ = ["run"]
