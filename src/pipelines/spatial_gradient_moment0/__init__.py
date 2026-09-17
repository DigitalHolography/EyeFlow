"""Hidden spatial-gradient prerequisite for waveform velocity."""

from pipeline_engine.imports import pipeline

from .runner import run_spatial_gradient_moment0


@pipeline(
    name="spatial_gradient_moment0",
    description=(
        "Reuse velocity segment geometry, extract and rotate M0ff segments, "
        "then apply selectable temporal filters before and after Sobel gradients."
    ),
    requires=["numpy", "scipy", "PIL"],
    dag_requires=["waveform_velocity_core"],
    dag_produces=["spatial_gradient_moment0"],
    input_slot="hd",
    visibility="hidden",
)
def run(ctx) -> None:
    run_spatial_gradient_moment0(ctx)


__all__ = ["run"]
