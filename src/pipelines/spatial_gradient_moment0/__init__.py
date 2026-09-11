"""Hidden spatial-gradient prerequisite for waveform velocity."""

from pipeline_engine.imports import pipeline

from .runner import run_spatial_gradient_moment0


@pipeline(
    name="spatial_gradient_moment0",
    description=(
        "Export a contrast-enhanced framewise Sobel spatial gradient of "
        "HoloDoppler moment0ff."
    ),
    requires=["numpy", "scipy", "PIL"],
    dag_produces=["spatial_gradient_moment0"],
    input_slot="hd",
    visibility="hidden",
)
def run(ctx) -> None:
    run_spatial_gradient_moment0(ctx)


__all__ = ["run"]
