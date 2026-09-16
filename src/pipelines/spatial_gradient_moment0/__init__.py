"""Independent spatial-gradient artifacts and vessel profile pipeline."""

from pipeline_engine.imports import pipeline

from .runner import run_spatial_gradient_moment0


@pipeline(
    name="spatial_gradient_moment0",
    description=(
        "Apply a centered temporal median and export a contrast-enhanced "
        "framewise Sobel spatial gradient and topology-aligned vessel profiles."
    ),
    requires=["numpy", "scipy", "PIL", "h5py", "skimage", "cv2"],
    dag_requires=["heartbeat"],
    dag_produces=["spatial_gradient_moment0"],
    input_slot="both",
)
def run(ctx) -> None:
    run_spatial_gradient_moment0(ctx)


__all__ = ["run"]
