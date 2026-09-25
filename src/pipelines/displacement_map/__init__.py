"""Dense retinal displacement-map pipeline."""

from pipeline_engine.imports import pipeline

from .runner import run_displacement_map


@pipeline(
    name="displacement_map",
    description=(
        "Estimate dense retinal displacement from a HoloDoppler moment sequence "
        "inside a DopplerView vessel mask."
    ),
    requires=["numpy", "h5py", "cv2", "SimpleITK", "tqdm"],
    dag_produces=["displacement_map"],
    input_slot="both",
)
def run(ctx) -> None:
    """Entrypoint called by the pipeline runtime."""

    run_displacement_map(ctx)


__all__ = ["run"]
