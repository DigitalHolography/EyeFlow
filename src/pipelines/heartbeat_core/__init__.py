"""Hidden shared heartbeat-boundary detection pipeline."""

from pipeline_engine.imports import pipeline

from .runner import run_heartbeat_core


@pipeline(
    name="heartbeat_core",
    description=(
        "Detect heartbeat-cycle boundaries once for independent downstream pipelines."
    ),
    requires=["numpy", "h5py", "scipy", "skimage"],
    dag_produces=["heartbeat"],
    input_slot="both",
    visibility="hidden",
)
def run(ctx) -> None:
    run_heartbeat_core(ctx)


__all__ = ["run"]
