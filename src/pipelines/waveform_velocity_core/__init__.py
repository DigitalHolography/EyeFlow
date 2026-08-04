"""Hidden shared foundation for waveform velocity pipelines."""

from pipeline_engine.imports import ProcessResult, pipeline

from .runner import run_waveform_velocity_core


@pipeline(
    name="waveform_velocity_core",
    description=(
        "Load waveform sources and run shared DopplerView with optional segment extraction."
    ),
    requires=["numpy", "h5py", "scipy", "skimage"],
    dag_produces=["dopplerview_analysis", "waveform_velocity_core"],
    input_slot="both",
    visibility="hidden",
)
def run(ctx) -> ProcessResult:
    metrics, attrs = run_waveform_velocity_core(ctx)
    return ProcessResult(metrics=metrics, attrs=attrs)


__all__ = ["run"]
