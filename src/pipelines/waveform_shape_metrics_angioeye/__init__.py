"""Registered AngioEye waveform-shape pipeline for EyeFlow."""

from pipeline_engine.imports import ProcessResult, pipeline

from .runner import run_waveform_shape_metrics_angioeye


@pipeline(
    name="waveform_shape_metrics_angioeye",
    description="AngioEye waveform-shape metrics over EyeFlow velocity outputs.",
    requires=["numpy"],
    dag_requires=["waveform_shape_metrics"],
    input_slot="both",
)
def run(ctx) -> ProcessResult:
    """Run after EyeFlow produces its per-beat waveform datasets."""

    return ProcessResult(metrics=run_waveform_shape_metrics_angioeye(ctx))


__all__ = ["run"]
