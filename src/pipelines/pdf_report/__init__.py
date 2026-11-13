"""PDF report generation pipeline."""

from pipeline_engine.imports import pipeline

from .runner import run_pdf_report


@pipeline(
    name="pdf_report",
    description="Generate moments velocity/BVR and shared spectra/gaps PDF reports.",
    requires=["numpy", "matplotlib", "PIL"],
    dag_requires=[
        "velocity_analysis",
        "waveform_shape_metrics",
        "blood_volume_rate",
    ],
    input_slot="both",
    wait_for_artifacts=True,
)
def run(ctx):
    """Entrypoint called by the pipeline runtime."""
    return run_pdf_report(ctx)


@pipeline(
    name="pdf_width_report",
    description="Summarize available displacement and spatial-gradient results in one PDF.",
    requires=["numpy", "matplotlib"],
    visibility="hidden",
    # Both displacement profiles and spatial gradients require velocity timing.
    # Waiting for artifacts keeps this shared report after their optional work.
    auto_run_after="velocity",
    run_once=True,
    wait_for_artifacts=True,
)
def run_width(ctx):
    """Write a shared third PDF without scheduling either optional analysis."""
    from .width_runner import run_width_pdf_report

    return run_width_pdf_report(ctx)


__all__ = ["run", "run_width"]
