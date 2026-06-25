"""PDF report generation pipeline."""

from pipeline_engine.imports import pipeline

from .runner import run_pdf_report


@pipeline(
    name="pdf_report_generator",
    description="Generate A4 PDF report from EyeFlow analysis outputs.",
    requires=["numpy", "matplotlib", "PIL"],
    dag_requires=["waveform_shape_metrics"],
    input_slot="both",
)
def run(ctx):
    """Entrypoint called by the pipeline runtime."""
    return run_pdf_report(ctx)


__all__ = ["run"]