"""Moments PDF and HTML report generation pipeline."""

from pipeline_engine import pipeline

from .runner import run_pdf_report, pdf_report_graph


@pipeline(
    name="pdf_report",
    description="Generate moments velocity/BVR PDF and HTML reports.",
    requires=["numpy", "matplotlib", "PIL"],
    dag_requires=[
        "velocity_analysis",
        "waveform_shape_metrics",
        "blood_volume_rate",
    ],
    graph_factory=pdf_report_graph,
    input_slot="both",
    wait_for_artifacts=True,
)
def run(ctx):
    """Entrypoint called by the pipeline runtime."""
    return run_pdf_report(ctx)


__all__ = ["run"]
