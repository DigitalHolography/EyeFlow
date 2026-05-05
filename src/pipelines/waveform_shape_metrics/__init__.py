"""Pipeline 1 MVP: compute DopplerView analysis, then AE waveform metrics."""

from __future__ import annotations

from pipeline_engine import pipeline

from .runner import run_waveform_shape_metrics


@pipeline(
    name="waveform_shape_metrics",
    description="Pipeline 1 MVP: global per-beat velocity outputs for AngioEye.",
    requires=["numpy", "h5py", "scipy", "skimage"],
    dag_produces=[
        "dopplerview_analysis",
        "velocity_per_beat",
        "waveform_shape_metrics",
    ],
    input_slot="both",
)
def run(ctx) -> None:
    metrics, attrs = run_waveform_shape_metrics(ctx)
    ctx.write_many(metrics)
    ctx.set_attrs(attrs)
