"""Selectable velocity analysis products."""

from pipeline_engine.imports import PipelineOption, ProcessResult, pipeline

from .builder import pack_velocity_analysis_meta_outputs, velocity_analysis
from .models import VelocityAnalysis
from .runner import run_velocity_analysis


@pipeline(
    name="velocity_analysis",
    description=(
        "Analyze continuous and per-beat velocity with optional spatial products."
    ),
    requires=["numpy", "h5py", "scipy", "skimage", "matplotlib"],
    dag_requires=["velocity", "prepared_topology"],
    dag_produces=["velocity_analysis"],
    options=[
        PipelineOption(
            "segments",
            "Segments",
            "Continuous and per-beat velocity signals for spatial vessel segments.",
        ),
        PipelineOption(
            "segment_velocity_maps",
            "Segment velocity maps",
            (
                "Per-beat segment velocity-map datasets and first-beat "
                "artery/vein mosaic movies."
            ),
            default_enabled=False,
            requires=("segments",),
        ),
        PipelineOption(
            "velocity_profiles",
            "Velocity profiles",
            "Per-beat cross-section velocity profiles.",
            default_enabled=False,
            requires=("segments",),
        ),
        PipelineOption(
            "velocity_profile_analysis",
            "Velocity profile analysis",
            "Weighted quadratic fits of artery and vein transverse profiles.",
            default_enabled=False,
            requires=("velocity_profiles",),
        ),
        PipelineOption(
            "velocity_profile_fft",
            "Velocity profile FFT",
            "Per-beat temporal FFTs of cross-section velocity profiles.",
            default_enabled=False,
            requires=("velocity_profiles",),
        ),
        PipelineOption(
            "quadrants",
            "Quadrants",
            "Four-quadrant velocity and per-beat velocity aggregates.",
            requires=("segments",),
        ),
    ],
    input_slot="both",
)
def run(ctx) -> ProcessResult:
    metrics = run_velocity_analysis(ctx)
    analysis = velocity_analysis(ctx)
    metrics.update(pack_velocity_analysis_meta_outputs(analysis))
    return ProcessResult(metrics=metrics, attrs=analysis.attrs)


__all__ = ["VelocityAnalysis", "run", "velocity_analysis"]
