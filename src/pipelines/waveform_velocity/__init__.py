"""Selectable waveform velocity products."""

from pipeline_engine.imports import PipelineOption, ProcessResult, pipeline

from .builder import pack_waveform_meta_outputs, waveform_velocity
from .models import WaveformVelocity
from .runner import run_waveform_velocity


@pipeline(
    name="waveform_velocity",
    description=(
        "Compute raw and band-limited waveform velocity with optional derived products."
    ),
    requires=["numpy", "h5py", "scipy", "skimage", "matplotlib"],
    dag_requires=["retinal_velocity", "prepared_topology"],
    dag_produces=["waveform_velocity"],
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
    metrics = run_waveform_velocity(ctx)
    waveform = waveform_velocity(ctx)
    metrics.update(pack_waveform_meta_outputs(waveform))
    return ProcessResult(metrics=metrics, attrs=waveform.attrs)


__all__ = ["WaveformVelocity", "run", "waveform_velocity"]
