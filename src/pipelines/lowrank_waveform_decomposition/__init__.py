"""Low-rank SVD decomposition for beat-aligned segment waveforms."""

from pipeline_engine.imports import PipelineOption, pipeline

from .runner import run_lowrank_waveform_decomposition


@pipeline(
    name="lowrank_waveform_decomposition",
    description=(
        "Joint and per-beat low-rank SVD decomposition of beat-aligned "
        "arterial and venous segment waveforms."
    ),
    requires=["numpy", "h5py"],
    dag_requires=["waveform_velocity"],
    dag_produces=["lowrank_waveform_decomposition"],
    options=[
        PipelineOption(
            "quadrants",
            "Quadrants",
            "Four-quadrant low-rank waveform decomposition metrics.",
        ),
    ],
    input_slot="both",
)
def run(ctx):
    return run_lowrank_waveform_decomposition(ctx)


__all__ = ["run_lowrank_waveform_decomposition"]
