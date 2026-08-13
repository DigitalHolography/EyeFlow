"""Low-rank SVD decomposition for beat-aligned segment waveforms."""

from pipeline_engine.imports import PipelineOption, pipeline

from .runner import run_lowrank_waveform_decomposition


@pipeline(
    name="lowrank_waveform_decomposition",
    description=(
        "Low-rank SVD decomposition of beat-aligned arterial segment "
        "waveforms, with optional venous processing. Joint SVD and "
        "per-beat SVD are independently selectable."
    ),
    requires=["numpy", "h5py"],
    dag_requires=["waveform_velocity"],
    dag_produces=["lowrank_waveform_decomposition"],
    options=[
        PipelineOption(
            "veins",
            "Veins",
            "Also decompose raw and band-limited venous segment waveforms.",
            default_enabled=False,
        ),
        PipelineOption(
            "joint_svd",
            "Joint SVD",
            "One SVD over all valid beat-location waveforms (primary endpoints).",
            default_enabled=True,
        ),
        PipelineOption(
            "per_beat_svd",
            "Per-beat SVD",
            "Independent SVD of each beat's spatial columns (robustness variant).",
            default_enabled=False,
        ),
    ],
    input_slot="both",
)
def run(ctx):
    return run_lowrank_waveform_decomposition(ctx)


__all__ = ["run_lowrank_waveform_decomposition"]
