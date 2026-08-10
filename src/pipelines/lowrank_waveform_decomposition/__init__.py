"""Low-rank SVD decomposition for beat-aligned segment waveforms."""

from pipeline_engine.imports import PipelineOption, pipeline

from .runner import run_lowrank_waveform_decomposition


@pipeline(
    name="lowrank_waveform_decomposition",
    description=(
        "Joint and per-beat low-rank SVD decomposition of beat-aligned "
        "arterial segment waveforms, with optional venous processing, "
        "reporting R0, modal/residual amplitudes, robust ratios, and "
        "singular-spectrum endpoints."
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
    ],
    input_slot="both",
)
def run(ctx):
    return run_lowrank_waveform_decomposition(ctx)


__all__ = ["run_lowrank_waveform_decomposition"]
