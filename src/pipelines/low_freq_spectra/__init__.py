"""Hidden shared low-frequency spectral analysis after velocity estimation."""

from pipeline_engine.imports import pipeline


@pipeline(
    name="low_freq_spectra",
    description="Compute spatial, gap, and beat-resolved spectra from optional HD spectra.",
    requires=["numpy", "h5py", "scipy", "matplotlib"],
    dag_requires=["velocity_workflows", "capillary_analysis_zones"],
    dag_produces=["gap_spectrograms", "gap_per_beat_spectrograms"],
    input_slot="hd",
    visibility="hidden",
    auto_run_after="velocity",
    run_once=True,
)
def run(ctx) -> None:
    from .runner import run_low_freq_spectra

    run_low_freq_spectra(ctx)
