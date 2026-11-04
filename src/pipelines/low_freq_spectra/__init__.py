"""Selectable, once-per-run low-frequency spectral analysis."""

from pipeline_engine import pipeline


from .runner import spectral_graph, run_low_freq_spectra


@pipeline(
    name="low_freq_spectra",
    description="Compute spatial, gap, and beat-resolved spectra from optional HD spectra.",
    requires=["numpy", "h5py", "scipy", "matplotlib"],
    dag_requires=["velocity_workflows", "capillary_analysis_zones"],
    dag_produces=["gap_spectrograms", "gap_per_beat_spectrograms"],
    input_slot="hd",
    run_once=True,
    graph_factory=spectral_graph,
)
def run(ctx):
    return run_low_freq_spectra(ctx)
