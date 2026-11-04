"""Build the executable, once-per-run spectral analysis."""


def run_low_freq_spectra(ctx):
    from pipeline_engine.context import require_managed_runtime
    require_managed_runtime(ctx, "low_freq_spectra")
    from .graph import spectral_graph

    return spectral_graph()
