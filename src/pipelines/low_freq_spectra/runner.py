"""Build the executable, once-per-run spectral analysis."""


def run_low_freq_spectra(ctx):
    from .graph import spectral_graph

    return spectral_graph()
