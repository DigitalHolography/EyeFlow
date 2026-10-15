"""Own once-per-run spectral reductions and cardiac-timing selection."""

from __future__ import annotations

from pipelines.topology_core.runner import optional_capillary_zones
from pipelines.velocity import velocity_workflows
from utils.logger import Logger

from .gap_analysis import run_gap_analysis
from .spectrograms import preprocess_spectrograms


def run_low_freq_spectra(ctx) -> None:
    ctx.require_inputs("hd")
    spectra = ctx.inputs.hd.as_holodoppler().optional_spectograms_dataset()
    if spectra is None:
        Logger.log("Skipping low_freq_spectra: HD /spectrograms video is unavailable.")
        return

    outputs = preprocess_spectrograms(spectra, optional_capillary_zones(ctx))
    ctx.output.h5.write_many(outputs)
    workflows = velocity_workflows(ctx)
    timing = None
    method = "doppler_moments"
    for method, candidate in (
        ("doppler_moments", workflows.moments),
        ("frequency_bands", workflows.bandratio),
    ):
        if candidate is not None and (
            candidate.cardiac_cycle_source != "none"
            and candidate.cycle_boundary_indexes.size >= 2
        ):
            timing = candidate
            break
    if timing is not None and method == "frequency_bands":
        Logger.log_warning(
            "Low-frequency spectra use band-ratio cardiac cycles: "
            "detected Doppler moments timing is unavailable."
        )
    ctx.output.h5.write_many(run_gap_analysis(ctx, timing, timing_method=method))
