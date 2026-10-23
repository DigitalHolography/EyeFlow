"""Own once-per-run spectral reductions and cardiac-timing selection."""

from __future__ import annotations

from input_output.schema import EyeFlowOutputPaths
from pipelines.topology_core.runner import optional_capillary_zones
from pipelines.velocity import velocity_workflows
from utils.logger import Logger

from .artifacts import export_spectrogram_figures
from .gap_analysis import run_gap_analysis
from .frequency import write_frequency_bin_table
from .qc import pack_gap_qc
from .spectrograms import preprocess_spectrograms


def run_low_freq_spectra(ctx) -> None:
    ctx.require_inputs("hd")
    hd = ctx.inputs.hd.as_holodoppler()
    spectra = hd.optional_spectograms_dataset()
    if spectra is None:
        Logger.log("Skipping low_freq_spectra: HD /spectrograms video is unavailable.")
        return

    schema = EyeFlowOutputPaths.active()
    if ctx.output.h5.get(schema.sampling_frequency_hz) is None:
        try:
            sampling_frequency = hd.sampling_frequency_hz()
        except (KeyError, TypeError, ValueError) as exc:
            Logger.log_warning(f"Cannot calibrate spectral figure frequencies: {exc}")
        else:
            ctx.output.h5.write(schema.sampling_frequency_hz, sampling_frequency, unit="Hz")
    outputs = preprocess_spectrograms(spectra, optional_capillary_zones(ctx))
    ctx.output.h5.write_many(outputs)
    write_frequency_bin_table(ctx.output, spectra.shape[1], hd=hd)
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
    gap_outputs = run_gap_analysis(ctx, timing, timing_method=method)
    ctx.output.h5.write_many(gap_outputs)
    gap_source = ctx.output.h5.get(schema.spectrograms_masked)
    if gap_source is not None:
        ctx.output.h5.write_many(pack_gap_qc(
            gap_source[()], source_attrs=dict(gap_source.attrs), gap_outputs=gap_outputs,
        ))
    if ctx.output.available:
        try:
            dt_seconds = hd.timing().dt_seconds
        except (KeyError, TypeError, ValueError, ZeroDivisionError) as exc:
            Logger.log_warning(f"Cannot calibrate spectrogram time in seconds: {exc}")
            dt_seconds = None
        export_spectrogram_figures(ctx.output, dt_seconds=dt_seconds)
