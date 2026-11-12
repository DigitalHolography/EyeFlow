"""Executable spectral stages and once-per-run orchestration."""

from functools import partial

from input_output.schema import EyeFlowOutputPaths
from pipeline_engine.stage_graph import StageGraph, StageSpec, ProductSpec
from pipeline_engine.stage_graph import InputAccess, InputUse
from pipelines.topology_core.runner import optional_capillary_zones
from pipelines.velocity import velocity_workflows
from utils.logger import Logger

from .artifacts import export_spectrogram_figures
from .figure_export import export_spectral_snapshot, snapshot_spectral_figures
from .frequency import pack_frequency_bin_table
from .gap_analysis import run_gap_analysis
from .qc import pack_gap_qc
from .spectrograms import preprocess_spectrograms


def _packed(value, ctx):
    return value


def _source(ctx):
    ctx.require_inputs("hd")
    try:
        spectra = ctx.inputs.hd.as_holodoppler().spectrograms()
    except (TypeError, ValueError) as exc:
        Logger.log_error(
            "Low-frequency spectra rejected the HD /spectrograms input: "
            f"{exc}"
        )
        raise
    if spectra is None:
        Logger.log_warning(
            "Skipping low-frequency spectra: HD /spectrograms is missing or "
            "uses the unsupported legacy-group layout."
        )
        ctx.dag.skip("HD /spectrograms video is unavailable.")
    ctx.dag.handoff("source", {"bin_count": spectra.shape[1]})


def _sampling(ctx):
    schema = EyeFlowOutputPaths.active()
    values = {}
    if ctx.output.h5.get(schema.sampling_frequency_hz) is None:
        try:
            fs = ctx.inputs.hd.as_holodoppler().sampling_frequency_hz()
        except (KeyError, TypeError, ValueError) as exc:
            Logger.log_warning(f"Cannot calibrate spectral figure frequencies: {exc}")
        else:
            values[schema.sampling_frequency_hz] = (fs, {"unit": "Hz"})
    ctx.dag.handoff("sampling", values)


def _calibration(ctx):
    ctx.dag.handoff("calibration", pack_frequency_bin_table(
        ctx.output, ctx.dag.require("source")["bin_count"], hd=ctx.inputs.hd.as_holodoppler(),
    ))


def _spatial(ctx):
    spectra = ctx.inputs.hd.as_holodoppler().spectrograms()
    try:
        values = preprocess_spectrograms(spectra, optional_capillary_zones(ctx))
    except (KeyError, TypeError, ValueError) as exc:
        Logger.log_error(
            "Low-frequency spatial spectra failed; check the HD spectrogram "
            "layout and DopplerView segmentation inputs: "
            f"{exc}"
        )
        raise
    ctx.dag.handoff("spatial", values)


def _beats(ctx):
    workflows = velocity_workflows(ctx)
    timing = None
    method = "doppler_moments"
    for method, candidate in (("doppler_moments", workflows.moments), ("frequency_bands", workflows.bandratio)):
        if candidate is not None and candidate.cardiac_cycle_source != "none" and candidate.cycle_boundary_indexes.size >= 2:
            timing = candidate
            break
    if timing is not None and method == "frequency_bands":
        Logger.log_warning("Low-frequency spectra use band-ratio cardiac cycles: detected Doppler moments timing is unavailable.")
    ctx.dag.handoff("beats", run_gap_analysis(ctx, timing, timing_method=method))


def _qc(ctx):
    source = ctx.output.h5.get(EyeFlowOutputPaths.active().spectrograms_masked)
    values = {} if source is None else pack_gap_qc(
        source[()], source_attrs=dict(source.attrs), gap_outputs=ctx.dag.get("beats", {}),
    )
    ctx.dag.handoff("qc", values)


def _figures(ctx):
    if not ctx.output.available:
        ctx.dag.skip("No artifact destination.")
    try:
        dt = ctx.inputs.hd.as_holodoppler().timing().dt_seconds
    except (KeyError, TypeError, ValueError, ZeroDivisionError) as exc:
        Logger.log_warning(f"Cannot calibrate spectrogram time in seconds: {exc}")
        dt = None
    if ctx.can_defer_artifacts:
        ctx.defer_artifacts(
            "spectral_figures", export_spectral_snapshot,
            partial(snapshot_spectral_figures, ctx.output), dt_seconds=dt,
        )
    else:
        export_spectrogram_figures(ctx.output, dt_seconds=dt)


def spectral_graph(*, options=()):
    product = lambda name: (ProductSpec(name, pack=_packed),)
    return StageGraph((
        StageSpec("source", _source, products=(ProductSpec("source"),),
                  inputs=(InputUse("hd.spectrograms", InputAccess.METADATA, optional=True),)),
        StageSpec("sampling", _sampling, requires=("source",), products=product("sampling")),
        StageSpec("calibration", _calibration, requires=("source", "sampling"), products=product("calibration")),
        StageSpec("spatial", _spatial, requires=("source",), products=product("spatial"),
                  inputs=(InputUse("hd.spectrograms", InputAccess.STREAM, optional=True),)),
        StageSpec("beats", _beats, requires=("spatial",), optional_requires=("calibration",), products=product("beats")),
        StageSpec("qc", _qc, requires=("spatial",), optional_requires=("beats",), products=product("qc")),
        StageSpec("figures", _figures, requires=("spatial",), optional_requires=("calibration", "beats", "qc"), artifact=True),
    ))


def run_low_freq_spectra(ctx):
    from pipeline_engine.context import require_managed_runtime
    require_managed_runtime(ctx, "low_freq_spectra")
    return spectral_graph(options=ctx.options_for("low_freq_spectra"))
