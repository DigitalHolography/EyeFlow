"""Export Matlab pulse-analysis style PNG artifacts."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from scipy import optimize, signal

from input_output.output_manager import OutputType

if TYPE_CHECKING:
    from calculations.blood_flow_velocity import PerBeatAnalysisResult
    from pipelines.waveform_shape_metrics.runner import WaveformShapeMetricsContext


PULSE_PNG_SUFFIXES = (
    "f_artery_graph.png",
    "f_vein_graph.png",
    "f_vessel_graph.png",
    "df_vessel_graph.png",
    "v_vessel_graph.png",
    "f_bkg_map.png",
    "f_bkg_colorBar.png",
    "df_map_vessel.png",
    "df_colorBar_vessel.png",
    "df_map.png",
    "df_colorBar.png",
    "f_map.png",
    "f_colorBar.png",
    "map_df_vessel.png",
    "find_systoles_indices_artery.png",
    "find_systoles_indices_vein.png",
    "RI_v_artery.png",
    "PI_v_artery.png",
    "RI_v_vein.png",
    "PI_v_vein.png",
    "ArterialWaveformAnalysis_v_artery.png",
    "VenousWaveformAnalysis_v_vein.png",
    "ArterialSpectralAnalysis_v_artery.png",
    "ArterialSpectralAnalysis_Phase_v_artery.png",
    "VenousSpectralAnalysis_v_vein.png",
    "VenousSpectralAnalysis_Phase_v_vein.png",
    "syntheticSpectralAnalysis.png",
    "arterial_venous_msc.png",
    "arterial_venous_correlation.png",
    "detrended_signals.png",
    "lags.png",
    "Transfer_function_Velocity_AV_mod.png",
    "Transfer_function_Velocity_AV_phase.png",
    "ArterialVenous_Delay.png",
    "map_v_Artery.png",
    "map_v_Vein.png",
    "map_v_Vessel.png",
    "colorbar_v_Artery.png",
    "colorbar_v_Vein.png",
    "v_map.png",
    "histogramVelocityartery.png",
    "histogramVelocityvein.png",
    "AVGflowVideoCombined.png",
)

FREQUENCY_DISPLAY_SCALE = np.float32(1.0 / 1000.0)
VELOCITY_DISPLAY_SCALE = np.float32(1.0 / 1000.0)


@dataclass(frozen=True)
class PulseFigureContext:
    output: object
    stem: str
    time: np.ndarray
    dt_seconds: float
    moment0_avg: np.ndarray
    artery_mask: np.ndarray
    vein_mask: np.ndarray
    section_mask: np.ndarray
    analysis: dict[str, object]
    per_beat_result: PerBeatAnalysisResult
    log: object | None = None

    @property
    def artery_section_mask(self) -> np.ndarray:
        return self.artery_mask & self.section_mask

    @property
    def vein_section_mask(self) -> np.ndarray:
        return self.vein_mask & self.section_mask

    @property
    def vessel_section_mask(self) -> np.ndarray:
        return (self.artery_mask | self.vein_mask) & self.section_mask


class FigureWriter:
    def __init__(self, output, stem: str) -> None:
        self.output = output
        self.stem = str(stem)

    def path(self, suffix: str) -> Path:
        filename = f"{self.stem}_{suffix}"
        path = self.output.path_for(OutputType.PNG, filename)
        path.parent.mkdir(parents=True, exist_ok=True)
        return path

    def savefig(self, fig, suffix: str, *, dpi: int = 150) -> Path:
        path = self.path(suffix)
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        _plt().close(fig)
        return path

    def save_image(self, image: np.ndarray, suffix: str) -> Path:
        filename = f"{self.stem}_{suffix}"
        return self.output.write_png(_finite_image(image), filename)


def export_pulse_pngs(
    output,
    context: WaveformShapeMetricsContext,
    per_beat_result: PerBeatAnalysisResult,
    *,
    log=None,
) -> list[str]:
    """Export core pulse-analysis PNGs for a waveform-shape run."""

    if not getattr(output, "available", False):
        return []
    _matplotlib()
    source_data = context.source_data
    analysis = context.dopplerview_analysis
    frame_count = int(np.asarray(analysis["retinal_artery_velocity_signal"]).size)
    pulse_context = PulseFigureContext(
        output=output,
        stem=_output_stem(output),
        time=np.arange(frame_count, dtype=np.float32) * np.float32(source_data.timing.dt_seconds),
        dt_seconds=float(source_data.timing.dt_seconds),
        moment0_avg=np.mean(source_data.moment0, axis=0, dtype=np.float32),
        artery_mask=np.asarray(source_data.retinal_artery_mask, dtype=bool),
        vein_mask=np.asarray(source_data.retinal_vein_mask, dtype=bool),
        section_mask=_section_mask(analysis, source_data.retinal_artery_mask.shape),
        analysis=analysis,
        per_beat_result=per_beat_result,
        log=log,
    )
    writer = FigureWriter(output, pulse_context.stem)
    paths: list[Path] = []
    for exporter in (
        _export_signal_plots,
        _export_maps,
        _export_systole_plots,
        _export_ri_pi_plots,
        _export_waveform_plots,
        _export_spectral_plots,
        _export_correlation_plots,
        _export_final_visualizations,
    ):
        paths.extend(exporter(writer, pulse_context))
    return [str(path) for path in paths]


def _export_signal_plots(writer: FigureWriter, ctx: PulseFigureContext) -> list[Path]:
    paths: list[Path] = []
    f_video = _array_or_none(ctx.analysis.get("fRMS"))
    f_bkg = _array_or_none(ctx.analysis.get("fRMS_bkg"))
    delta = _array_or_none(ctx.analysis.get("deltafRMS"))
    if f_video is not None and f_bkg is not None:
        f_artery = _display_frequency(_masked_signal(f_video, ctx.artery_section_mask))
        f_artery_bkg = _display_frequency(_masked_signal(f_bkg, ctx.artery_section_mask))
        f_vein = _display_frequency(_masked_signal(f_video, ctx.vein_section_mask))
        f_vein_bkg = _display_frequency(_masked_signal(f_bkg, ctx.vein_section_mask))
        f_vessel_bkg = _display_frequency(_masked_signal(f_bkg, ctx.vessel_section_mask))
        paths.append(
            _line_plot(
                writer,
                "f_artery_graph.png",
                ctx.time,
                [(f_artery, "-", "tab:red", "arteries"), (f_artery_bkg, "--", "k", "background")],
                xlabel="Time(s)",
                ylabel="frequency (kHz)",
            )
        )
        paths.append(
            _line_plot(
                writer,
                "f_vein_graph.png",
                ctx.time,
                [(f_vein, "-", "tab:blue", "veins"), (f_vein_bkg, "--", "k", "background")],
                xlabel="Time(s)",
                ylabel="frequency (kHz)",
            )
        )
        paths.append(
            _line_plot(
                writer,
                "f_vessel_graph.png",
                ctx.time,
                [
                    (f_artery, "-", "tab:red", "arteries"),
                    (f_vein, "-", "tab:blue", "veins"),
                    (f_vessel_bkg, "--", "k", "background"),
                ],
                xlabel="Time(s)",
                ylabel="frequency (kHz)",
            )
        )
    else:
        _log(ctx, "Skipping fRMS signal PNGs; fRMS intermediates are unavailable.")

    if delta is not None:
        paths.append(
            _line_plot(
                writer,
                "df_vessel_graph.png",
                ctx.time,
                [
                    (
                        _display_frequency(_masked_signal(delta, ctx.artery_section_mask)),
                        "-",
                        "tab:red",
                        "arteries",
                    ),
                    (
                        _display_frequency(_masked_signal(delta, ctx.vein_section_mask)),
                        "-",
                        "tab:blue",
                        "veins",
                    ),
                ],
                xlabel="Time(s)",
                ylabel="frequency (kHz)",
            )
        )
    else:
        _log(ctx, "Skipping df_vessel_graph.png; deltafRMS is unavailable.")

    paths.append(
        _line_plot(
            writer,
            "v_vessel_graph.png",
            ctx.time,
            [
                (
                    _display_velocity(_vector(ctx.analysis["retinal_artery_velocity_signal"])),
                    "-",
                    "tab:red",
                    "arteries",
                ),
                (
                    _display_velocity(_vector(ctx.analysis["retinal_vein_velocity_signal"])),
                    "-",
                    "tab:blue",
                    "veins",
                ),
            ],
            title="average velocity in arteries and veins",
            xlabel="Time(s)",
            ylabel="Velocity (mm/s)",
        )
    )
    return paths


def _export_maps(writer: FigureWriter, ctx: PulseFigureContext) -> list[Path]:
    paths: list[Path] = []
    f_bkg_avg = _array_or_none(ctx.analysis.get("fRMS_bkg_avg"))
    f_avg = _array_or_none(ctx.analysis.get("fRMS_avg"))
    delta = _array_or_none(ctx.analysis.get("deltafRMS"))
    if f_bkg_avg is not None:
        f_bkg_avg = _display_frequency(f_bkg_avg)
        paths.extend(
            _heatmap_with_colorbar(
                writer,
                f_bkg_avg,
                "f_bkg_map.png",
                "f_bkg_colorBar.png",
                cmap="gray",
                label="background RMS frequency (kHz)",
            )
        )
    if delta is not None:
        df_mean = _display_frequency(np.nanmean(delta, axis=0))
        paths.extend(
            _heatmap_with_colorbar(
                writer,
                np.where(ctx.vessel_section_mask, df_mean, np.nan),
                "df_map_vessel.png",
                "df_colorBar_vessel.png",
                cmap="gray",
                label="Delta Doppler RMS frequency (kHz)",
            )
        )
        paths.extend(
            _heatmap_with_colorbar(
                writer,
                df_mean,
                "df_map.png",
                "df_colorBar.png",
                cmap="gray",
                label="Delta Doppler RMS frequency (kHz)",
            )
        )
        paths.append(
            _image_map(
                writer,
                np.where(ctx.vessel_section_mask, df_mean, 0.0),
                "map_df_vessel.png",
                cmap="turbo",
                colorbar=True,
                label="kHz",
            )
        )
    if f_avg is not None:
        f_avg = _display_frequency(f_avg)
        paths.extend(
            _heatmap_with_colorbar(
                writer,
                f_avg,
                "f_map.png",
                "f_colorBar.png",
                cmap="gray",
                label="RMS frequency (kHz)",
            )
        )
    return paths


def _export_systole_plots(writer: FigureWriter, ctx: PulseFigureContext) -> list[Path]:
    peaks = _safe_indexes(ctx.analysis.get("beat_indices"))
    if peaks.size == 0:
        _log(ctx, "Skipping systole index PNGs; no beat indices are available.")
        return []
    artery = _display_velocity(_vector(ctx.analysis.get("retinal_artery_velocity_signal_filtered")))
    artery_deriv = _display_velocity(
        _vector(ctx.analysis.get("retinal_artery_velocity_signal_derivative"))
    )
    vein = _display_velocity(_vector(ctx.analysis.get("retinal_vein_velocity_signal_filtered")))
    vein_deriv = _display_velocity(
        _vector(ctx.analysis.get("retinal_vein_velocity_signal_derivative"))
    )
    maxima, minima = _cycle_extrema(artery, peaks)
    return [
        _systole_plot(
            writer,
            "find_systoles_indices_artery.png",
            ctx.time,
            artery,
            artery_deriv,
            peaks,
            maxima,
            minima,
        ),
        _systole_plot(
            writer,
            "find_systoles_indices_vein.png",
            ctx.time,
            vein,
            vein_deriv,
            peaks,
            maxima,
            minima,
        ),
    ]


def _export_ri_pi_plots(writer: FigureWriter, ctx: PulseFigureContext) -> list[Path]:
    peaks = _safe_indexes(ctx.analysis.get("beat_indices"))
    paths: list[Path] = []
    for suffix_prefix, signal_values in (
        ("v_artery", _display_velocity(_vector(ctx.analysis["retinal_artery_velocity_signal"]))),
        ("v_vein", _display_velocity(_vector(ctx.analysis["retinal_vein_velocity_signal"]))),
    ):
        cycle = _average_cycle(signal_values, peaks, 60)
        if cycle is None:
            _log(ctx, f"Skipping RI/PI plots for {suffix_prefix}; insufficient beats.")
            continue
        paths.append(
            _ri_pi_plot(
                writer,
                f"RI_{suffix_prefix}.png",
                ctx.time,
                signal_values,
                cycle,
                "RI",
            )
        )
        paths.append(
            _ri_pi_plot(
                writer,
                f"PI_{suffix_prefix}.png",
                ctx.time,
                signal_values,
                cycle,
                "PI",
            )
        )
    return paths


def _export_waveform_plots(writer: FigureWriter, ctx: PulseFigureContext) -> list[Path]:
    peaks = _safe_indexes(ctx.analysis.get("beat_indices"))
    artery = _average_cycle(
        _display_velocity(_vector(ctx.analysis["retinal_artery_velocity_signal"])), peaks, 128
    )
    vein = _average_cycle(
        _display_velocity(_vector(ctx.analysis["retinal_vein_velocity_signal"])), peaks, 128
    )
    paths: list[Path] = []
    if artery is not None:
        paths.append(_arterial_waveform_plot(writer, ctx, artery))
    else:
        _log(ctx, "Skipping arterial waveform PNG; insufficient beats.")
    if vein is not None:
        paths.append(_venous_waveform_plot(writer, ctx, vein))
    else:
        _log(ctx, "Skipping venous waveform PNG; insufficient beats.")
    return paths


def _export_spectral_plots(writer: FigureWriter, ctx: PulseFigureContext) -> list[Path]:
    paths: list[Path] = []
    artery = _display_velocity(_vector(ctx.analysis["retinal_artery_velocity_signal"]))
    vein = _display_velocity(_vector(ctx.analysis["retinal_vein_velocity_signal"]))
    for prefix, vessel_name, values, color in (
        ("Arterial", "artery", artery, "tab:red"),
        ("Venous", "vein", vein, "tab:blue"),
    ):
        spectrum = _spectrum(values, ctx.dt_seconds)
        paths.append(
            _spectrum_plot(
                writer,
                f"{prefix}SpectralAnalysis_v_{vessel_name}.png",
                spectrum,
                color,
            )
        )
        paths.append(
            _phase_spectrum_plot(
                writer,
                f"{prefix}SpectralAnalysis_Phase_v_{vessel_name}.png",
                spectrum,
                color,
            )
        )
    peaks = _safe_indexes(ctx.analysis.get("beat_indices"))
    cycle = _average_cycle(vein, peaks, 128)
    if cycle is None:
        cycle = _average_cycle(artery, peaks, 128)
    if cycle is not None:
        paths.append(_synthetic_spectral_plot(writer, cycle, ctx.dt_seconds, peaks))
    else:
        _log(ctx, "Skipping syntheticSpectralAnalysis.png; insufficient beats.")
    return paths


def _export_correlation_plots(writer: FigureWriter, ctx: PulseFigureContext) -> list[Path]:
    artery = _display_velocity(_vector(ctx.analysis["retinal_artery_velocity_signal"]))
    vein = -_display_velocity(_vector(ctx.analysis["retinal_vein_velocity_signal"]))
    if artery.size != vein.size or artery.size < 3:
        _log(ctx, "Skipping arterial/venous correlation PNGs; signals are incompatible.")
        return []
    corr = _correlation_data(artery, vein, ctx.dt_seconds)
    spectrum = _transfer_function(artery, vein, ctx.dt_seconds)
    paths = [
        _msc_plot(writer, corr),
        _correlation_panel_plot(writer, ctx.time, corr),
        _detrended_plot(writer, ctx.time, corr),
        _lags_plot(writer, corr),
        _transfer_plot(writer, "Transfer_function_Velocity_AV_mod.png", spectrum, phase=False),
        _transfer_plot(writer, "Transfer_function_Velocity_AV_phase.png", spectrum, phase=True),
    ]
    peaks = _safe_indexes(ctx.analysis.get("beat_indices"))
    artery_cycle = _average_cycle(artery, peaks, 128)
    vein_cycle = _average_cycle(vein, peaks, 128)
    if artery_cycle is not None and vein_cycle is not None:
        paths.append(_delay_plot(writer, ctx, artery_cycle, vein_cycle))
    else:
        _log(ctx, "Skipping ArterialVenous_Delay.png; insufficient beats.")
    return paths


def _export_final_visualizations(writer: FigureWriter, ctx: PulseFigureContext) -> list[Path]:
    velocity_map = _array_or_none(ctx.analysis.get("retinal_vessel_velocity"))
    velocity_avg = _array_or_none(ctx.analysis.get("velocity_map_avg"))
    if velocity_avg is None and velocity_map is not None:
        velocity_avg = np.nanmean(velocity_map, axis=0)
    if velocity_avg is None:
        _log(ctx, "Skipping final velocity visualizations; velocity map is unavailable.")
        return []
    velocity_map_display = _display_velocity(velocity_map) if velocity_map is not None else None
    velocity_avg_display = _display_velocity(velocity_avg)
    vmax = float(np.nanmax(np.abs(velocity_avg_display)))
    paths = [
        _image_map(
            writer,
            np.where(ctx.artery_section_mask, velocity_avg_display, 0.0),
            "map_v_Artery.png",
            cmap="Reds",
            colorbar=True,
            label="mm/s",
        ),
        _image_map(
            writer,
            np.where(ctx.vein_section_mask, velocity_avg_display, 0.0),
            "map_v_Vein.png",
            cmap="Blues",
            colorbar=True,
            label="mm/s",
        ),
        _image_map(
            writer,
            np.where(ctx.vessel_section_mask, velocity_avg_display, 0.0),
            "map_v_Vessel.png",
            cmap="turbo",
            colorbar=True,
            label="mm/s",
        ),
        _colorbar(
            writer,
            "colorbar_v_Artery.png",
            cmap="Reds",
            vmin=0.0,
            vmax=vmax,
            label="mm/s",
        ),
        _colorbar(
            writer,
            "colorbar_v_Vein.png",
            cmap="Blues",
            vmin=0.0,
            vmax=vmax,
            label="mm/s",
        ),
    ]
    flow_rgb = _flow_rgb(ctx, velocity_avg_display)
    writer.save_image(flow_rgb, "v_map.png")
    paths.append(writer.path("v_map.png"))
    hist_artery = _histogram_plot(
        writer,
        ctx,
        velocity_map_display,
        ctx.artery_section_mask,
        "histogramVelocityartery.png",
        "artery",
    )
    hist_vein = _histogram_plot(
        writer,
        ctx,
        velocity_map_display,
        ctx.vein_section_mask,
        "histogramVelocityvein.png",
        "vein",
    )
    paths.extend([hist_artery, hist_vein])
    paths.append(_combined_plot(writer, ctx, flow_rgb, velocity_map_display))
    return paths


def _line_plot(
    writer: FigureWriter,
    suffix: str,
    time: np.ndarray,
    series: list[tuple[np.ndarray, str, str, str]],
    *,
    xlabel: str,
    ylabel: str,
    title: str | None = None,
) -> Path:
    fig, ax = _plt().subplots(figsize=(7.5, 3.0))
    for values, style, color, label in series:
        ax.plot(time[: values.size], values, style, color=color, linewidth=2, label=label)
    if title:
        ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if any(label for _, _, _, label in series):
        ax.legend(loc="best")
    _style_axes(ax)
    return writer.savefig(fig, suffix)


def _heatmap_with_colorbar(
    writer: FigureWriter,
    image: np.ndarray,
    image_suffix: str,
    colorbar_suffix: str,
    *,
    cmap: str,
    label: str,
) -> list[Path]:
    values = _finite_image(image)
    writer.save_image(values, image_suffix)
    return [
        writer.path(image_suffix),
        _colorbar(
            writer,
            colorbar_suffix,
            cmap=cmap,
            vmin=float(np.nanmin(values)),
            vmax=float(np.nanmax(values)),
            label=label,
        ),
    ]


def _image_map(
    writer: FigureWriter,
    image: np.ndarray,
    suffix: str,
    *,
    cmap: str,
    colorbar: bool,
    label: str,
) -> Path:
    fig, ax = _plt().subplots(figsize=(4.0, 4.0))
    im = ax.imshow(_finite_image(image), cmap=cmap)
    ax.axis("off")
    ax.set_aspect("equal")
    if colorbar:
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label(label)
    return writer.savefig(fig, suffix, dpi=180)


def _colorbar(
    writer: FigureWriter,
    suffix: str,
    *,
    cmap: str,
    vmin: float,
    vmax: float,
    label: str,
) -> Path:
    fig, ax = _plt().subplots(figsize=(4.5, 0.55))
    fig.subplots_adjust(bottom=0.45)
    norm = _matplotlib().colors.Normalize(
        vmin=vmin,
        vmax=vmax if vmax > vmin else vmin + 1.0,
    )
    cbar = fig.colorbar(
        _matplotlib().cm.ScalarMappable(norm=norm, cmap=cmap),
        cax=ax,
        orientation="horizontal",
    )
    cbar.set_label(label)
    return writer.savefig(fig, suffix, dpi=150)


def _systole_plot(
    writer: FigureWriter,
    suffix: str,
    time: np.ndarray,
    filtered: np.ndarray,
    derivative: np.ndarray,
    peaks: np.ndarray,
    maxima: np.ndarray,
    minima: np.ndarray,
) -> Path:
    fig, ax = _plt().subplots(figsize=(7.5, 4.0))
    ax.plot(time[: derivative.size], derivative, color="0.55", linewidth=1.5)
    ax.plot(time[: filtered.size], filtered, color="k", linewidth=1.5)
    for indexes, color in ((peaks, "k"), (maxima, "tab:red"), (minima, "tab:blue")):
        valid = indexes[(indexes >= 0) & (indexes < filtered.size)]
        if valid.size:
            ax.scatter(time[valid], filtered[valid], color=color, s=28, zorder=3)
            for idx in valid:
                ax.axvline(time[idx], color=color, linestyle="--", linewidth=1.0, alpha=0.75)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Velocity (mm/s)")
    _style_axes(ax)
    return writer.savefig(fig, suffix)


def _ri_pi_plot(
    writer: FigureWriter,
    suffix: str,
    time: np.ndarray,
    values: np.ndarray,
    cycle: np.ndarray,
    metric_name: str,
) -> Path:
    v_max = float(np.nanmax(cycle))
    v_min = float(np.nanmin(cycle))
    v_mean = float(np.nanmean(cycle))
    if metric_name == "RI" and v_max:
        metric = (v_max - v_min) / v_max
    else:
        metric = (v_max - v_min) / v_mean if v_mean else np.nan
    fig, ax = _plt().subplots(figsize=(7.0, 4.0))
    ax.plot(time[: values.size], values, color="k", linewidth=2)
    for y, label in ((v_max, f"{v_max:.1f} mm/s"), (v_min, f"{v_min:.1f} mm/s")):
        ax.axhline(y, color="0.45", linestyle="--", linewidth=1.5)
        ax.text(time[-1], y, label, ha="right", va="bottom", backgroundcolor="w")
    if metric_name == "PI":
        ax.axhline(v_mean, color="0.45", linestyle="--", linewidth=1.5)
    ax.text(
        0.40,
        0.84,
        f"{metric_name} = {metric:.2f}",
        transform=ax.transAxes,
        backgroundcolor="w",
    )
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Velocity (mm/s)")
    _style_axes(ax)
    return writer.savefig(fig, suffix)


def _arterial_waveform_plot(
    writer: FigureWriter,
    ctx: PulseFigureContext,
    cycle: np.ndarray,
) -> Path:
    pulse_time = _cycle_time(ctx, cycle)
    period = _mean_period_seconds(ctx)
    gradient = np.gradient(cycle)
    peak_indexes, properties = signal.find_peaks(
        cycle,
        height=float(np.nanmax(cycle)) * 0.3,
        distance=max(1, cycle.size // 4),
    )
    if peak_indexes.size == 0:
        peak_indexes = np.asarray([int(np.nanargmax(cycle))])
    else:
        order = np.argsort(properties["peak_heights"])[::-1][:2]
        peak_indexes = peak_indexes[order]
    end_min = int(np.nanargmin(cycle[int(0.75 * cycle.size):]) + int(0.75 * cycle.size))
    notch_idx = _dicrotic_notch_index(cycle, peak_indexes)
    padded_time, padded_signal = _padded_cycle(pulse_time, cycle, period)
    _, padded_gradient = _padded_cycle(pulse_time, gradient, period)
    fig, ax = _plt().subplots(figsize=(6.5, 4.0))
    ax.plot(padded_time, padded_gradient, color="0.85", linewidth=2)
    ax.plot(padded_time, padded_signal, color="0.85", linewidth=2)
    ax.plot(pulse_time, gradient, color="0.70", linewidth=2)
    ax.plot(pulse_time, cycle, color="k", linewidth=2)
    primary_peak = int(peak_indexes[0])
    _annotated_vline(ax, pulse_time[primary_peak], f"{pulse_time[primary_peak]:.2f} s")
    _annotated_vline(ax, period, f"{period:.2f} s")
    _annotated_hline(ax, cycle[primary_peak], f"{cycle[primary_peak]:.1f} mm/s")
    _annotated_vline(ax, pulse_time[end_min], f"{pulse_time[end_min]:.2f} s")
    _annotated_hline(
        ax,
        cycle[end_min],
        f"{cycle[end_min]:.1f} mm/s",
        color="tab:blue",
        vertical_alignment="bottom",
    )
    if notch_idx is not None:
        _annotated_vline(ax, pulse_time[notch_idx], f"{pulse_time[notch_idx]:.2f} s")
        _annotated_hline(ax, cycle[notch_idx], f"{cycle[notch_idx]:.1f} mm/s")
    if peak_indexes.size > 1:
        secondary_peak = int(peak_indexes[1])
        _annotated_vline(
            ax,
            pulse_time[secondary_peak],
            f"{pulse_time[secondary_peak]:.2f} s",
        )
        _annotated_hline(
            ax,
            cycle[secondary_peak],
            f"{cycle[secondary_peak]:.1f} mm/s",
            vertical_alignment="top",
        )
    ax.scatter(
        pulse_time[peak_indexes],
        cycle[peak_indexes],
        color="tab:red",
        s=55,
        edgecolor="k",
    )
    ax.scatter(
        [pulse_time[end_min]],
        [cycle[end_min]],
        color="lightcoral",
        s=55,
        edgecolor="k",
    )
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Velocity (mm/s)")
    ax.margins(y=0.18)
    _style_axes(ax)
    return writer.savefig(fig, "ArterialWaveformAnalysis_v_artery.png", dpi=180)


def _venous_waveform_plot(
    writer: FigureWriter,
    ctx: PulseFigureContext,
    cycle: np.ndarray,
) -> Path:
    pulse_time = _cycle_time(ctx, cycle)
    period = _mean_period_seconds(ctx)
    loc_peak = int(np.nanargmax(cycle))
    loc_trough = int(np.nanargmin(cycle))
    padded_time, padded_signal = _padded_cycle(pulse_time, cycle, period)
    fig, ax = _plt().subplots(figsize=(6.5, 4.0))
    ax.plot(padded_time, padded_signal, color="0.85", linewidth=2)
    ax.plot(pulse_time, cycle, color="k", linewidth=2)
    _annotated_vline(ax, pulse_time[loc_peak], f"{pulse_time[loc_peak]:.2f} s")
    _annotated_vline(ax, period, f"{period:.2f} s")
    _annotated_hline(ax, cycle[loc_peak], f"{cycle[loc_peak]:.1f} mm/s")
    _annotated_vline(ax, pulse_time[loc_trough], f"{pulse_time[loc_trough]:.2f} s")
    _annotated_hline(
        ax,
        cycle[loc_trough],
        f"{cycle[loc_trough]:.1f} mm/s",
        color="tab:blue",
        vertical_alignment="bottom",
    )
    ax.scatter(
        [pulse_time[loc_peak]],
        [cycle[loc_peak]],
        color="tab:blue",
        s=55,
        edgecolor="k",
    )
    ax.scatter(
        [pulse_time[loc_trough]],
        [cycle[loc_trough]],
        color="lightskyblue",
        s=55,
        edgecolor="k",
    )
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Velocity (mm/s)")
    ax.margins(y=0.18)
    _style_axes(ax)
    return writer.savefig(fig, "VenousWaveformAnalysis_v_vein.png", dpi=180)


@dataclass(frozen=True)
class SpectrumData:
    frequencies: np.ndarray
    magnitude: np.ndarray
    phase: np.ndarray
    peak_indexes: np.ndarray
    fundamental_hz: float
    heart_rate_bpm: float
    heart_rate_se_bpm: float


def _spectrum(values: np.ndarray, dt_seconds: float) -> SpectrumData:
    clean = _nan_to_mean(values)
    window = np.hamming(clean.size)
    windowed = clean * window
    padded = np.pad(windowed, (0, clean.size * 2))
    fft = np.fft.rfft(padded)
    freq = np.fft.rfftfreq(padded.size, dt_seconds)
    mag = np.abs(fft) / max(float(np.mean(window)), np.finfo(np.float32).eps)
    mag = mag / max(float(np.nanmax(mag)), np.finfo(np.float32).eps)
    peaks, _ = signal.find_peaks(
        mag,
        distance=max(1, int(0.25 / dt_seconds)),
        prominence=max(float(np.nanmax(mag)) * 0.05, np.finfo(np.float32).eps),
    )
    if peaks.size == 0:
        peaks, _ = signal.find_peaks(
            mag,
            distance=max(1, int(0.25 / dt_seconds)),
        )
    peaks = peaks[np.argsort(freq[peaks])]
    fundamental, heart_rate, heart_rate_se = _spectral_heart_rate(freq, peaks)
    return SpectrumData(
        freq.astype(np.float32),
        mag.astype(np.float32),
        np.angle(fft).astype(np.float32),
        peaks.astype(np.int32),
        fundamental,
        heart_rate,
        heart_rate_se,
    )


def _spectrum_plot(
    writer: FigureWriter,
    suffix: str,
    spectrum: SpectrumData,
    color: str,
) -> Path:
    fig, ax = _plt().subplots(figsize=(6.5, 4.0))
    ax.plot(spectrum.frequencies, spectrum.magnitude, color="k", linewidth=2)
    peaks = spectrum.peak_indexes[:8]
    if peaks.size:
        ax.scatter(
            spectrum.frequencies[peaks],
            spectrum.magnitude[peaks],
            color=color,
            s=45,
            edgecolor="k",
        )
        _annotate_spectral_peaks(ax, spectrum, peaks, phase=False)
    ax.set_xlim(0, 10)
    ax.set_ylim(0, max(1.0, float(np.nanmax(spectrum.magnitude)) * 1.15))
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Normalized Magnitude")
    if np.isfinite(spectrum.heart_rate_bpm):
        ax.text(
            0.52,
            0.62,
            (
                f"HR : {spectrum.heart_rate_bpm:.1f} BPM "
                f"+/- {spectrum.heart_rate_se_bpm:.1f}"
            ),
            transform=ax.transAxes,
            fontsize=10,
            backgroundcolor="w",
        )
    _style_axes(ax, grid=True)
    return writer.savefig(fig, suffix, dpi=180)


def _phase_spectrum_plot(
    writer: FigureWriter,
    suffix: str,
    spectrum: SpectrumData,
    color: str,
) -> Path:
    fig, ax = _plt().subplots(figsize=(6.5, 4.0))
    ax.plot(spectrum.frequencies, spectrum.phase, color="k", linewidth=2)
    peaks = spectrum.peak_indexes[:8]
    if peaks.size:
        ax.scatter(
            spectrum.frequencies[peaks],
            spectrum.phase[peaks],
            color=color,
            s=45,
            edgecolor="k",
        )
        _annotate_spectral_peaks(ax, spectrum, peaks, phase=True)
    ax.set_xlim(0, 10)
    ax.set_ylim(-np.pi, np.pi)
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Phase (rad)")
    _style_axes(ax, grid=True)
    return writer.savefig(fig, suffix, dpi=180)


def _dicrotic_notch_index(cycle: np.ndarray, peak_indexes: np.ndarray) -> int | None:
    if peak_indexes.size < 2:
        return None
    first, second = sorted(int(index) for index in peak_indexes[:2])
    if second <= first:
        return None
    notch = first + int(np.nanargmin(cycle[first : second + 1]))
    return notch if second - notch > cycle.size * 0.05 else None


def _padded_cycle(
    pulse_time: np.ndarray,
    values: np.ndarray,
    period_seconds: float,
) -> tuple[np.ndarray, np.ndarray]:
    half = values.size // 2
    if half == 0:
        return pulse_time, values
    left_time = np.linspace(-period_seconds / 2.0, 0.0, half, endpoint=False)
    right_count = values.size - half
    right_time = np.linspace(
        period_seconds,
        period_seconds * 1.5,
        right_count,
        endpoint=False,
    )
    return (
        np.concatenate((left_time, pulse_time, right_time)).astype(np.float32),
        np.concatenate((values[-half:], values, values[:right_count])).astype(np.float32),
    )


def _annotated_vline(ax, x: float, label: str) -> None:
    ax.axvline(x, color="0.25", linestyle="--", linewidth=1.5, zorder=1)
    ax.text(
        x,
        0.04,
        label,
        transform=ax.get_xaxis_transform(),
        color="0.25",
        fontsize=9,
        ha="center",
        va="bottom",
        backgroundcolor="w",
    )


def _annotated_hline(
    ax,
    y: float,
    label: str,
    *,
    color: str = "0.25",
    vertical_alignment: str = "top",
) -> None:
    ax.axhline(y, color=color, linestyle="--", linewidth=1.5, zorder=1)
    ax.text(
        0.98,
        y,
        label,
        transform=ax.get_yaxis_transform(),
        color="0.25",
        fontsize=9,
        ha="right",
        va=vertical_alignment,
        backgroundcolor="w",
    )


def _spectral_heart_rate(
    frequencies: np.ndarray,
    peak_indexes: np.ndarray,
) -> tuple[float, float, float]:
    if peak_indexes.size == 0:
        return np.nan, np.nan, np.nan
    fundamental = float(frequencies[peak_indexes[0]])
    if not np.isfinite(fundamental) or fundamental <= 0:
        return fundamental, np.nan, np.nan
    peak_freqs = frequencies[peak_indexes]
    harmonic_numbers = np.rint(peak_freqs / fundamental)
    valid = harmonic_numbers > 0
    if not np.any(valid):
        return fundamental, np.nan, np.nan
    numbers = harmonic_numbers[valid]
    locations = peak_freqs[valid]
    heart_rate_hz = float(np.dot(numbers, locations) / np.dot(numbers, numbers))
    residuals = locations - numbers * heart_rate_hz
    heart_rate_se_hz = float(np.sqrt(np.mean(residuals**2)) / np.sqrt(locations.size))
    return fundamental, heart_rate_hz * 60.0, heart_rate_se_hz * 60.0


def _annotate_spectral_peaks(
    ax,
    spectrum: SpectrumData,
    peak_indexes: np.ndarray,
    *,
    phase: bool,
) -> None:
    if not np.isfinite(spectrum.fundamental_hz) or spectrum.fundamental_hz <= 0:
        return
    y_values = spectrum.phase if phase else spectrum.magnitude
    y_offset = 0.30 if phase else max(0.06, float(np.nanmax(spectrum.magnitude)) * 0.08)
    for index in peak_indexes:
        freq = float(spectrum.frequencies[index])
        harmonic = int(round(freq / spectrum.fundamental_hz))
        if not phase:
            ax.axvline(freq, color="0.5", linestyle="--", linewidth=1.0, zorder=1)
            ax.text(
                freq,
                0.04,
                f"{harmonic}x",
                transform=ax.get_xaxis_transform(),
                color="0.35",
                fontsize=9,
                ha="center",
                va="bottom",
                backgroundcolor="w",
            )
        ax.text(
            freq,
            float(y_values[index]) + y_offset,
            f"{freq:.2f}",
            fontsize=9,
            ha="center",
            va="bottom",
            backgroundcolor="w",
        )


def _synthetic_spectral_plot(
    writer: FigureWriter,
    cycle: np.ndarray,
    dt_seconds: float,
    beat_indexes: np.ndarray,
) -> Path:
    repeated = np.tile(cycle, 512)
    if beat_indexes.size > 1:
        period = float(np.nanmean(np.diff(beat_indexes)) * dt_seconds)
    else:
        period = cycle.size * dt_seconds
    fs = cycle.size / max(period, np.finfo(np.float32).eps)
    frequencies = np.fft.rfftfreq(repeated.size, 1.0 / fs)
    fft = np.fft.rfft(repeated)
    magnitude = np.abs(fft)
    magnitude = (
        magnitude
        / max(float(magnitude[0]), np.finfo(np.float32).eps)
        * float(np.nanmean(repeated))
    )
    phase = np.angle(fft)
    peaks, _ = signal.find_peaks(magnitude, height=0.001)
    fig, axes = _plt().subplots(2, 1, figsize=(7.0, 5.0), sharex=True)
    axes[0].plot(frequencies, magnitude, color="k", linewidth=2)
    axes[0].scatter(frequencies[peaks], magnitude[peaks], s=25)
    axes[0].set_ylabel("Magnitude |Y|")
    axes[1].plot(frequencies, phase / np.pi, color="k", linewidth=2)
    axes[1].scatter(frequencies[peaks], phase[peaks] / np.pi, s=25)
    axes[1].set_ylabel("Phase / pi")
    axes[1].set_xlabel("Frequency (Hz)")
    for ax in axes:
        ax.set_xlim(0, 10)
        _style_axes(ax)
    axes[1].set_ylim(-1, 1)
    return writer.savefig(fig, "syntheticSpectralAnalysis.png", dpi=180)


@dataclass(frozen=True)
class CorrelationData:
    artery: np.ndarray
    vein: np.ndarray
    lags_seconds: np.ndarray
    cross_corr: np.ndarray
    time_lag: float
    max_corr: float
    coherence_freq: np.ndarray
    coherence: np.ndarray
    heart_rate_hz: float
    gamma_0: float


def _correlation_data(artery: np.ndarray, vein: np.ndarray, dt_seconds: float) -> CorrelationData:
    a = _standardize(artery)
    v = _standardize(vein)
    cross = signal.correlate(a, v, mode="full")
    denom = max(float(a.size * np.nanstd(a) * np.nanstd(v)), np.finfo(np.float32).eps)
    cross = cross / denom
    lags = signal.correlation_lags(a.size, v.size, mode="full")
    max_idx = int(np.nanargmax(cross))
    fs = 1.0 / dt_seconds
    nperseg = min(64, a.size)
    freqs, coherence = signal.coherence(
        a,
        v,
        fs=fs,
        window=np.hamming(nperseg),
        nperseg=nperseg,
        noverlap=None,
        nfft=max(256, nperseg),
    )
    spectrum = _spectrum(a, dt_seconds)
    heart_rate = (
        float(spectrum.frequencies[spectrum.peak_indexes[0]])
        if spectrum.peak_indexes.size
        else np.nan
    )
    gamma_0 = _gamma_0(coherence, freqs, heart_rate)
    return CorrelationData(
        a.astype(np.float32),
        v.astype(np.float32),
        (lags * dt_seconds).astype(np.float32),
        cross.astype(np.float32),
        float(lags[max_idx] * dt_seconds),
        float(cross[max_idx]),
        freqs.astype(np.float32),
        coherence.astype(np.float32),
        heart_rate,
        gamma_0,
    )


def _msc_plot(writer: FigureWriter, corr: CorrelationData) -> Path:
    fig, ax = _plt().subplots(figsize=(6.5, 4.0))
    ax.plot(corr.coherence_freq, corr.coherence, color="k", linewidth=2)
    if np.isfinite(corr.heart_rate_hz):
        ax.axvline(corr.heart_rate_hz, color="k", linestyle="--", linewidth=1.5)
        ax.text(
            corr.heart_rate_hz,
            0.04,
            f"{corr.heart_rate_hz:.2f} Hz",
            transform=ax.get_xaxis_transform(),
            fontsize=9,
            ha="center",
            va="bottom",
            backgroundcolor="w",
        )
    if np.isfinite(corr.gamma_0):
        ax.text(
            0.62,
            0.62,
            rf"$\Gamma_0 = {corr.gamma_0:.2f}$",
            transform=ax.transAxes,
            fontsize=12,
            backgroundcolor="w",
        )
    ax.set_xlim(0, 10)
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Magnitude-squared coherence")
    _style_axes(ax)
    return writer.savefig(fig, "arterial_venous_msc.png")


def _gamma_0(coherence: np.ndarray, frequencies: np.ndarray, heart_rate_hz: float) -> float:
    if not np.isfinite(heart_rate_hz):
        return np.nan
    valid = (frequencies < heart_rate_hz + 0.3) & (frequencies > heart_rate_hz - 0.3)
    if not np.any(valid):
        return np.nan
    return float(np.nanmean(coherence[valid]))


def _correlation_panel_plot(writer: FigureWriter, time: np.ndarray, corr: CorrelationData) -> Path:
    fig, axes = _plt().subplots(2, 1, figsize=(7.0, 5.2))
    axes[0].plot(time[: corr.artery.size], corr.artery, color="tab:red", linewidth=2)
    axes[0].plot(time[: corr.vein.size], -corr.vein, color="tab:blue", linewidth=2)
    axes[0].set_ylabel("Amplitude")
    axes[1].plot(corr.lags_seconds, corr.cross_corr, color="k", linewidth=1.5)
    axes[1].scatter([corr.time_lag], [corr.max_corr], color="tab:red", s=55)
    axes[1].legend([f"Peak Lag: {corr.time_lag:.3f} s", f"Peak Corr: {corr.max_corr:.2f}"])
    axes[1].set_xlabel("Lag (s)")
    axes[1].set_ylabel("Cross-Correlation")
    for ax in axes:
        _style_axes(ax, grid=True)
    return writer.savefig(fig, "arterial_venous_correlation.png")


def _detrended_plot(writer: FigureWriter, time: np.ndarray, corr: CorrelationData) -> Path:
    return _line_plot(
        writer,
        "detrended_signals.png",
        time,
        [(corr.artery, "-", "tab:red", "artery"), (-corr.vein, "-", "tab:blue", "vein")],
        xlabel="Time (s)",
        ylabel="Amplitude",
    )


def _lags_plot(writer: FigureWriter, corr: CorrelationData) -> Path:
    fig, ax = _plt().subplots(figsize=(6.5, 4.0))
    ax.plot(corr.lags_seconds, corr.cross_corr, color="k", linewidth=1.5)
    ax.scatter([corr.time_lag], [corr.max_corr], color="tab:red", s=55)
    ax.legend([f"Peak Corr: {corr.max_corr:.2f}", f"Peak Lag: {corr.time_lag:.3f} s"])
    ax.set_xlabel("Lag (s)")
    ax.set_ylabel("Cross-Correlation")
    _style_axes(ax, grid=True)
    return writer.savefig(fig, "lags.png")


@dataclass(frozen=True)
class TransferData:
    frequencies: np.ndarray
    transfer: np.ndarray


def _transfer_function(artery: np.ndarray, vein: np.ndarray, dt_seconds: float) -> TransferData:
    n = max(artery.size, 1) * 10
    fft_a = np.fft.fft(_nan_to_mean(artery), n)
    fft_v = np.fft.fft(_nan_to_mean(vein), n)
    transfer = fft_v / np.where(np.abs(fft_a) == 0, np.nan + 0j, fft_a)
    freqs = np.fft.fftfreq(n, d=dt_seconds)
    order = np.argsort(freqs)
    return TransferData(freqs[order].astype(np.float32), transfer[order].astype(np.complex64))


def _transfer_plot(writer: FigureWriter, suffix: str, data: TransferData, *, phase: bool) -> Path:
    fig, ax = _plt().subplots(figsize=(6.5, 4.0))
    positive = data.frequencies >= 0
    x = data.frequencies[positive]
    if phase:
        y = np.angle(data.transfer[positive])
        ylabel = "transfer function angle"
        ax.plot(x, y, color="k", linewidth=2)
    else:
        y = np.abs(data.transfer[positive])
        ylabel = "transfer function"
        ax.semilogy(x, y, color="k", linewidth=2)
    ax.set_xlim(0, 10)
    ax.set_xlabel("Freq (Hz)")
    ax.set_ylabel(ylabel)
    _style_axes(ax, grid=True)
    return writer.savefig(fig, suffix)


def _delay_plot(
    writer: FigureWriter,
    ctx: PulseFigureContext,
    artery_cycle: np.ndarray,
    vein_cycle: np.ndarray,
) -> Path:
    artery = _rescale(artery_cycle)
    vein = _rescale(vein_cycle)
    amin = int(np.nanargmin(vein_cycle))
    t = np.arange(artery.size, dtype=np.float32)

    def model_tau(tau):
        kernel = np.exp(-t / max(tau, 1e-6)) / max(tau, 1e-6)
        result = np.convolve(artery, kernel, mode="full")[: artery.size]
        return _rescale(np.roll(result, amin))

    target = np.zeros(vein.size, dtype=np.float32)
    target[amin:] = 1.0

    def objective(tau):
        return float(np.sum((_rescale(vein) - model_tau(float(tau[0]))) ** 2 * target))

    fitted = optimize.minimize(objective, x0=np.asarray([40.0]), method="Nelder-Mead")
    tau = float(fitted.x[0]) if fitted.success else 40.0
    period = _mean_period_seconds(ctx)
    ti = np.linspace(0, period, artery.size)
    fig, ax = _plt().subplots(figsize=(6.5, 4.0))
    ax.plot(ti, artery, color="tab:red", linewidth=2)
    ax.plot(ti, model_tau(tau), color="tab:blue", linestyle="--", linewidth=2)
    ax.plot(ti, vein * float(np.nanmax(model_tau(tau))), color="tab:blue", linewidth=2)
    ax.legend(
        [
            "Artery (normalized)",
            f"Vein model fit (tau_RC = {tau / artery.size * period * 1000:.2f} ms)",
            f"Vein shifted ({amin / artery.size * period * 1000:.2f} ms)",
        ],
        loc="best",
    )
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Arterio-venous decay fit")
    _style_axes(ax)
    return writer.savefig(fig, "ArterialVenous_Delay.png", dpi=180)


def _histogram_plot(
    writer: FigureWriter,
    ctx: PulseFigureContext,
    velocity_map: np.ndarray | None,
    mask: np.ndarray,
    suffix: str,
    vessel: str,
) -> Path:
    if velocity_map is None:
        velocity_map = _display_velocity(
            np.asarray(ctx.analysis["velocity_map_avg"], dtype=np.float32)
        )[None, :, :]
    histo = _histogram_matrix(velocity_map, mask)
    fig, ax = _plt().subplots(figsize=(6.0, 2.75))
    extent = [0, velocity_map.shape[0] * ctx.dt_seconds, histo.vmin, histo.vmax]
    ax.imshow(
        histo.counts,
        aspect="auto",
        origin="lower",
        extent=extent,
        cmap=_vessel_histogram_colormap(vessel),
        interpolation="nearest",
        vmin=0,
        vmax=histo.count_max,
    )
    ax.set_facecolor("black")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Velocity (mm/s)")
    _style_axes(ax)
    return writer.savefig(fig, suffix)


@dataclass(frozen=True)
class HistogramData:
    counts: np.ndarray
    vmin: float
    vmax: float
    count_max: float


def _histogram_matrix(velocity_map: np.ndarray, mask: np.ndarray, bins: int = 256) -> HistogramData:
    data = np.asarray(velocity_map, dtype=np.float32)
    selected = data[:, mask]
    selected = selected[np.isfinite(selected)]
    if selected.size == 0:
        return HistogramData(
            np.zeros((bins, data.shape[0]), dtype=np.float32),
            0.0,
            1.0,
            1.0,
        )
    vmin = float(np.nanmin(selected))
    vmax = float(np.nanmax(selected))
    if vmax <= vmin:
        vmax = vmin + 1.0
    counts = np.zeros((bins, data.shape[0]), dtype=np.float32)
    edges = np.linspace(vmin, vmax, bins + 1, dtype=np.float32)
    for frame_idx, frame in enumerate(data):
        counts[:, frame_idx] = np.histogram(frame[mask], bins=edges)[0]
    count_max = float(np.nanmax(counts))
    return HistogramData(counts, vmin, vmax, count_max if count_max > 0 else 1.0)


def _vessel_histogram_colormap(vessel: str):
    if vessel == "vein":
        stops = (
            ([0.0, 0.0, 0.0], 0.0),
            ([0.0, 0.0, 1.0], 1.0 / 3.0),
            ([0.0, 1.0, 1.0], 2.0 / 3.0),
            ([1.0, 1.0, 1.0], 1.0),
        )
    else:
        stops = (
            ([0.0, 0.0, 0.0], 0.0),
            ([1.0, 0.0, 0.0], 1.0 / 3.0),
            ([1.0, 1.0, 0.0], 2.0 / 3.0),
            ([1.0, 1.0, 1.0], 1.0),
        )
    return _matplotlib().colors.ListedColormap(_cmap_lab(256, stops))


def _cmap_lab(n: int, stops) -> np.ndarray:
    from skimage.color import lab2rgb, rgb2lab

    lab = np.zeros((n, 3), dtype=np.float32)
    for (rgb1, pos1), (rgb2, pos2) in zip(stops[:-1], stops[1:]):
        first = int(round(n * pos1))
        last = max(first, int(round(n * pos2)) - 1)
        lab1 = rgb2lab(np.asarray(rgb1, dtype=np.float32).reshape(1, 1, 3))[0, 0]
        lab2 = rgb2lab(np.asarray(rgb2, dtype=np.float32).reshape(1, 1, 3))[0, 0]
        x = np.linspace(0.0, 1.0, last - first + 1, dtype=np.float32)
        lab[first : last + 1] = lab1 + (lab2 - lab1) * x[:, None]
    cmap = lab2rgb(lab.reshape(n, 1, 3)).reshape(n, 3)
    cmap = np.clip(cmap, 0.0, 1.0)
    cmap[~np.isfinite(cmap)] = 0.0
    return cmap.astype(np.float32)


def _combined_plot(
    writer: FigureWriter,
    ctx: PulseFigureContext,
    flow_rgb: np.ndarray,
    velocity_map: np.ndarray | None,
) -> Path:
    fig, axes = _plt().subplots(2, 2, figsize=(8.0, 5.5))
    axes[0, 0].imshow(flow_rgb)
    axes[0, 0].axis("off")
    axes[1, 0].axis("off")
    for ax, mask, title, vessel in (
        (axes[0, 1], ctx.artery_section_mask, "Artery", "artery"),
        (axes[1, 1], ctx.vein_section_mask, "Vein", "vein"),
    ):
        histogram_source = (
            velocity_map
            if velocity_map is not None
            else _display_velocity(np.asarray(ctx.analysis["velocity_map_avg"]))[None, :, :]
        )
        histo = _histogram_matrix(
            histogram_source,
            mask,
        )
        extent = [0, histo.counts.shape[1] * ctx.dt_seconds, histo.vmin, histo.vmax]
        ax.imshow(
            histo.counts,
            aspect="auto",
            origin="lower",
            extent=extent,
            cmap=_vessel_histogram_colormap(vessel),
            interpolation="nearest",
            vmin=0,
            vmax=histo.count_max,
        )
        ax.set_facecolor("black")
        ax.set_title(title)
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Velocity (mm/s)")
        _style_axes(ax)
    return writer.savefig(fig, "AVGflowVideoCombined.png", dpi=150)


def _flow_rgb(ctx: PulseFigureContext, velocity_avg: np.ndarray) -> np.ndarray:
    background = _rescale(ctx.moment0_avg)
    rgb = np.dstack([background, background, background])
    value = _rescale(np.abs(velocity_avg))
    artery_color = _plt().get_cmap("Reds")(value)[..., :3]
    vein_color = _plt().get_cmap("Blues")(value)[..., :3]
    vessel_color = _plt().get_cmap("Purples")(value)[..., :3]
    mask_av = ctx.artery_section_mask & ctx.vein_section_mask
    artery_only = ctx.artery_section_mask & ~mask_av
    vein_only = ctx.vein_section_mask & ~mask_av
    rgb[artery_only] = artery_color[artery_only]
    rgb[vein_only] = vein_color[vein_only]
    rgb[mask_av] = vessel_color[mask_av]
    return np.clip(rgb, 0, 1)


def _masked_signal(video: np.ndarray, mask: np.ndarray) -> np.ndarray:
    data = np.asarray(video, dtype=np.float32)
    count = np.count_nonzero(mask)
    if count == 0:
        return np.full((data.shape[0],), np.nan, dtype=np.float32)
    return np.nanmean(np.where(mask[None, :, :], data, np.nan), axis=(1, 2)).astype(np.float32)


def _cycle_extrema(values: np.ndarray, peaks: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    maxima: list[int] = []
    minima: list[int] = []
    if peaks.size == 0:
        return np.asarray(maxima, dtype=np.int32), np.asarray(minima, dtype=np.int32)
    for first, second in zip(peaks[:-1], peaks[1:]):
        if second <= first + 2:
            continue
        midpoint = first + int(round((second - first) / 2))
        maxima.append(first + int(np.nanargmax(values[first:midpoint])))
        minima.append(midpoint + int(np.nanargmin(values[midpoint:second])))
    if peaks[0] > 0:
        minima.insert(0, int(np.nanargmin(values[: peaks[0]])))
    if peaks[-1] < values.size:
        maxima.append(peaks[-1] + int(np.nanargmax(values[peaks[-1]:])))
    return np.asarray(maxima, dtype=np.int32), np.asarray(minima, dtype=np.int32)


def _average_cycle(values: np.ndarray, peaks: np.ndarray, samples: int) -> np.ndarray | None:
    if peaks.size < 2:
        return None
    cycles = []
    base_x = np.linspace(0, 1, samples, dtype=np.float32)
    for start, stop in zip(peaks[:-1], peaks[1:]):
        if stop <= start + 1:
            continue
        beat = values[start:stop]
        cycles.append(np.interp(base_x, np.linspace(0, 1, beat.size), beat))
    if not cycles:
        return None
    return np.nanmean(np.asarray(cycles, dtype=np.float32), axis=0).astype(np.float32)


def _cycle_time(ctx: PulseFigureContext, cycle: np.ndarray) -> np.ndarray:
    return np.linspace(0, _mean_period_seconds(ctx), cycle.size, dtype=np.float32)


def _mean_period_seconds(ctx: PulseFigureContext) -> float:
    peaks = _safe_indexes(ctx.analysis.get("beat_indices"))
    if peaks.size > 1:
        return float(np.nanmean(np.diff(peaks)) * ctx.dt_seconds)
    return float(ctx.dt_seconds * 128)


def _standardize(values: np.ndarray) -> np.ndarray:
    clean = _nan_to_mean(values)
    std = float(np.nanstd(clean))
    if std <= 0:
        return np.zeros_like(clean)
    return ((clean - float(np.nanmean(clean))) / std).astype(np.float32)


def _nan_to_mean(values: np.ndarray) -> np.ndarray:
    clean = np.asarray(values, dtype=np.float32).copy()
    finite = np.isfinite(clean)
    if not np.any(finite):
        return np.zeros_like(clean)
    clean[~finite] = float(np.nanmean(clean[finite]))
    return clean


def _rescale(values: np.ndarray) -> np.ndarray:
    data = _finite_image(values).astype(np.float32)
    vmin = float(np.nanmin(data))
    vmax = float(np.nanmax(data))
    if vmax <= vmin:
        return np.zeros_like(data)
    return (data - vmin) / (vmax - vmin)


def _finite_image(image: np.ndarray) -> np.ndarray:
    data = np.asarray(image, dtype=np.float32)
    finite = np.isfinite(data)
    if not np.any(finite):
        return np.zeros_like(data)
    clean = data.copy()
    clean[~finite] = float(np.nanmin(data[finite]))
    return clean


def _safe_indexes(values) -> np.ndarray:
    if values is None:
        return np.asarray([], dtype=np.int32)
    indexes = np.asarray(values, dtype=np.int32).reshape(-1)
    return indexes[indexes >= 0]


def _vector(values) -> np.ndarray:
    if values is None:
        return np.asarray([], dtype=np.float32)
    return np.asarray(values, dtype=np.float32).reshape(-1)


def _display_frequency(values) -> np.ndarray:
    return np.asarray(values, dtype=np.float32) * FREQUENCY_DISPLAY_SCALE


def _display_velocity(values) -> np.ndarray:
    return np.asarray(values, dtype=np.float32) * VELOCITY_DISPLAY_SCALE


def _array_or_none(values) -> np.ndarray | None:
    if values is None:
        return None
    return np.asarray(values, dtype=np.float32)


def _section_mask(analysis: dict[str, object], shape: tuple[int, int]) -> np.ndarray:
    section = analysis.get("velocity_section_mask")
    if section is not None:
        return np.asarray(section, dtype=bool)
    return np.ones(shape, dtype=bool)


def _output_stem(output) -> str:
    manager = getattr(output, "manager", None)
    layout = getattr(manager, "layout", None)
    stem = getattr(layout, "stem", None)
    return str(stem or "eyeflow")


def _style_axes(ax, *, grid: bool = False) -> None:
    ax.set_box_aspect(1 / 1.618)
    ax.grid(grid)
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)


def _log(ctx: PulseFigureContext, message: str) -> None:
    if callable(ctx.log):
        ctx.log(message)


def _matplotlib():
    import matplotlib

    matplotlib.use("Agg", force=True)
    return matplotlib


def _plt():
    _matplotlib()
    import matplotlib.pyplot as plt

    return plt
