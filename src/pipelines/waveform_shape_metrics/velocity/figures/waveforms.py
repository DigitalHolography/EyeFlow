"""Per-cycle waveform, RI, and PI PNG exporters."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from calculations.blood_flow_velocity.signal_analysis.waveform import (
    ArterialWaveformAnalysis,
    VenousWaveformAnalysis,
    arterial_waveform_analysis,
    pulse_metric,
    vessel_cycle_analysis,
    venous_waveform_analysis,
)
from input_output.writers.png import PngArtifactWriter as FigureWriter

from .common import (
    PulseFigureContext,
    _log,
    _plt,
    _safe_indexes,
    _vector,
    display_velocity as _display_velocity,
)
from .plotting import _annotated_hline, _annotated_vline, _style_axes


def _export_ri_pi_plots(writer: FigureWriter, ctx: PulseFigureContext) -> list[Path]:
    peaks = _safe_indexes(ctx.analysis.get("beat_indices"))
    cycles = vessel_cycle_analysis(
        _display_velocity(_vector(ctx.analysis["retinal_artery_velocity_signal"])),
        _display_velocity(_vector(ctx.analysis["retinal_vein_velocity_signal"])),
        peaks,
        60,
    )
    paths: list[Path] = []
    for suffix_prefix, signal_values in (
        ("v_artery", _display_velocity(_vector(ctx.analysis["retinal_artery_velocity_signal"]))),
        ("v_vein", _display_velocity(_vector(ctx.analysis["retinal_vein_velocity_signal"]))),
    ):
        cycle = cycles.artery if suffix_prefix == "v_artery" else cycles.vein
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
    cycles = vessel_cycle_analysis(
        _display_velocity(_vector(ctx.analysis["retinal_artery_velocity_signal"])),
        _display_velocity(_vector(ctx.analysis["retinal_vein_velocity_signal"])),
        peaks,
        128,
    )
    paths: list[Path] = []
    if cycles.artery is not None:
        paths.append(_arterial_waveform_plot(writer, ctx, cycles.artery))
    else:
        _log(ctx, "Skipping arterial waveform PNG; insufficient beats.")
    if cycles.vein is not None:
        paths.append(_venous_waveform_plot(writer, ctx, cycles.vein))
    else:
        _log(ctx, "Skipping venous waveform PNG; insufficient beats.")
    return paths

def _ri_pi_plot(
    writer: FigureWriter,
    suffix: str,
    time: np.ndarray,
    values: np.ndarray,
    cycle: np.ndarray,
    metric_name: str,
) -> Path:
    metric = pulse_metric(cycle, metric_name)
    fig, ax = _plt().subplots(figsize=(7.0, 4.0))
    ax.plot(time[: values.size], values, color="k", linewidth=2)
    for y, label in (
        (metric.maximum, f"{metric.maximum:.1f} mm/s"),
        (metric.minimum, f"{metric.minimum:.1f} mm/s"),
    ):
        ax.axhline(y, color="0.45", linestyle="--", linewidth=1.5)
        ax.text(time[-1], y, label, ha="right", va="bottom", backgroundcolor="w")
    if metric_name == "PI":
        ax.axhline(metric.mean, color="0.45", linestyle="--", linewidth=1.5)
    ax.text(
        0.40,
        0.84,
        f"{metric_name} = {metric.value:.2f}",
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
    data = arterial_waveform_analysis(
        cycle,
        _safe_indexes(ctx.analysis.get("beat_indices")),
        ctx.dt_seconds,
    )
    return _arterial_waveform_analysis_plot(writer, cycle, data)


def _arterial_waveform_analysis_plot(
    writer: FigureWriter,
    cycle: np.ndarray,
    data: ArterialWaveformAnalysis,
) -> Path:
    fig, ax = _plt().subplots(figsize=(6.5, 4.0))
    ax.plot(data.padded_time, data.padded_gradient, color="0.85", linewidth=2)
    ax.plot(data.padded_time, data.padded_signal, color="0.85", linewidth=2)
    ax.plot(data.pulse_time, data.gradient, color="0.70", linewidth=2)
    ax.plot(data.pulse_time, cycle, color="k", linewidth=2)
    primary_peak = int(data.peak_indexes[0])
    _annotated_vline(
        ax,
        data.pulse_time[primary_peak],
        f"{data.pulse_time[primary_peak]:.2f} s",
    )
    _annotated_vline(ax, data.period_seconds, f"{data.period_seconds:.2f} s")
    _annotated_hline(ax, cycle[primary_peak], f"{cycle[primary_peak]:.1f} mm/s")
    _annotated_vline(
        ax,
        data.pulse_time[data.end_min_index],
        f"{data.pulse_time[data.end_min_index]:.2f} s",
    )
    _annotated_hline(
        ax,
        cycle[data.end_min_index],
        f"{cycle[data.end_min_index]:.1f} mm/s",
        color="tab:blue",
        vertical_alignment="bottom",
    )
    if data.notch_index is not None:
        _annotated_vline(
            ax,
            data.pulse_time[data.notch_index],
            f"{data.pulse_time[data.notch_index]:.2f} s",
        )
        _annotated_hline(ax, cycle[data.notch_index], f"{cycle[data.notch_index]:.1f} mm/s")
    if data.peak_indexes.size > 1:
        secondary_peak = int(data.peak_indexes[1])
        _annotated_vline(
            ax,
            data.pulse_time[secondary_peak],
            f"{data.pulse_time[secondary_peak]:.2f} s",
        )
        _annotated_hline(
            ax,
            cycle[secondary_peak],
            f"{cycle[secondary_peak]:.1f} mm/s",
            vertical_alignment="top",
        )
    ax.scatter(
        data.pulse_time[data.peak_indexes],
        cycle[data.peak_indexes],
        color="tab:red",
        s=55,
        edgecolor="k",
    )
    ax.scatter(
        [data.pulse_time[data.end_min_index]],
        [cycle[data.end_min_index]],
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
    data = venous_waveform_analysis(
        cycle,
        _safe_indexes(ctx.analysis.get("beat_indices")),
        ctx.dt_seconds,
    )
    return _venous_waveform_analysis_plot(writer, cycle, data)


def _venous_waveform_analysis_plot(
    writer: FigureWriter,
    cycle: np.ndarray,
    data: VenousWaveformAnalysis,
) -> Path:
    fig, ax = _plt().subplots(figsize=(6.5, 4.0))
    ax.plot(data.padded_time, data.padded_signal, color="0.85", linewidth=2)
    ax.plot(data.pulse_time, cycle, color="k", linewidth=2)
    _annotated_vline(
        ax,
        data.pulse_time[data.peak_index],
        f"{data.pulse_time[data.peak_index]:.2f} s",
    )
    _annotated_vline(ax, data.period_seconds, f"{data.period_seconds:.2f} s")
    _annotated_hline(ax, cycle[data.peak_index], f"{cycle[data.peak_index]:.1f} mm/s")
    _annotated_vline(
        ax,
        data.pulse_time[data.trough_index],
        f"{data.pulse_time[data.trough_index]:.2f} s",
    )
    _annotated_hline(
        ax,
        cycle[data.trough_index],
        f"{cycle[data.trough_index]:.1f} mm/s",
        color="tab:blue",
        vertical_alignment="bottom",
    )
    ax.scatter(
        [data.pulse_time[data.peak_index]],
        [cycle[data.peak_index]],
        color="tab:blue",
        s=55,
        edgecolor="k",
    )
    ax.scatter(
        [data.pulse_time[data.trough_index]],
        [cycle[data.trough_index]],
        color="lightskyblue",
        s=55,
        edgecolor="k",
    )
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Velocity (mm/s)")
    ax.margins(y=0.18)
    _style_axes(ax)
    return writer.savefig(fig, "VenousWaveformAnalysis_v_vein.png", dpi=180)
