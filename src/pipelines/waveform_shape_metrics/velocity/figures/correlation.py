"""Arterial/venous correlation PNG exporters."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from calculations.blood_flow_velocity.signal_analysis.spectrum import (
    CorrelationData,
    DelayFitData,
    TransferData,
    paired_spectrum_analysis,
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
from .plotting import _line_plot, _style_axes


def _export_correlation_plots(writer: FigureWriter, ctx: PulseFigureContext) -> list[Path]:
    artery = _display_velocity(_vector(ctx.analysis["retinal_artery_velocity_signal"]))
    vein = -_display_velocity(_vector(ctx.analysis["retinal_vein_velocity_signal"]))
    if artery.size != vein.size or artery.size < 3:
        _log(ctx, "Skipping arterial/venous correlation PNGs; signals are incompatible.")
        return []
    analysis = paired_spectrum_analysis(
        artery,
        vein,
        ctx.dt_seconds,
        _safe_indexes(ctx.analysis.get("beat_indices")),
    )
    corr = analysis.correlation
    spectrum = analysis.transfer
    paths = [
        _msc_plot(writer, corr),
        _correlation_panel_plot(writer, ctx.time, corr),
        _detrended_plot(writer, ctx.time, corr),
        _lags_plot(writer, corr),
        _transfer_plot(writer, "Transfer_function_Velocity_AV_mod.png", spectrum, phase=False),
        _transfer_plot(writer, "Transfer_function_Velocity_AV_phase.png", spectrum, phase=True),
    ]
    if analysis.delay is not None:
        paths.append(_delay_plot(writer, analysis.delay))
    else:
        _log(ctx, "Skipping ArterialVenous_Delay.png; insufficient beats.")
    return paths

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

def _correlation_panel_plot(writer: FigureWriter, time: np.ndarray, corr: CorrelationData) -> Path:
    fig, axes = _plt().subplots(2, 1, figsize=(7.0, 5.2))
    axes[0].plot(time[: corr.first.size], corr.first, color="tab:red", linewidth=2)
    axes[0].plot(time[: corr.second.size], -corr.second, color="tab:blue", linewidth=2)
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
        [(corr.first, "-", "tab:red", "artery"), (-corr.second, "-", "tab:blue", "vein")],
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
    data: DelayFitData,
) -> Path:
    fig, ax = _plt().subplots(figsize=(6.5, 4.0))
    ax.plot(data.time, data.first, color="tab:red", linewidth=2)
    ax.plot(data.time, data.model, color="tab:blue", linestyle="--", linewidth=2)
    ax.plot(data.time, data.second_scaled, color="tab:blue", linewidth=2)
    ax.legend(
        [
            "Artery (normalized)",
            f"Vein model fit (tau_RC = {data.tau_milliseconds:.2f} ms)",
            f"Vein shifted ({data.shift_milliseconds:.2f} ms)",
        ],
        loc="best",
    )
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Arterio-venous decay fit")
    _style_axes(ax)
    return writer.savefig(fig, "ArterialVenous_Delay.png", dpi=180)
