"""Spectral-analysis PNG exporters."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .spectrum import (
    SpectrumData,
    SyntheticSpectrumData,
    spectrum_signal_analysis,
    synthetic_spectrum_from_signals,
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
from .plotting import _style_axes


def _export_spectral_plots(writer: FigureWriter, ctx: PulseFigureContext) -> list[Path]:
    paths: list[Path] = []
    artery = _display_velocity(_vector(ctx.analysis["retinal_artery_velocity_signal"]))
    vein = _display_velocity(_vector(ctx.analysis["retinal_vein_velocity_signal"]))
    for prefix, vessel_name, values, color in (
        ("Arterial", "artery", artery, "tab:red"),
        ("Venous", "vein", vein, "tab:blue"),
    ):
        spectrum = spectrum_signal_analysis(values, ctx.dt_seconds)
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
    synthetic = synthetic_spectrum_from_signals(vein, artery, peaks, ctx.dt_seconds)
    if synthetic is not None:
        paths.append(_synthetic_spectral_plot(writer, synthetic))
    else:
        _log(ctx, "Skipping syntheticSpectralAnalysis.png; insufficient beats.")
    return paths

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
    spectrum: SyntheticSpectrumData,
) -> Path:
    fig, axes = _plt().subplots(2, 1, figsize=(7.0, 5.0), sharex=True)
    axes[0].plot(spectrum.frequencies, spectrum.magnitude, color="k", linewidth=2)
    axes[0].scatter(
        spectrum.frequencies[spectrum.peak_indexes],
        spectrum.magnitude[spectrum.peak_indexes],
        s=25,
    )
    axes[0].set_ylabel("Magnitude |Y|")
    axes[1].plot(spectrum.frequencies, spectrum.phase / np.pi, color="k", linewidth=2)
    axes[1].scatter(
        spectrum.frequencies[spectrum.peak_indexes],
        spectrum.phase[spectrum.peak_indexes] / np.pi,
        s=25,
    )
    axes[1].set_ylabel("Phase / pi")
    axes[1].set_xlabel("Frequency (Hz)")
    for ax in axes:
        ax.set_xlim(0, 10)
        _style_axes(ax)
    axes[1].set_ylim(-1, 1)
    return writer.savefig(fig, "syntheticSpectralAnalysis.png", dpi=180)
