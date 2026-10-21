"""Export shared gap-spectrum figures from their persisted HDF5 products."""

from __future__ import annotations

from dataclasses import dataclass
from functools import partial
from pathlib import Path

import numpy as np

from input_output.schema import EyeFlowOutputPaths
from input_output.writers.eps import EpsArtifactWriter
from input_output.writers.png import PngArtifactWriter
from utils.logger import Logger

from .frequency import frequency_ticks, read_frequency_bin_table
from .metrics import PHASE_METRIC_NAMES, spectrum_valley_height

FIGURE_DPI = 320
FIGURE_WIDTH_INCHES = 8.0 / 2.54
FIGURE_ASPECT_RATIO = 1.618
SPECTROGRAM_WIDTH_INCHES = 12.0 / 2.54
SPECTROGRAM_ASPECT_RATIO = 2.65
FONT_FAMILY = "Times New Roman"
# Shared point sizes across exported figures, suitable for an article column.
FONT_SIZE = 10
TICK_FONT_SIZE = 8
BOX_LINE_WIDTH = 0.4
SIGNAL_LINE_WIDTH = 0.8


def export_spectrogram_figures(output, *, dt_seconds: float | None = None) -> list[Path]:
    """Export shared figures, calibrating heatmap time with the HD frame interval."""
    schema = EyeFlowOutputPaths.active()
    manager = output.manager.for_artifact_namespace("spectrograms")
    paths: list[Path] = []
    frequency_table = read_frequency_bin_table(output)
    for source_path, name, make_figure in (
        (schema.spectrograms_masked, "spectrograms_gaps", _gap_spectrogram_figure),
        (
            schema.spectrograms_masked,
            "spectrograms_gaps_monochrome",
            partial(_gap_spectrogram_figure, monochrome=True),
        ),
        (
            schema.spectrograms_masked_time_aggregated_median,
            "spectrogram_gaps_median",
            _gap_median_figure,
        ),
    ):
        source = output.h5.get(source_path)
        if source is None:
            continue
        if frequency_table is None or source.shape[-1] != frequency_table.size:
            Logger.log_warning(f"Skipping {name}: a matching calibrated FrequencyBinTable is required for kHz axes.")
            continue
        if source_path == schema.spectrograms_masked:
            if dt_seconds is None or not np.isfinite(dt_seconds) or dt_seconds <= 0:
                Logger.log_warning(f"Skipping {name}: a positive finite HD frame interval is required for seconds.")
                continue
            make_figure = partial(make_figure, dt_seconds=dt_seconds)
        values = source[()]
        paths.extend(_save_figure_pair(manager, name, partial(
            make_figure, values, frequencies_khz=frequency_table[()],
        )))
    mask = output.h5.get(schema.segmentation.gaps)
    if mask is not None:
        paths.extend(_save_figure_pair(manager, "gaps_mask", partial(_gaps_mask_figure, mask[()])))
    for name in PHASE_METRIC_NAMES:
        source = output.h5.get(f"{schema.spectrograms_masked_metrics_root}/{name}_t")
        if source is None:
            continue
        paths.extend(_save_figure_pair(manager, f"gaps_{name}", partial(
            _gap_metric_figure, source[()], np.asarray(source.attrs["phase_samples"]),
            metric=name, unit=source.attrs["unit"],
        )))
    if frequency_table is not None:
        paths.extend(_export_positive_spectra(output, manager, frequency_table[()]))
    return paths


def _export_positive_spectra(output, manager, frequencies_khz: np.ndarray) -> list[Path]:
    schema = EyeFlowOutputPaths.active()
    sources = (
        ("Median", schema.spectrograms_masked_time_aggregated_median,
         schema.spectrograms_masked_time_aggregated_std, "-"),
        ("Early", schema.spectrograms_masked_time_aggregated_median_early,
         schema.spectrograms_masked_time_aggregated_std_early, "--"),
        ("Late", schema.spectrograms_masked_time_aggregated_median_late,
         schema.spectrograms_masked_time_aggregated_std_late, ":"),
    )
    available = [
        (label, output.h5.get(center_path)[()], output.h5.get(std_path)[()], style)
        for label, center_path, std_path, style in sources
        if output.h5.get(center_path) is not None and output.h5.get(std_path) is not None
    ]
    if not available:
        return []
    peak_source = output.h5.get(f"{schema.spectrograms_masked_metrics_root}/index_max1")
    try:
        peak = float(np.asarray(peak_source[()]).reshape(())) if peak_source is not None else np.nan
    except (TypeError, ValueError):
        peak = np.nan
    if not np.isfinite(peak) or peak != int(peak) or peak < 0:
        Logger.log_warning("Skipping positive spectra: stored index_max1 must be a nonnegative integer.")
        return []
    curves = {}
    for label, center, std, style in available:
        if center.shape != frequencies_khz.shape:
            Logger.log_warning(f"Skipping {label.lower()} spectrum: FrequencyBinTable shape does not match.")
            continue
        curve = _positive_spectrum_curve(
            center, std, frequencies_khz=frequencies_khz,
            label=label, linestyle=style, index_max1=peak,
        )
        if curve is None:
            Logger.log_warning(
                f"Skipping normalized {label.lower()} spectrum: need a valid positive-bin cutoff "
                "and a positive finite Hv from the central 25% of all bins."
            )
        else:
            curves[label] = curve
    paths = []
    for suffix, labels in (
        ("early_late", ("Early", "Late")),
        ("beat", ("Median",)),
        ("all", ("Median", "Early", "Late")),
    ):
        if all(label in curves for label in labels):
            paths.extend(_save_figure_pair(
                manager, f"gaps_positive_freq_median_phase_{suffix}",
                partial(_positive_spectra_figure, [curves[label] for label in labels]),
            ))
    return paths


def _save_figure_pair(manager, name: str, make_figure) -> list[Path]:
    return [
        PngArtifactWriter(manager).save_figure(
            make_figure(), f"{name}.png", dpi=FIGURE_DPI, bbox_inches=None,
        ),
        EpsArtifactWriter(manager).save_figure(
            make_figure(), f"{name}.eps", dpi=FIGURE_DPI,
        ),
    ]


def _figure(
    *, square: bool = False, width: float = FIGURE_WIDTH_INCHES,
    aspect_ratio: float = FIGURE_ASPECT_RATIO,
):
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    height = width / (1.0 if square else aspect_ratio)
    fig = Figure(figsize=(width, height), constrained_layout=not square)
    FigureCanvasAgg(fig)
    ax = fig.add_axes((0, 0, 1, 1)) if square else fig.subplots()
    return fig, ax


def _style_figure(fig) -> None:
    from matplotlib.text import Text

    for text in fig.findobj(match=Text):
        text.set_fontfamily(FONT_FAMILY)
        text.set_fontsize(FONT_SIZE)
    for ax in fig.axes:
        ax.tick_params(axis="both", labelsize=TICK_FONT_SIZE, width=BOX_LINE_WIDTH)
        ax.tick_params(axis="x", labelrotation=90)
        ax.tick_params(axis="y", labelrotation=0)
        ax.xaxis.get_offset_text().set_fontsize(TICK_FONT_SIZE)
        ax.yaxis.get_offset_text().set_fontsize(TICK_FONT_SIZE)
        ax.xaxis.label.set_fontsize(FONT_SIZE)
        ax.yaxis.label.set_fontsize(FONT_SIZE)
        for spine in ax.spines.values():
            spine.set_color("black")
            spine.set_linewidth(BOX_LINE_WIDTH)


def _gap_spectrogram_figure(
    values: np.ndarray, *, frequencies_khz: np.ndarray, dt_seconds: float,
    monochrome: bool = False,
):
    from matplotlib.colors import PowerNorm
    from matplotlib.ticker import MaxNLocator

    if not np.isfinite(dt_seconds) or dt_seconds <= 0:
        raise ValueError("The HD frame interval must be positive and finite.")
    fig, ax = _figure(width=SPECTROGRAM_WIDTH_INCHES, aspect_ratio=SPECTROGRAM_ASPECT_RATIO)
    # HDF5 is (time, freq); display frequency vertically without changing bin order.
    # Keep frame zero at t=0 and frequency-bin centers at their original indexes.
    image = ax.imshow(
        np.asarray(values).T, origin="lower", aspect="auto",
        interpolation="nearest", cmap="gray_r" if monochrome else "viridis",
        # Slightly darken monochrome midtones without clipping the power range.
        norm=PowerNorm(gamma=0.9) if monochrome else None,
        extent=(-0.5 * dt_seconds, (values.shape[0] - 0.5) * dt_seconds,
                -0.5, values.shape[1] - 0.5),
    )
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Frequency (kHz)", labelpad=6)
    last_time = (values.shape[0] - 1) * dt_seconds
    time_ticks = MaxNLocator(nbins=7, steps=[1, 2, 5, 10]).tick_values(0, last_time)
    ax.set_xticks(time_ticks[(time_ticks >= 0) & (time_ticks <= last_time)])
    # Space ticks within each signed block without crowding the upper edge or
    # the FFT wraparound with an extra endpoint tick.
    blocks = np.split(np.arange(frequencies_khz.size), np.flatnonzero(frequencies_khz < 0)[:1])
    tick_step = max(1, int(np.ceil(frequencies_khz.size / 8)))
    positions = []
    for block in blocks:
        if not block.size:
            continue
        if positions and block[0] - positions[-1] < tick_step:
            positions.pop()
        positions.extend(block[::tick_step])
    ax.set_yticks(positions, labels=[f"{value:.3g}" for value in frequencies_khz[positions]])
    fig.colorbar(image, ax=ax, orientation="vertical", label="Power (a.u.)", pad=0.02)
    _style_figure(fig)
    ax.tick_params(axis="x", labelrotation=0)
    return fig


def _gap_median_figure(values: np.ndarray, *, frequencies_khz: np.ndarray):
    fig, ax = _figure()
    bins = np.arange(values.size)
    ax.plot(bins, values, color="black", linewidth=SIGNAL_LINE_WIDTH)
    ax.axhline(0, color="black", linestyle=":", linewidth=BOX_LINE_WIDTH, zorder=1)
    positions, labels = frequency_ticks(frequencies_khz)
    ax.set_xticks(positions, labels=labels)
    ax.set_xlabel("Frequency (kHz)")
    ax.set_ylabel("Power (a.u.)")
    _style_figure(fig)
    return fig


@dataclass(frozen=True)
class _PositiveSpectrum:
    frequencies_khz: np.ndarray
    power: np.ndarray
    std: np.ndarray
    label: str
    linestyle: str


def _positive_spectrum_curve(
    values, std, *, frequencies_khz: np.ndarray, label: str, linestyle: str,
    index_max1: float,
) -> _PositiveSpectrum | None:
    values = np.asarray(values, dtype=np.float64)
    std = np.asarray(std, dtype=np.float64)
    frequencies_khz = np.asarray(frequencies_khz)
    if values.ndim != 1 or std.shape != values.shape or frequencies_khz.shape != values.shape:
        raise ValueError("Spectrum, beat STD, and frequency centers must share one frequency axis.")
    negative_bins = np.flatnonzero(frequencies_khz < 0)
    stop = int(negative_bins[0]) if negative_bins.size else values.size
    if (
        stop < 1 or not np.isfinite(index_max1) or index_max1 != int(index_max1)
        or not 0 <= index_max1 < stop
    ):
        return None
    divisor = _positive_spectrum_baseline(values)
    if divisor is None:
        return None
    if not np.isfinite(divisor) or divisor <= 0:
        return None
    start = int(index_max1)
    return _PositiveSpectrum(
        frequencies_khz[start:stop],
        values[start:stop] / divisor, std[start:stop] / divisor, label, linestyle,
    )


def _positive_spectrum_baseline(values) -> float | None:
    """Use the same central-quarter Hv as A0, computed before positive-bin cropping."""
    values = np.asarray(values, dtype=np.float64)
    if values.ndim != 1 or values.size < 2:
        return None
    positive = values[:values.size // 2]
    if np.any(np.isinf(positive)):
        return None
    hv = float(spectrum_valley_height(values))
    return hv if np.isfinite(hv) else None


def _positive_spectra_figure(curves: list[_PositiveSpectrum]):
    fig, ax = _figure()
    limits = np.concatenate([np.array([0.0])] + [
        np.r_[curve.power, curve.power - curve.std, curve.power + curve.std] for curve in curves
    ])
    finite = limits[np.isfinite(limits)]
    spacing = max(float(np.ptp(finite)), 1e-6) / 24
    for curve in curves:
        lower, upper = curve.power - curve.std, curve.power + curve.std
        _std_envelope(ax, curve.frequencies_khz, lower, upper, curve.linestyle, spacing)
        ax.plot(
            curve.frequencies_khz, curve.power, color="black", linestyle=curve.linestyle,
            linewidth=SIGNAL_LINE_WIDTH, label=curve.label, zorder=3,
            marker="." if curve.power.size == 1 else None,
        )
    ax.axhline(0, color="black", linestyle=":", linewidth=BOX_LINE_WIDTH, zorder=2)
    ax.set_xlabel("Frequency (kHz)")
    ax.set_ylabel("Power (a.u.)")
    ax.margins(x=0.05, y=0.08)
    first = min(curve.frequencies_khz[0] for curve in curves)
    last = max(curve.frequencies_khz[-1] for curve in curves)
    if last > first:
        padding = (last - first) * 0.05
        ax.set_xlim(first - padding, last + padding)
    if len(curves) > 1:
        ax.legend(loc="upper right", frameon=False, labelspacing=0.25, borderpad=0.1)
    _style_figure(fig)
    ax.tick_params(axis="x", labelrotation=0)
    return fig


def _std_envelope(ax, x, lower, upper, style: str, spacing: float) -> None:
    """Use vector light-gray/dashed/dotted fills consistently in PNG and EPS."""
    from matplotlib.collections import LineCollection
    from matplotlib.path import Path as MplPath

    band = ax.fill_between(
        x, lower, upper, facecolor="0.93" if style == "-" else "none",
        edgecolor="none", zorder=0,
    )
    valid = np.isfinite(lower) & np.isfinite(upper)
    if style != "-" and np.any(valid) and band.get_paths():
        levels = np.arange(
            np.ceil(np.min(lower[valid]) / spacing) * spacing,
            np.max(upper[valid]) + spacing / 2, spacing,
        )
        segments = [[(x[0], level), (x[-1], level)] for level in levels]
        pattern = LineCollection(segments, colors="0.65", linewidths=0.3, linestyles=style, zorder=1)
        pattern.set_clip_path(MplPath.make_compound_path(*band.get_paths()), ax.transData)
        ax.add_collection(pattern, autolim=False)
    for bound in (lower, upper):
        ax.plot(x, bound, color="0.5", linestyle=style, linewidth=BOX_LINE_WIDTH, zorder=2)


def _gap_metric_figure(values: np.ndarray, phase: np.ndarray, *, metric: str, unit: str):
    fig, ax = _figure()
    ax.plot(phase, values, color="black", linewidth=SIGNAL_LINE_WIDTH)
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 1 / 3, 2 / 3, 1], labels=["0", "1/3", "2/3", "1"])
    ax.set_xlabel("Cardiac phase (t/T)")
    display_unit = "a.u." if unit.replace(" ", "").lower() == "a.u.*khz" else unit
    ax.set_ylabel(metric if display_unit == "1" else f"{metric} ({display_unit})")
    _style_figure(fig)
    ax.tick_params(axis="x", labelrotation=0)
    return fig


def _gaps_mask_figure(values: np.ndarray):
    fig, ax = _figure(square=True)
    # Persisted segmentation is lower-left (x, y); imshow rows represent y.
    ax.imshow(
        np.asarray(values, dtype=bool).T, origin="lower", aspect="equal",
        interpolation="nearest", cmap="gray", vmin=0, vmax=1,
    )
    ax.set_xticks([])
    ax.set_yticks([])
    _style_figure(fig)
    fig.set_facecolor("black")
    ax.set_axis_off()
    return fig
