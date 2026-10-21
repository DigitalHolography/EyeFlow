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
from .metrics import PHASE_METRIC_NAMES, _valley_neighborhood_indices, spectrum_valley_height

FIGURE_DPI = 320
FIGURE_WIDTH_INCHES = 8.0 / 2.54
FIGURE_ASPECT_RATIO = 1.618
SPECTRUM_HEIGHT_SCALE = 1.30
SPECTROGRAM_WIDTH_INCHES = 12.0 / 2.54
SPECTROGRAM_ASPECT_RATIO = 2.65
FONT_FAMILY = "Times New Roman"
# Scale text with canvas width to match at the same displayed width.
FONT_SIZE = 10
TICK_FONT_SIZE = 8
SPECTRUM_FONT_SIZE = FONT_SIZE * (2 * FIGURE_WIDTH_INCHES / SPECTROGRAM_WIDTH_INCHES)
SPECTRUM_TICK_FONT_SIZE = TICK_FONT_SIZE * (2 * FIGURE_WIDTH_INCHES / SPECTROGRAM_WIDTH_INCHES)
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
        paths.extend(_export_spectra(output, manager, frequency_table[()]))
    return paths


def _export_spectra(output, manager, frequencies_khz: np.ndarray) -> list[Path]:
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
    peak_sources = [
        output.h5.get(f"{schema.spectrograms_masked_metrics_root}/index_max{half}")
        for half in (1, 2)
    ]
    try:
        peaks = [
            float(np.asarray(source[()]).reshape(())) if source is not None else np.nan
            for source in peak_sources
        ]
    except (TypeError, ValueError):
        peaks = [np.nan, np.nan]
    if any(not np.isfinite(peak) or peak != int(peak) or peak < 0 for peak in peaks):
        Logger.log_warning("Skipping spectra: stored peak indexes must be nonnegative integers.")
        return []
    valley_length = peak_sources[0].attrs.get("valley_neighborhood_length")
    try:
        _valley_neighborhood_indices(valley_length, frequencies_khz.size)
    except ValueError as exc:
        Logger.log_warning(f"Skipping spectra: invalid Hv neighborhood: {exc}")
        return []
    curves = {}
    for label, center, std, style in available:
        if center.shape != frequencies_khz.shape:
            Logger.log_warning(f"Skipping {label.lower()} spectrum: FrequencyBinTable shape does not match.")
            continue
        curve = _spectrum_curve(
            center, std, frequencies_khz=frequencies_khz,
            label=label, linestyle=style, index_max1=peaks[0], index_max2=peaks[1],
            valley_neighborhood_length=valley_length,
        )
        if curve is None:
            Logger.log_warning(
                f"Skipping normalized {label.lower()} spectrum: need valid peak indexes "
                "and a positive finite Hv from the recorded central neighborhood."
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
                partial(_spectra_figure, [curves[label] for label in labels]),
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


def _style_figure(fig, *, font_size: float = FONT_SIZE, tick_font_size: float = TICK_FONT_SIZE) -> None:
    from matplotlib.text import Text

    for text in fig.findobj(match=Text):
        text.set_fontfamily(FONT_FAMILY)
        text.set_fontsize(font_size)
    for ax in fig.axes:
        ax.tick_params(axis="both", labelsize=tick_font_size, width=BOX_LINE_WIDTH)
        ax.tick_params(axis="x", labelrotation=90)
        ax.tick_params(axis="y", labelrotation=0)
        ax.xaxis.get_offset_text().set_fontsize(tick_font_size)
        ax.yaxis.get_offset_text().set_fontsize(tick_font_size)
        ax.xaxis.label.set_fontsize(font_size)
        ax.yaxis.label.set_fontsize(font_size)
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
    # Space ticks within each signed block while retaining the final bin as an
    # upper-bound label. Replace its neighbor if there is too little room.
    blocks = np.split(np.arange(frequencies_khz.size), np.flatnonzero(frequencies_khz < 0)[:1])
    tick_step = max(1, int(np.ceil(frequencies_khz.size / 8)))
    positions = []
    for block in blocks:
        if not block.size:
            continue
        if positions and block[0] - positions[-1] < tick_step:
            positions.pop()
        positions.extend(block[::tick_step])
    last_bin = frequencies_khz.size - 1
    if positions[-1] != last_bin:
        if last_bin - positions[-1] < 0.75 * tick_step:
            positions[-1] = last_bin
        else:
            positions.append(last_bin)
    ax.set_yticks(positions, labels=[str(round(float(value))) for value in frequencies_khz[positions]])
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
class _Spectrum:
    frequencies_khz: np.ndarray
    power: np.ndarray
    std: np.ndarray
    label: str
    linestyle: str
    peak_bins: tuple[int, int]
    valley_bins: np.ndarray


def _spectrum_curve(
    values, std, *, frequencies_khz: np.ndarray, label: str, linestyle: str,
    index_max1: float, index_max2: float, valley_neighborhood_length: int | None = None,
) -> _Spectrum | None:
    values = np.asarray(values, dtype=np.float64)
    std = np.asarray(std, dtype=np.float64)
    frequencies_khz = np.asarray(frequencies_khz)
    if values.ndim != 1 or std.shape != values.shape or frequencies_khz.shape != values.shape:
        raise ValueError("Spectrum, beat STD, and frequency centers must share one frequency axis.")
    middle = values.size // 2
    if (
        values.size < 2
        or any(not np.isfinite(peak) or peak != int(peak) for peak in (index_max1, index_max2))
        or not 0 <= index_max1 < middle <= index_max2 < values.size
    ):
        return None
    divisor = _spectrum_baseline(values, valley_neighborhood_length=valley_neighborhood_length)
    if divisor is None:
        return None
    if not np.isfinite(divisor) or divisor <= 0:
        return None
    return _Spectrum(
        frequencies_khz,
        values / divisor, std / divisor, label, linestyle,
        (int(index_max1), int(index_max2)),
        _valley_neighborhood_indices(valley_neighborhood_length, values.size),
    )


def _spectrum_baseline(values, *, valley_neighborhood_length: int | None = None) -> float | None:
    """Use the same recorded central-neighborhood Hv as A0 over the full spectrum."""
    values = np.asarray(values, dtype=np.float64)
    if values.ndim != 1 or values.size < 2:
        return None
    if np.any(np.isinf(values)):
        return None
    hv = float(spectrum_valley_height(values, valley_neighborhood_length=valley_neighborhood_length))
    return hv if np.isfinite(hv) else None


def _spectra_figure(curves: list[_Spectrum]):
    from matplotlib.patches import Rectangle
    from matplotlib.ticker import FormatStrFormatter, MaxNLocator

    fig, ax = _figure(
        width=2 * FIGURE_WIDTH_INCHES,
        aspect_ratio=2 * FIGURE_ASPECT_RATIO / SPECTRUM_HEIGHT_SCALE,
    )
    # Composite translucent fills on white in one raster layer for identical
    # overlap shading in PNG/EPS; signals, axes, and text remain vector in EPS.
    ax.add_patch(Rectangle(
        (0, 0), 1, 1, transform=ax.transAxes, facecolor="white",
        edgecolor="none", zorder=-2,
    ))
    ax.set_rasterization_zorder(1)
    bins = np.arange(curves[0].power.size)
    valley = curves[0].valley_bins
    ax.axvspan(
        valley[0] - 0.5, valley[-1] + 0.5, facecolor="0.5", alpha=0.1,
        edgecolor="none", zorder=-1,
    )
    for curve in curves:
        lower, upper = curve.power - curve.std, curve.power + curve.std
        _std_envelope(ax, bins, lower, upper)
        ax.plot(
            bins, curve.power, color="black", linestyle=curve.linestyle,
            linewidth=SIGNAL_LINE_WIDTH, label=curve.label.lower(), zorder=3,
            marker="." if curve.power.size == 1 else None,
        )
    ax.axhline(0, color="0.5", linestyle=":", linewidth=BOX_LINE_WIDTH, zorder=2)
    ax.set_xlabel("Frequency (kHz)")
    ax.set_ylabel("Power (a.u.)")
    ax.margins(x=0, y=0.08)
    ax.set_xlim(0, bins[-1])
    positions, _ = frequency_ticks(curves[0].frequencies_khz)
    positions = np.union1d(positions, curves[0].peak_bins)
    ax.set_xticks(positions, labels=[
        str(round(float(curves[0].frequencies_khz[position]))) for position in positions
    ])
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.0f"))
    ax.legend(
        loc="upper center", ncols=len(curves), frameon=False,
        labelspacing=0.25, borderpad=0.1, fontsize=SPECTRUM_FONT_SIZE,
    )
    _style_figure(fig, font_size=SPECTRUM_FONT_SIZE, tick_font_size=SPECTRUM_TICK_FONT_SIZE)
    ax.tick_params(axis="x", labelrotation=0)
    _prune_spectrum_tick_labels(fig, ax, curves[0].peak_bins)
    return fig


def _prune_spectrum_tick_labels(fig, ax, peak_bins: tuple[int, int]) -> None:
    """Keep peak labels first, dropping other ticks that crowd them or each other."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    ticks = {
        int(position): (label.get_text(), label.get_window_extent(renderer))
        for position, label in zip(ax.get_xticks(), ax.get_xticklabels())
    }
    bounds = {min(ticks), max(ticks)}
    priority = list(peak_bins) + [
        position for position in ticks if position not in peak_bins and position not in bounds
    ] + [position for position in ticks if position in bounds and position not in peak_bins]
    keep = []
    gap = 2 * fig.dpi / 72  # Two points of space between rendered labels.
    for position in priority:
        box = ticks[position][1]
        if position not in peak_bins and any(
            box.x0 < ticks[other][1].x1 + gap and box.x1 > ticks[other][1].x0 - gap
            for other in keep
        ):
            continue
        keep.append(position)
    keep.sort()
    ax.set_xticks(keep, labels=[ticks[position][0] for position in keep])


def _std_envelope(ax, x, lower, upper) -> None:
    """Keep signals prominent over soft solid bands whose overlaps blend."""
    ax.fill_between(
        x, lower, upper, facecolor="0.35", alpha=0.45,
        edgecolor="none", zorder=0,
    )


def _gap_metric_figure(values: np.ndarray, phase: np.ndarray, *, metric: str, unit: str):
    fig, ax = _figure()
    ax.plot(phase, values, color="black", linewidth=SIGNAL_LINE_WIDTH)
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 1 / 3, 2 / 3, 1], labels=["0", "1/3", "2/3", "1"])
    ax.set_xlabel(r"Cardiac phase ($t/T$)", math_fontfamily="stix")
    display_unit = "a.u." if unit.replace(" ", "").lower() == "a.u.*khz" else unit
    metric_label = f"${metric[0]}_{{{metric[1:]}}}$"
    ax.set_ylabel(
        metric_label if display_unit == "1" else f"{metric_label} ({display_unit})",
        math_fontfamily="stix",
    )
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
