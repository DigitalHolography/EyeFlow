"""Export shared spectral and Doppler CVI figures and movies from HDF5 products."""

from __future__ import annotations

from dataclasses import dataclass, replace
from functools import partial
from pathlib import Path

import numpy as np

from plotting.article import (
    ARTICLE_ASPECT_RATIO as FIGURE_ASPECT_RATIO,
    ARTICLE_DPI as FIGURE_DPI,
    ARTICLE_STYLE,
    ARTICLE_WIDTH_INCHES as FIGURE_WIDTH_INCHES,
    SPECTRUM_ARTICLE_STYLE,
    SPECTROGRAM_ARTICLE_STYLE,
    apply_article_style,
    create_article_figure as _figure,
    prepare_article_fill_layer,
)
from plotting.formatting import TickFormat, add_colorbar, add_legend, apply_tick_format
from input_output.schema import EyeFlowOutputPaths
from input_output.writers.article import save_article_figure_pair
from input_output.writers.png import FigureArtifactWriter
from input_output.writers.avi import AviArtifactWriter
from utils.logger import Logger

from .frequency import frequency_ticks, read_frequency_bin_table
from .metrics import PHASE_METRIC_OUTPUT_NAMES, _valley_neighborhood_indices, spectrum_valley_height
from .qc import valley_repeatability
from .shape_metrics import SHAPE_PHASE_METRIC_NAMES
from .cvi import mean_frequency_occupancy

SPECTRUM_HEIGHT_SCALE = 1.30
SPECTROGRAM_WIDTH_INCHES = 12.0 / 2.54
SPECTROGRAM_ASPECT_RATIO = 2.65
SPECTRUM_STYLE = SPECTRUM_ARTICLE_STYLE.scaled_fonts(
    2 * FIGURE_WIDTH_INCHES / SPECTROGRAM_WIDTH_INCHES,
)
SPECTRUM_FONT_SIZE = SPECTRUM_STYLE.label_size
SPECTRUM_TICK_FONT_SIZE = SPECTRUM_STYLE.tick_size
BOX_LINE_WIDTH = ARTICLE_STYLE.reference_line_width
SIGNAL_LINE_WIDTH = ARTICLE_STYLE.signal_line_width


def export_spectrogram_figures(output, *, dt_seconds: float | None = None) -> list[Path]:
    """Export shared figures, calibrating heatmap time with the HD frame interval."""
    schema = EyeFlowOutputPaths.active()
    manager = output.manager.for_artifact_namespace("spectrograms")
    paths: list[Path] = []
    frequency_table = read_frequency_bin_table(output)
    for source_path, name, make_figure in (
        (
            schema.spectrograms_masked,
            "spectrograms_gaps_monochrome",
            _gap_spectrogram_figure,
        ),
        (
            schema.spectrograms_masked_time_aggregated_median,
            "spectrogram_gaps_median",
            _gap_median_figure,
        ),
        (
            schema.spectrograms_masked_time_aggregated_median_beat,
            "spectrograms_gaps_monochrome_onebeat",
            _gap_onebeat_spectrogram_figure,
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
        elif source_path == schema.spectrograms_masked_time_aggregated_median_beat:
            phase = np.asarray(source.attrs.get("phase_samples", []))
            if (phase.shape != source.shape[:1] or phase.dtype.kind not in "iuf"
                    or not np.all(np.isfinite(phase))
                    or np.any(np.diff(phase) <= 0) or np.any(phase < 0) or np.any(phase >= 1)):
                Logger.log_warning(f"Skipping {name}: valid stored cardiac-phase coordinates are required.")
                continue
            make_figure = partial(make_figure, phase_samples=phase)
        values = source[()]
        paths.extend(_save_figure_pair(manager, name, partial(
            make_figure, values, frequencies_khz=frequency_table[()],
        )))
    mask = output.h5.get(schema.segmentation.gaps)
    if mask is not None:
        paths.extend(_save_figure_pair(manager, "gaps_mask", partial(_gaps_mask_figure, mask[()])))
    for curve_name in (*PHASE_METRIC_OUTPUT_NAMES.values(), *SHAPE_PHASE_METRIC_NAMES):
        source = output.h5.get(f"{schema.spectrograms_masked_metrics_root}/{curve_name}")
        if source is None:
            continue
        artifact_name = curve_name.removesuffix("_t")
        paths.extend(_save_figure_pair(manager, f"gaps_{artifact_name}", partial(
            _gap_metric_figure, source[()], np.asarray(source.attrs["phase_samples"]),
            metric=artifact_name, unit=source.attrs["unit"],
        )))
    if frequency_table is not None:
        paths.extend(_export_spectra(output, manager, frequency_table[()]))
    paths.extend(_export_distribution_figures(output, manager))
    if output.h5.get(schema.spectrograms_masked_qc.root) is not None:
        paths.append(FigureArtifactWriter(manager).save_figure(
            _gap_qc_figure(output, dt_seconds=dt_seconds), "gaps_qc.png", dpi=150,
        ))
    paths.extend(export_cvi_figures(output))
    paths.extend(export_cvi_movies(output))
    return paths


def _gap_qc_figure(output, *, dt_seconds=None):
    """Build a diagnostic panel from persisted QC; export as PNG only."""
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    schema = EyeFlowOutputPaths.active()
    paths = schema.spectrograms_masked_qc

    def read(family, name):
        source = output.h5.get(f"{family}/{name}")
        return None if source is None else np.asarray(source[()])

    def number(value, suffix=""):
        return "unavailable" if value is None or not np.isfinite(value) else f"{float(value):.3g}{suffix}"

    # Pipeline exports run in a worker thread; never create a GUI figure manager.
    fig = Figure(figsize=(12, 10), layout="constrained")
    FigureCanvasAgg(fig)
    axes = fig.subplots(3, 2)
    fig.suptitle("Gap-spectrum quality indicators", fontsize=16)
    support, motion, valley, repeatability, influence, mask_axis = axes.flat
    spatial = output.h5.get(schema.spectrograms_masked)
    counts = np.isfinite(spatial[()]).sum(axis=1)
    calibrated = dt_seconds is not None and np.isfinite(dt_seconds) and dt_seconds > 0
    times = np.arange(counts.size) * (dt_seconds if calibrated else 1)
    time_label = "Time (s)" if calibrated else "Acquisition frame"
    support.step(times, counts, where="mid", color="black")
    support.set(xlabel=time_label, ylabel="Finite frequency bins", title="Data completeness")
    support.set_ylim(bottom=0)
    raw_source = output.h5.get(schema.spectrograms_masked_per_beat)
    beats = raw_source.shape[2] if raw_source is not None else 0
    fraction = read(paths.data_completeness, "ValidSEPhaseFraction")
    support.text(.02, .03,
                 f"Complete beats: {number(beats)}\nValid SE phases: {number(None if fraction is None else 100 * fraction, '%')}",
                 transform=support.transAxes, va="bottom", fontsize=9,
                 bbox={"facecolor": "white", "edgecolor": "none", "alpha": .85})
    score = read(paths.potential_motion_contamination, "RelativeBroadbandChange")
    motion.plot(times, score, color="black", linewidth=1)
    motion.set(xlabel=time_label, ylabel="Relative broadband change", title="Potential motion contamination")
    motion.set_ylim(bottom=0)
    motion.text(.02, .97, "Transient indicator; not specific to saccades", transform=motion.transAxes,
                va="top", fontsize=9)

    relative_source = output.h5.get(f"{paths.valley_repeatability}/RelativeMAD")
    if relative_source is not None and raw_source is not None:
        bins = np.asarray(relative_source.attrs["valley_neighborhood_indices"], dtype=int)
        phase = np.asarray(relative_source.attrs["phase_samples"])
        heights, median, mad, _ = valley_repeatability(raw_source[()], bins)
        valley.plot(phase, heights, color="0.65", linewidth=.8, alpha=.7)
        valley.plot(phase, median, color="black", label="beat median", linewidth=1.5)
        valley.fill_between(phase, median - mad, median + mad, color="0.4", alpha=.25, label="± MAD")
        add_legend(valley, fontsize=9)
        valley.set_ylabel(f"Valley height ({raw_source.attrs.get('unit', 'a.u.')})")
        relative = relative_source[()]
        repeatability.plot(phase, relative, color="black")
        if not np.any(np.isfinite(relative)):
            repeatability.text(.5, .5, "Needs ≥2 valid beats and nonzero median height",
                               transform=repeatability.transAxes, ha="center", fontsize=9)
    else:
        for ax in (valley, repeatability):
            ax.text(.5, .5, "Registered valley spectra unavailable", transform=ax.transAxes, ha="center")
    valley.set(title="Valley heights across beats", xlabel=r"Cardiac phase $t/T$")
    repeatability.set(title="Valley repeatability", xlabel=r"Cardiac phase $t/T$", ylabel="MAD / |median height|")
    repeatability.set_ylim(bottom=0)
    for ax in (valley, repeatability):
        ax.set_xlim(0, 1)

    delta = read(paths.sei_beat_sensitivity, "DeltaSEI")
    influence.axhline(0, color="0.5", linewidth=.8)
    if delta is not None and np.any(np.isfinite(delta)):
        indexes = np.arange(1, delta.size + 1)
        influence.bar(indexes, delta, color="0.4")
        influence.set_xticks(indexes)
        reference = read(schema.spectrograms_masked_metrics_root, "SEI")
        influence.text(.02, .97, f"Full SEI: {number(reference)}", transform=influence.transAxes,
                       va="top", fontsize=9)
    else:
        influence.text(.5, .5, "Needs ≥2 beats and evaluable SEI", transform=influence.transAxes, ha="center")
    influence.set(title="SEI beat sensitivity", xlabel="Omitted beat (chronological)",
                  ylabel="SEI without beat − full SEI")
    mask = output.h5.get(schema.segmentation.gaps)
    if mask is not None:
        mask_axis.imshow(np.asarray(mask[()], dtype=bool).T, origin="lower", cmap="gray",
                         interpolation="nearest", vmin=0, vmax=1)
    else:
        mask_axis.text(.5, .5, "Mask unavailable", transform=mask_axis.transAxes, ha="center")
    coverage = read(paths.spatial_support_quality, "RetainedFieldFraction")
    mask_axis.set_title("Spatial support quality")
    mask_axis.set_xlabel(f"Retained field: {number(None if coverage is None else 100 * coverage, '%')}\n"
                         "White: retained mask; anatomical review required", fontsize=9)
    mask_axis.set_xticks([])
    mask_axis.set_yticks([])
    return fig


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
    return save_article_figure_pair(manager, make_figure(), name, dpi=FIGURE_DPI)


def signed_frequency_heatmap(values, samples, frequencies_khz, *, xlabel, colorbar_label,
                             cmap, norm=None, vmin=None, vmax=None, colorbar_ticks=None):
    """Sort only the display rows; all exported signed bin centers stay distinct."""
    frequencies = np.asarray(frequencies_khz)
    order = np.argsort(frequencies, kind="stable")
    fig, ax = _figure(width=SPECTROGRAM_WIDTH_INCHES,
                                    aspect_ratio=SPECTROGRAM_ASPECT_RATIO)
    image = ax.pcolormesh(samples, frequencies[order], np.asarray(values)[:, order].T,
                         shading="nearest", cmap=cmap, norm=norm, vmin=vmin, vmax=vmax,
                         rasterized=True)
    ax.set(xlabel=xlabel, ylabel="Signed frequency (kHz)")
    colorbar = add_colorbar(fig, image, ax=ax, label=colorbar_label,
                            layout=SPECTROGRAM_ARTICLE_STYLE.colorbar_layout)
    if colorbar_ticks is not None:
        colorbar.set_ticks(colorbar_ticks)
    apply_article_style(fig, SPECTROGRAM_ARTICLE_STYLE)
    return fig


def cvi_heatmap_figure(values, phase, frequencies_khz):
    return signed_frequency_heatmap(
        np.asarray(values) * 100, phase, frequencies_khz, xlabel="Cardiac phase",
        colorbar_label="Doppler CVI proxy (%)", cmap="gray", vmin=0, vmax=100,
        colorbar_ticks=[0, 25, 50, 75, 100],
    )


def cvi_frequency_mean_figure(values, support):
    """Article spatial occupancy map in the frequency movie's native HD frame."""
    from matplotlib.patches import Circle

    fig, ax = _figure(width=12 / 2.54, aspect_ratio=1.)
    image = ax.imshow(np.where(support, values, 0.) * 100, cmap="gray", vmin=0, vmax=100,
                       origin="upper", interpolation="nearest", aspect="equal")
    height, width = np.shape(values)
    ax.add_patch(Circle(((width - 1) / 2, (height - 1) / 2), .4 * width,
                        fill=False, edgecolor="0.6", linewidth=.4))
    ax.set(xlabel="x", ylabel="y",
           xlim=(-.5, width - .5), ylim=(height - .5, -.5))
    colorbar = add_colorbar(fig, image, ax=ax, label="Mean binary occupancy (%)",
                            layout=SPECTROGRAM_ARTICLE_STYLE.colorbar_layout)
    colorbar.set_ticks([0, 25, 50, 75, 100])
    apply_article_style(fig, SPECTROGRAM_ARTICLE_STYLE)
    return fig


def export_cvi_frequency_mean(output):
    root = EyeFlowOutputPaths.active().spectrograms_cvi_root
    mean = output.h5.get(f"{root}/Movies/FrequencyMeanOverBins")
    if mean is None:
        source = output.h5.get(f"{root}/Movies/FrequencyMean")
        if source is None:
            return []
        values = mean_frequency_occupancy(source)
    else:
        values = mean[()]
    figure = cvi_frequency_mean_figure(values, output.h5.get(f"{root}/DiscMask")[()].astype(bool))
    return save_article_figure_pair(
        output.manager.for_artifact_namespace("spectrograms/cvi"), figure,
        "doppler_cvi_binary_frequency_mean",
    )


def export_cvi_figures(output):
    root = EyeFlowOutputPaths.active().spectrograms_cvi_root
    source = output.h5.get(f"{root}/CVI")
    if source is None:
        return []
    figure = cvi_heatmap_figure(source[()], output.h5.get(f"{root}/Phase")[()],
                               output.h5.get(f"{root}/Frequency")[()])
    paths = save_article_figure_pair(
        output.manager.for_artifact_namespace("spectrograms/cvi"), figure, "doppler_cvi_heatmap",
    )
    paths.extend(export_cvi_frequency_mean(output))
    return paths


def export_cvi_movies(output, *, fps=10.):
    """Stream the persisted occupancy maps through the shared AVI writer."""
    from matplotlib import colormaps
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    from matplotlib.patches import Circle

    schema = EyeFlowOutputPaths.active()
    root = schema.spectrograms_cvi_root
    if output.h5.get(f"{root}/Movies/FrequencyMean") is None:
        return []
    frequencies = output.h5.get(f"{root}/Frequency")[()]
    phases = output.h5.get(f"{root}/Phase")[()]
    support = output.h5.get(f"{root}/DiscMask")[()].astype(bool)
    height, width = support.shape
    writer = AviArtifactWriter(output.manager.for_artifact_namespace("spectrograms/cvi"))
    paths = []
    for axis, dataset, order, coordinates, unit in (
        ("phase", "PhaseMean", np.arange(len(phases)), phases, "1"),
        ("frequency", "FrequencyMean", np.argsort(frequencies, kind="stable"), frequencies, "kHz"),
    ):
        source = output.h5.get(f"{root}/Movies/{dataset}")
        if source is None:
            continue
        fig = Figure(figsize=(6.4, 6.4), dpi=100, constrained_layout=True)
        canvas = FigureCanvasAgg(fig)
        ax = fig.subplots()
        cmap = colormaps["gray"].with_extremes(bad="magenta")
        image = ax.imshow(np.zeros(support.shape), cmap=cmap, vmin=0, vmax=1,
                          origin="upper", interpolation="nearest")
        colorbar = fig.colorbar(image, ax=ax, shrink=.8, ticks=[0, .25, .5, .75, 1])
        colorbar.ax.set_yticklabels(["0", "25", "50", "75", "100"])
        colorbar.set_label("Mean binary occupancy (%)")
        ax.add_patch(Circle(((width - 1) / 2, (height - 1) / 2), .4 * width,
                            fill=False, edgecolor="tab:red", linewidth=1))
        ax.set(xlim=(-.5, width - .5), ylim=(height - .5, -.5),
               xlabel="x", ylabel="y")
        fig.text(.5, .01, "Doppler CVI proxy | invalid included pixels: magenta",
                 ha="center", fontsize=8)
        metadata = dict(
            quantity="mean_binary_occupancy", display_range=[0., 1.], fps=fps,
            scientific_dataset=f"/{root}/Movies/{dataset}", threshold_dataset=f"/{root}/Thresholds",
            support_dataset=f"/{root}/DiscMask", frame_axis=axis, coordinate_unit=unit,
            frame_source_indices=order.tolist(), frame_coordinates=coordinates[order].tolist(),
            beat_aggregation="mean_valid_beats_at_each_phase_frequency",
            post_aggregation_threshold="none", image_coordinates="native_HD_yx",
        )
        if axis == "phase":
            selected = np.asarray(source.attrs["frequency_bin_indices"], dtype=int)
            summary = f"Mean binary images over bins {selected[0]} through {selected[-1]} (inclusive)"
            metadata["frequency_bin_indices"] = selected.tolist()
        else:
            summary = "Mean binary images over cardiac phase (all signed bins)"
        with writer.open_mjpeg(f"doppler_cvi_binary_{axis}.avi", width=640, height=640,
                               fps=fps, metadata=metadata, jpeg_quality=95) as video:
            for index in order:
                frame = np.asarray(source[int(index)], dtype=float)
                image.set_data(np.ma.masked_invalid(np.where(support, frame, 0.)))
                label = (f"Cardiac phase {phases[index]:.4f}" if axis == "phase" else
                         f"Source bin {index}: {frequencies[index]:.5g} kHz")
                ax.set_title(f"{label}\n{summary}", fontsize=10)
                canvas.draw()
                video.write_frame(np.asarray(canvas.buffer_rgba())[:, :, :3])
        paths.append(video.path)
        fig.clear()
    return paths


def _gap_onebeat_spectrogram_figure(values, *, phase_samples, frequencies_khz):
    """Display the persisted phase/frequency beat median with its own power range."""
    from matplotlib.colors import PowerNorm

    fig = signed_frequency_heatmap(
        values, phase_samples, frequencies_khz, xlabel="Cardiac phase",
        colorbar_label="Power (a.u.)", cmap="gray", norm=PowerNorm(gamma=0.9),
    )
    fig.axes[0].set_xlim(0, 1)
    return fig


def _gap_spectrogram_figure(values: np.ndarray, *, frequencies_khz: np.ndarray, dt_seconds: float):
    from matplotlib.colors import PowerNorm
    from matplotlib.ticker import MaxNLocator

    if not np.isfinite(dt_seconds) or dt_seconds <= 0:
        raise ValueError("The HD frame interval must be positive and finite.")
    fig = signed_frequency_heatmap(
        values, np.arange(values.shape[0]) * dt_seconds, frequencies_khz,
        xlabel="Time (s)", colorbar_label="Power (a.u.)", cmap="gray", norm=PowerNorm(gamma=0.9),
    )
    ax = fig.axes[0]
    last_time = (values.shape[0] - 1) * dt_seconds
    ticks = MaxNLocator(nbins=7, steps=[1, 2, 5, 10]).tick_values(0, last_time)
    ax.set_xticks(ticks[(ticks >= 0) & (ticks <= last_time)])
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
    apply_article_style(fig, replace(ARTICLE_STYLE, x_tick_rotation=90))
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
    fig, ax = _figure(
        width=2 * FIGURE_WIDTH_INCHES,
        aspect_ratio=2 * FIGURE_ASPECT_RATIO / SPECTRUM_HEIGHT_SCALE,
        signal_padding=(False, True),
    )
    # Composite translucent fills on white in one raster layer for identical
    # overlap shading in PNG/EPS; signals, axes, and text remain vector in EPS.
    prepare_article_fill_layer(ax)
    bins = np.arange(curves[0].power.size)
    valley = curves[0].valley_bins
    ax.axvspan(
        valley[0] - 0.5, valley[-1] + 0.5, **SPECTRUM_STYLE.region_kwargs(), zorder=-1,
    )
    for curve in curves:
        lower, upper = curve.power - curve.std, curve.power + curve.std
        _std_envelope(ax, bins, lower, upper)
        ax.plot(
            bins, curve.power, color="black", linestyle=curve.linestyle,
            linewidth=SPECTRUM_STYLE.signal_line_width, label=curve.label.lower(), zorder=3,
            marker="." if curve.power.size == 1 else None,
        )
    ax.axhline(0, color="0.5", linestyle=":", linewidth=SPECTRUM_STYLE.reference_line_width, zorder=2)
    ax.set_xlabel("Frequency (kHz)")
    ax.set_ylabel("Power (a.u.)")
    ax.margins(y=0.08)
    ax.set_xlim(0, bins[-1])
    positions, _ = frequency_ticks(curves[0].frequencies_khz)
    positions = np.union1d(positions, curves[0].peak_bins)
    apply_tick_format(ax.xaxis, TickFormat(decimals=0),
                      positions=positions, values=curves[0].frequencies_khz[positions])
    add_legend(ax, SPECTRUM_STYLE.legend_layout, ncols=len(curves))
    apply_article_style(fig, SPECTRUM_STYLE)
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
        x, lower, upper, **SPECTRUM_STYLE.envelope_kwargs(), zorder=0,
    )


def _gap_metric_figure(values: np.ndarray, phase: np.ndarray, *, metric: str, unit: str):
    fig, ax = _figure()
    ax.plot(phase, values, color="black", linewidth=SIGNAL_LINE_WIDTH)
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 1 / 3, 2 / 3, 1], labels=["0", "1/3", "2/3", "1"])
    ax.set_xlabel(r"Cardiac phase ($t/T$)")
    display_unit = "a.u." if unit.replace(" ", "").lower() == "a.u.*khz" else unit
    if metric == "SpectralEntropy":
        metric_label = "$H$"
    elif "_" in metric:
        base, subscript = metric.split("_", 1)
        metric_label = f"${base}_{{{subscript}}}$"
    else:
        metric_label = f"${metric}$" if metric in {"A", "SE"} else f"${metric[0]}_{{{metric[1:]}}}$"
    ax.set_ylabel(
        metric_label if display_unit == "1" else f"{metric_label} ({display_unit})",
    )
    apply_article_style(fig)
    ax.tick_params(axis="x", labelrotation=0)
    return fig


def _export_distribution_figures(output, manager) -> list[Path]:
    root = EyeFlowOutputPaths.active().spectrograms_masked_metrics_root
    paths = []
    percent = output.h5.get(f"{root}/PowerQuantilePercent")
    quantiles = [(label, output.h5.get(f"{root}/{name}"), style) for label, name, style in (
        ("Whole cycle", "PowerQuantileFrequency", "-"),
        ("Early", "PowerQuantileFrequency_early", "--"),
        ("Late", "PowerQuantileFrequency_late", ":"),
    )]
    if percent is not None and percent.ndim == 1 and any(ds is not None for _, ds, _ in quantiles):
        curves = [(label, ds[()], style) for label, ds, style in quantiles
                  if ds is not None and ds.shape == percent.shape and ds.attrs.get("unit") == "kHz"]
        if curves:
            paths.extend(_save_figure_pair(manager, "gaps_power_quantiles", partial(
                _power_quantile_figure, percent[()], curves,
            )))
    sources = [output.h5.get(f"{root}/{name}") for name in (
        "WassersteinFrequencyMagnitude", "WassersteinEarlyCDF", "WassersteinLateCDF", "WassersteinEarlyLate",
    )]
    if all(ds is not None for ds in sources):
        frequency, early, late, distance = sources
        if frequency.ndim == 1 and early.shape == late.shape == frequency.shape and frequency.attrs.get("unit") == "kHz":
            paths.extend(_save_figure_pair(manager, "gaps_wasserstein_early_late", partial(
                _spectral_cdf_figure, frequency[()], early[()], late[()], float(distance[()]),
            )))
    return paths


def _power_quantile_figure(percent, curves):
    fig, ax = _figure(width=2 * FIGURE_WIDTH_INCHES, aspect_ratio=2)
    for label, values, style in curves:
        ax.step(percent, values, where="post", color="black", linestyle=style,
                linewidth=SIGNAL_LINE_WIDTH, label=label)
    ax.set(xlim=(0, 100), xlabel="Cumulative spectral power (%)", ylabel=r"$F_q$ (kHz)")
    add_legend(ax, SPECTRUM_STYLE.legend_layout, ncols=len(curves))
    apply_article_style(fig)
    return fig


def _spectral_cdf_figure(frequency, early, late, distance):
    fig, ax = _figure(width=2 * FIGURE_WIDTH_INCHES, aspect_ratio=2)
    if frequency[0] > 0:
        frequency, early, late = np.r_[0., frequency], np.r_[0., early], np.r_[0., late]
    ax.fill_between(frequency, early, late, step="post", color="0.9", linewidth=0)
    for values, style, label in ((early, "-", "Early"), (late, "--", "Late")):
        ax.step(frequency, values, where="post", color="black", linestyle=style,
                linewidth=SIGNAL_LINE_WIDTH, label=label)
    ax.set(xlim=(0, frequency[-1]), ylim=(0, 1), xlabel=r"$|f|$ (kHz)", ylabel="Cumulative power fraction")
    ax.text(.98, .08, f"$W_1$ = {distance:.3g} kHz", transform=ax.transAxes, ha="right", fontsize=9)
    add_legend(ax, SPECTRUM_STYLE.legend_layout, loc="lower right", bbox_to_anchor=(1, .15))
    apply_article_style(fig)
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
    apply_article_style(fig)
    fig.set_facecolor("black")
    ax.set_axis_off()
    return fig
