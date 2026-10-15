"""Export shared gap-spectrum figures from their persisted HDF5 products."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from input_output.schema import EyeFlowOutputPaths
from input_output.writers.eps import EpsArtifactWriter
from input_output.writers.png import PngArtifactWriter

FIGURE_DPI = 320
FIGURE_WIDTH_INCHES = 8.0 / 2.54
FIGURE_ASPECT_RATIO = 1.618
FONT_FAMILY = "Times New Roman"
FONT_SIZE = 8
TICK_FONT_SIZE = 5
BOX_LINE_WIDTH = 0.4
SIGNAL_LINE_WIDTH = 0.8


def export_spectrogram_figures(output) -> list[Path]:
    """Write each available product once under the shared spectrograms namespace."""
    schema = EyeFlowOutputPaths.active()
    manager = output.manager.for_artifact_namespace("spectrograms")
    paths: list[Path] = []
    for source_path, name, make_figure in (
        (schema.spectrograms_masked, "spectrograms_gaps", _gap_spectrogram_figure),
        (
            schema.spectrograms_masked_time_aggregated_median,
            "spectrogram_gaps_median",
            _gap_median_figure,
        ),
        (schema.segmentation.gaps, "gaps_mask", _gaps_mask_figure),
    ):
        source = output.h5.get(source_path)
        if source is None:
            continue
        values = source[()]
        paths.append(PngArtifactWriter(manager).save_figure(
            make_figure(values), f"{name}.png", dpi=FIGURE_DPI, bbox_inches=None,
        ))
        paths.append(EpsArtifactWriter(manager).save_figure(
            make_figure(values), f"{name}.eps", dpi=FIGURE_DPI,
        ))
    return paths


def _figure(*, square: bool = False):
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    height = FIGURE_WIDTH_INCHES / (1.0 if square else FIGURE_ASPECT_RATIO)
    fig = Figure(figsize=(FIGURE_WIDTH_INCHES, height), constrained_layout=not square)
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


def _gap_spectrogram_figure(values: np.ndarray):
    fig, ax = _figure()
    # HDF5 is (time, freq); display frequency vertically without changing bin order.
    image = ax.imshow(
        np.asarray(values).T, origin="lower", aspect="auto",
        interpolation="nearest", cmap="viridis",
    )
    ax.set_xlabel("Time (frames)")
    ax.set_ylabel("Frequency band")
    ax.set_xticks(np.arange(0, values.shape[0], 50))
    ax.set_yticks(np.arange(0, values.shape[1], 4))
    fig.colorbar(image, ax=ax, orientation="vertical", label="Power (a.u.)")
    _style_figure(fig)
    ax.tick_params(axis="x", labelrotation=0)
    return fig


def _gap_median_figure(values: np.ndarray):
    fig, ax = _figure()
    bins = np.arange(values.size)
    ax.plot(bins, values, color="black", linewidth=SIGNAL_LINE_WIDTH)
    ax.set_xticks(bins[::2])
    ax.set_xlabel("Frequency band")
    ax.set_ylabel("Power (a.u.)")
    _style_figure(fig)
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
