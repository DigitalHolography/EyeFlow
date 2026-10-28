"""Matplotlib plotting primitives for velocity PNG figure exporters."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import numpy as np

from plotting.article import (
    ARTICLE_ASPECT_RATIO,
    DIAGNOSTIC_ARTICLE_STYLE,
    apply_article_style,
)
from input_output.writers.png import FigureArtifactWriter as FigureWriter
from plotting.formatting import add_colorbar, add_legend

from .common import _finite_image, _matplotlib, _plt


def _matlab_vessel_colormap(vessel: str):
    if vessel == "vein":
        stops = (
            ([0.0, 0.0, 0.0], 0.0),
            ([0.0, 0.0, 1.0], 1.0 / 3.0),
            ([0.0, 1.0, 1.0], 2.0 / 3.0),
            ([1.0, 1.0, 1.0], 1.0),
        )
    elif vessel == "av":
        stops = (
            ([0.0, 0.0, 0.0], 0.0),
            ([0.5, 0.0, 0.5], 1.0 / 3.0),
            ([1.0, 0.0, 1.0], 2.0 / 3.0),
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
        ax.plot(time[: values.size], values, style, color=color,
                linewidth=DIAGNOSTIC_ARTICLE_STYLE.signal_line_width, label=label)
    if title:
        ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if any(label for _, _, _, label in series):
        add_legend(ax, DIAGNOSTIC_ARTICLE_STYLE.legend_layout)
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
        add_colorbar(fig, im, ax=ax, label=label,
                     layout=DIAGNOSTIC_ARTICLE_STYLE.colorbar_layout)
    apply_article_style(fig)
    return writer.savefig(fig, suffix, dpi=180)

def _mask_background_map(
    writer: FigureWriter,
    background: np.ndarray,
    mask: np.ndarray,
    suffix: str,
    *,
    color: str | None = None,
    cmap: str | None = None,
    values: np.ndarray | None = None,
) -> Path:
    fig, ax = _plt().subplots(figsize=(4.0, 4.0))
    ax.imshow(_finite_image(background), cmap="gray")
    mask = mask.astype(bool)
    if values is not None and cmap is not None:
        colormap = _plt().get_cmap(cmap) if isinstance(cmap, str) else cmap
        overlay = colormap(np.clip(_finite_image(values), 0.0, 1.0))
        overlay[..., 3] = np.where(mask, 0.65, 0.0)
        ax.imshow(overlay)
    else:
        overlay = np.ma.masked_where(~mask, mask)
        ax.imshow(overlay, cmap=_matplotlib().colors.ListedColormap([color or "red"]), alpha=0.45)
    ax.axis("off")
    ax.set_aspect("equal")
    apply_article_style(fig)
    return writer.savefig(fig, suffix, dpi=180)

def _velocity_gradient_values(
    velocity_avg: np.ndarray,
    section_mask: np.ndarray,
    velocity_map: np.ndarray | None = None,
) -> np.ndarray:
    del velocity_map
    velocity = np.abs(np.asarray(velocity_avg, dtype=np.float32))
    section_values = velocity[section_mask]
    finite = np.isfinite(section_values)
    if not np.any(finite):
        return np.zeros_like(velocity, dtype=np.float32)
    vmax = float(np.nanmax(section_values[finite]))
    if vmax <= 0.0:
        return np.zeros_like(velocity, dtype=np.float32)
    return np.clip(velocity / vmax, 0.0, 1.0)

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
    add_colorbar(
        fig,
        _matplotlib().cm.ScalarMappable(norm=norm, cmap=cmap),
        cax=ax, label=label,
        layout=replace(DIAGNOSTIC_ARTICLE_STYLE.colorbar_layout, orientation="horizontal"),
    )
    apply_article_style(fig)
    return writer.savefig(fig, suffix, dpi=150)

def _annotated_vline(ax, x: float, label: str) -> None:
    ax.axvline(x, color="0.25", linestyle="--",
               linewidth=DIAGNOSTIC_ARTICLE_STYLE.reference_line_width, zorder=1)
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
    ax.axhline(y, color=color, linestyle="--",
               linewidth=DIAGNOSTIC_ARTICLE_STYLE.reference_line_width, zorder=1)
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
def _style_axes(ax, *, grid: bool = False) -> None:
    ax.set_box_aspect(1 / ARTICLE_ASPECT_RATIO)
    ax.grid(grid)
    apply_article_style(ax.figure, DIAGNOSTIC_ARTICLE_STYLE, axes=[ax])
