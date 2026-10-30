"""Opt-in paired PNG/EPS export for article figures."""

from __future__ import annotations

from pathlib import Path

from plotting.article import ARTICLE_DPI

from .eps import EpsArtifactWriter
from .png import PngArtifactWriter
from .deferred import defer_write


def save_article_figure_pair(
    output, fig, name: str | Path, *, stem: str | None = None, dpi: int = ARTICLE_DPI,
) -> list[Path]:
    """Save one prepared article figure as PNG and EPS, then close it.

    ``name`` is the shared basename (subfolders and existing PNG/EPS extensions
    are accepted). Existing artifact writers own acquisition prefixes and paths.
    Draw once at export DPI using Agg and freeze layout for both formats. Both
    exports retain the complete canvas, including configured outer padding;
    global savefig cropping/transparency defaults cannot change that policy.

    This function owns figure cleanup, including on export failure. Callers
    apply article styling before export. Generic PNG/EPS writers are unchanged.
    """
    from matplotlib import pyplot as plt

    try:
        basename = Path(str(name).replace("\\", "/"))
        if basename.suffix.lower() in {".png", ".eps"}:
            basename = basename.with_suffix("")
        paths = [
            PngArtifactWriter(output, stem).path(f"{basename}.png"),
            EpsArtifactWriter(output, stem).path(f"{basename}.eps"),
        ]
        if not defer_write("article_pair", fig, paths, dpi=dpi):
            render_article_figure_pair(fig, paths, dpi=dpi)
        return paths
    finally:
        plt.close(fig)


def render_article_figure_pair(fig, paths, *, dpi):
    """Apply the identical frozen layout policy in immediate and background saves."""
    from matplotlib.backends.backend_agg import FigureCanvasAgg

    FigureCanvasAgg(fig)
    fig.set_dpi(dpi)
    fig.canvas.draw()
    fig.set_layout_engine("none")
    canvas_bounds = fig.bbox_inches.frozen()
    for path, format_name in zip(paths, ("png", "eps")):
        fig.savefig(
            path, format=format_name, dpi=dpi, bbox_inches=canvas_bounds,
            pad_inches=0, transparent=False,
            facecolor=fig.get_facecolor(), edgecolor=fig.get_edgecolor(),
        )
