"""Backend-independent Matplotlib figures for file artifacts."""

from __future__ import annotations

from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure


def create_agg_figure(*args, **kwargs) -> Figure:
    """Create a figure with an Agg canvas without loading an interactive backend."""

    figure = Figure(*args, **kwargs)
    FigureCanvasAgg(figure)
    return figure


def close_figure(figure: Figure) -> None:
    """Release artists owned by a backend-independent figure."""

    figure.clear()


__all__ = ["close_figure", "create_agg_figure"]
