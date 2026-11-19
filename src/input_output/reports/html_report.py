"""Write self-contained HTML reports and encode their image assets."""

from __future__ import annotations

import base64
from io import BytesIO
from pathlib import Path

import numpy as np
from PIL import Image

from input_output.artifact_ownership import register_artifact_path
from plotting.headless import close_figure


def figure_png_data_uri(figure, *, dpi: int = 180) -> str:
    """Render and close a Matplotlib figure as an embedded PNG URI."""
    buffer = BytesIO()
    try:
        figure.savefig(
            buffer,
            format="png",
            dpi=dpi,
            bbox_inches="tight",
            facecolor=figure.get_facecolor(),
            edgecolor=figure.get_edgecolor(),
        )
    finally:
        close_figure(figure)
    return _png_data_uri(buffer.getvalue())


def image_png_data_uri(image: np.ndarray) -> str:
    """Encode an RGB-like array as an embedded PNG URI."""
    array = np.asarray(image)
    if array.dtype != np.uint8:
        array = np.clip(array, 0, 255).astype(np.uint8)
    buffer = BytesIO()
    Image.fromarray(array).save(buffer, format="PNG")
    return _png_data_uri(buffer.getvalue())


def write_html_report(path: Path, document: str) -> Path:
    """Write one complete UTF-8 HTML document as an owned artifact."""
    target = Path(path)
    register_artifact_path(target, writing=True)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(document, encoding="utf-8", newline="\n")
    return target


def _png_data_uri(payload: bytes) -> str:
    encoded = base64.b64encode(payload).decode("ascii")
    return f"data:image/png;base64,{encoded}"


__all__ = [
    "figure_png_data_uri",
    "image_png_data_uri",
    "write_html_report",
]
