"""Read report images and write finished PDF pages."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.figure import Figure
from PIL import Image

from input_output.artifact_ownership import register_artifact_path


def load_first_image(paths: Iterable[Path]) -> np.ndarray | None:
    """Read the first usable RGB image from candidate paths."""
    for path in paths:
        try:
            with Image.open(path) as image:
                return np.asarray(image.convert("RGB"))
        except OSError:
            continue
    return None


def write_pdf_report(path: Path, pages: Iterable[Figure]) -> Path:
    """Save and close completed figures in one PDF file."""
    path = Path(path)
    register_artifact_path(path, writing=True)
    path.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(path) as pdf:
        for page in pages:
            try:
                pdf.savefig(page, dpi=300)
            finally:
                plt.close(page)
    return path
