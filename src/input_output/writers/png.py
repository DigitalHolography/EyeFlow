"""Write figure artifacts for EyeFlow runs."""

from pathlib import Path

import numpy as np
from skimage.io import imsave

from input_output.output_manager import OutputType, artifact_path, output_stem

REPORT_IMAGES_STATE = "report_image_paths"
ReportImagePaths = dict[tuple[str, str], Path]


class FigureArtifactWriter:
    """Write stem-prefixed PNG figures for one output namespace."""

    def __init__(self, output, stem: str | None = None) -> None:
        self.output = output
        self.stem = str(stem) if stem else output_stem(output)
        self.artifacts: dict[tuple[str, str], Path] = {}

    def register_artifact(self, kind: str, vessel: str, path: Path) -> None:
        """Record the actual PNG path for a downstream consumer."""
        self.artifacts[(kind, vessel)] = Path(path)

    def path(self, suffix: str, *, subfolder: str | None = None) -> Path:
        return artifact_path(
            self.output, OutputType.PNG, suffix, stem=self.stem, subfolder=subfolder
        )

    def save_array(self, image, suffix: str) -> Path:
        filename = f"{self.stem}_{suffix}"
        return self.output.write_png(image, filename)

    def save_figure(
        self,
        fig,
        suffix: str,
        *,
        dpi: int = 150,
        bbox_inches="tight",
        subfolder: str | None = None,
    ) -> Path:
        return write_png_figure(
            self.path(suffix, subfolder=subfolder),
            fig,
            dpi=dpi,
            bbox_inches=bbox_inches,
            close=True,
        )

    def save_image(self, image, suffix: str) -> Path:
        return self.save_array(image, suffix)

    def savefig(
        self,
        fig,
        suffix: str,
        *,
        dpi: int = 150,
        bbox_inches="tight",
        subfolder: str | None = None,
    ) -> Path:
        return self.save_figure(
            fig,
            suffix,
            dpi=dpi,
            bbox_inches=bbox_inches,
            subfolder=subfolder,
        )


def write_png_figure(
    path: str | Path,
    fig,
    *,
    dpi: int = 150,
    bbox_inches="tight",
    close: bool = False,
) -> Path:
    """Save a Matplotlib figure to PNG at an explicitly chosen path."""
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    try:
        fig.savefig(target, dpi=dpi, bbox_inches=bbox_inches)
    finally:
        if close:
            _close_figure(fig)
    return target


def write_png_file(path: str | Path, image) -> Path:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    imsave(target, _uint8_image(image), check_contrast=False)
    return target


def _uint8_image(image) -> np.ndarray:
    array = np.asarray(image)
    if array.dtype == np.uint8:
        return array
    if array.dtype == bool:
        return (array.astype(np.uint8) * 255)
    if np.issubdtype(array.dtype, np.floating):
        return _normalize_float(array)
    clipped = np.clip(array, 0, 255)
    return clipped.astype(np.uint8, copy=False)


def _normalize_float(array: np.ndarray) -> np.ndarray:
    finite = np.isfinite(array)
    if not np.any(finite):
        return np.zeros(array.shape, dtype=np.uint8)
    values = array[finite]
    min_value = np.min(values)
    span = np.max(values) - min_value
    if span <= 0:
        return np.where(finite, 255, 0).astype(np.uint8)
    scaled = np.zeros(array.shape, dtype=np.float32)
    scaled[finite] = (array[finite] - min_value) / span
    return np.rint(scaled * 255).astype(np.uint8)


# Backward-compatible name for callers that only use the PNG methods.
PngArtifactWriter = FigureArtifactWriter


def _close_figure(fig) -> None:
    import matplotlib.pyplot as plt

    plt.close(fig)
