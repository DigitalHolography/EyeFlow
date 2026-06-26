"""Write PNG artifacts for EyeFlow runs."""

from pathlib import Path

import numpy as np
from skimage.io import imsave


class PngArtifactWriter:
    """Write stem-prefixed PNG arrays and matplotlib figures for one output namespace."""

    def __init__(self, output, stem: str | None = None) -> None:
        self.output = output
        self.stem = str(stem) if stem else _output_stem(output)

    def path(self, suffix: str) -> Path:
        filename = f"{self.stem}_{suffix}"
        path = self.output.path_for(_png_output_type(), filename)
        path.parent.mkdir(parents=True, exist_ok=True)
        return path

    def save_array(self, image, suffix: str) -> Path:
        filename = f"{self.stem}_{suffix}"
        return self.output.write_png(image, filename)

    def save_figure(self, fig, suffix: str, *, dpi: int = 150) -> Path:
        path = self.path(suffix)
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        _close_figure(fig)
        return path

    def save_image(self, image, suffix: str) -> Path:
        return self.save_array(image, suffix)

    def savefig(self, fig, suffix: str, *, dpi: int = 150) -> Path:
        return self.save_figure(fig, suffix, dpi=dpi)


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


def _output_stem(output) -> str:
    manager = getattr(output, "manager", None)
    layout = getattr(manager, "layout", None)
    stem = getattr(layout, "stem", None)
    return str(stem or "eyeflow")


def _close_figure(fig) -> None:
    import matplotlib.pyplot as plt

    plt.close(fig)


def _png_output_type():
    from input_output.output_manager import OutputType

    return OutputType.PNG
