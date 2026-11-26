"""Write figure artifacts for EyeFlow runs."""

from pathlib import Path

import numpy as np
from skimage.io import imsave

from .artifact_names import (
    acquisition_stem, labeled_artifact_path, prefixed_artifact_path,
    writer_artifact_path,
)
from .deferred import defer_write, save_figure_file

REPORT_IMAGES_STATE = "report_image_paths"
ReportImagePaths = dict[tuple[str, str], Path]


class FigureArtifactWriter:
    """Write stem-prefixed figure and image artifacts for one output namespace."""

    def __init__(self, output, stem: str | None = None) -> None:
        self.output = output
        self.stem = str(stem) if stem else acquisition_stem(output)
        self.acquisition_stem = acquisition_stem(output, Path(self.stem).name)
        self.artifacts: ReportImagePaths = {}

    def register_artifact(self, kind: str, vessel: str, path: Path) -> None:
        """Record the actual PNG path for downstream report assembly."""
        self.artifacts[(kind, vessel)] = Path(path)

    def path(self, suffix: str, *, subfolder: str | None = None) -> Path:
        return writer_artifact_path(
            self.output, _png_output_type(), suffix, label=self.stem, subfolder=subfolder
        )

    def save_image(self, image, suffix: str) -> Path:
        filename = str(labeled_artifact_path(suffix, self.acquisition_stem, self.stem))
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
            self.path(suffix, subfolder=subfolder), fig,
            dpi=dpi, bbox_inches=bbox_inches, close=True,
        )

    def save_figure_pair(
        self,
        fig,
        suffix: str,
        *,
        dpi: int = 150,
        bbox_inches="tight",
        subfolder: str | None = None,
    ) -> list[Path]:
        """Save one figure as PNG and EPS with the diagnostic crop policy."""

        from .eps import EpsArtifactWriter

        name = Path(subfolder) / suffix if subfolder else Path(suffix)
        png_path = self.path(str(name.with_suffix(".png")))
        eps_path = EpsArtifactWriter(self.output, self.stem).path(
            str(name.with_suffix(".eps"))
        )
        paths = [png_path, eps_path]
        missing = [path for path in paths if not path.exists()]
        try:
            if missing and not defer_write(
                "figure_pair", fig, missing, dpi=dpi, bbox_inches=bbox_inches,
            ):
                render_figure_pair(
                    fig, missing, dpi=dpi, bbox_inches=bbox_inches,
                )
            return paths
        finally:
            _close_figure(fig)

    def save_image_pair(
        self,
        image,
        suffix: str,
        *,
        dpi: int = 150,
    ) -> list[Path]:
        """Save a raster source as PNG plus an EPS image in matching folders."""

        from matplotlib.backends.backend_agg import FigureCanvasAgg
        from matplotlib.figure import Figure

        from .eps import EpsArtifactWriter, write_eps_file

        png_path = self.path(suffix)
        eps_suffix = str(Path(suffix).with_suffix(".eps"))
        eps_path = EpsArtifactWriter(self.output, self.stem).path(eps_suffix)
        array = np.asarray(image)
        if not png_path.exists():
            write_png_file(png_path, array)
        if not eps_path.exists():
            height, width = array.shape[:2]
            fig = Figure(
                figsize=(max(width, 1) / dpi, max(height, 1) / dpi),
                dpi=dpi,
            )
            FigureCanvasAgg(fig)
            ax = fig.add_axes((0, 0, 1, 1))
            ax.imshow(array, aspect="auto")
            ax.set_axis_off()
            ax.set_xlim(-0.5, width - 0.5)
            ax.set_ylim(height - 0.5, -0.5)
            write_eps_file(eps_path, fig, dpi=dpi, close=True)
        return [png_path, eps_path]


def write_png_figure(
    path: str | Path, fig, *, dpi: int = 150,
    bbox_inches="tight", close: bool = False,
) -> Path:
    """Save a figure to an explicit PNG path and optionally close it."""
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    try:
        save_figure_file(target, fig, dpi=dpi, bbox_inches=bbox_inches)
    finally:
        if close:
            _close_figure(fig)
    return target


def render_figure_pair(fig, paths, **kwargs) -> None:
    """Render any missing PNG/EPS diagnostic formats from one figure."""

    for path in paths:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        format_name = target.suffix.lower().lstrip(".")
        if format_name not in {"png", "eps"}:
            raise ValueError(f"Unsupported figure format: {target}")
        fig.savefig(target, format=format_name, **kwargs)


def write_png_file(path: str | Path, image, *, stem: str | None = None) -> Path:
    target = prefixed_artifact_path(path, stem) if stem is not None else Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    if not defer_write("array", np.asarray(image), [target]):
        imsave(target, _uint8_image(image), check_contrast=False)
    return target


def write_cv_png_file(path: str | Path, image) -> Path:
    """Write OpenCV diagnostic images with their existing channel and byte policy."""
    import cv2

    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    from input_output.artifact_ownership import register_artifact_path

    register_artifact_path(target, writing=True)
    if not cv2.imwrite(str(target), image):
        raise OSError(f"Could not write PNG image: {target}")
    return target


def _uint8_image(image) -> np.ndarray:
    array = np.asarray(image)
    if array.dtype == np.uint8:
        return array
    if array.dtype == bool:
        return array.astype(np.uint8) * 255
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


def _png_output_type():
    from input_output.output_manager import OutputType

    return OutputType.PNG
