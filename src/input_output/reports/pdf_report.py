"""Render the EyeFlow A4 PDF report from measurements and registered images."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.figure import Figure
from PIL import Image

from .pdf_data import extract_parameters_from_h5, format_parameters_for_display


def generate_a4_report(
    output_h5_path: Path,
    output_path: Path,
    folder_name: str,
    report_images: Mapping[tuple[str, str], Path] | None = None,
    hd_png_dir: Path | None = None,
    mask_dir: Path | None = None,
) -> Path:
    """Create the established single-page A4 report without changing its layout."""
    pdf_path = Path(output_path)
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    parameters = extract_parameters_from_h5(Path(output_h5_path))
    images = report_images or {}
    with PdfPages(pdf_path) as pdf:
        figure = _report_page(
            folder_name,
            parameters,
            images,
            Path(hd_png_dir) if hd_png_dir is not None else None,
            Path(mask_dir) if mask_dir is not None else None,
        )
        pdf.savefig(figure, dpi=300)
    return pdf_path


def _report_page(folder_name, parameters, images, hd_png_dir, mask_dir) -> Figure:
    figure = Figure(figsize=(21 / 2.54, 29.7 / 2.54), dpi=150)
    figure.patch.set_facecolor("white")
    width, height = figure.get_size_inches()
    grid = figure.add_gridspec(
        16,
        2,
        left=1.5 / 2.54 / width,
        right=1 - 1.5 / 2.54 / width,
        bottom=0.5 / 2.54 / height,
        top=1 - 1.5 / 2.54 / height,
        hspace=0.05,
        wspace=0.1,
    )
    title = figure.add_subplot(grid[0, :])
    title.axis("off")
    title.text(
        0.5, 0.5, folder_name, fontsize=16, fontweight="bold",
        ha="center", va="center", transform=title.transAxes, wrap=True,
    )
    for column, vessel in enumerate(("artery", "vein")):
        for rows, image, options in (
            (
                slice(1, 5),
                _try_load_vessel_image(images, mask_dir, hd_png_dir, folder_name, vessel),
                {},
            ),
            (slice(5, 8), _try_load_ri_image(images, vessel),
             {"zoom": 1.02, "right_pad_fraction": 0.04}),
            (slice(8, 11), _try_load_systole_image(images, vessel),
             {"zoom": 1.02, "right_pad_fraction": 0.04}),
        ):
            axis = figure.add_subplot(grid[rows, column])
            axis.axis("off")
            if image is not None:
                _show_report_image(axis, image, **options)
    parameter_axis = figure.add_subplot(grid[11:16, :])
    parameter_axis.axis("off")
    _add_parameters_section(parameter_axis, parameters)
    return figure


def _show_report_image(
    axis,
    image: np.ndarray,
    *,
    zoom: float = 1.04,
    right_pad_fraction: float = 0.0,
) -> None:
    """Draw an image without distorting its pixel aspect ratio."""
    height, width = np.asarray(image).shape[:2]
    axis.imshow(image, aspect="equal")
    axis.set_anchor("C")
    axis.set_adjustable("box")
    axis.margins(0)
    if zoom > 1.0 and width > 0 and height > 0:
        x_center = (width - 1) / 2
        y_center = (height - 1) / 2
        half_width = width / (2 * zoom)
        half_height = height / (2 * zoom)
        right_pad = width * max(right_pad_fraction, 0.0)
        axis.set_xlim(x_center - half_width, x_center + half_width + right_pad)
        axis.set_ylim(y_center + half_height, y_center - half_height)


def _add_parameters_section(axis, parameters: dict[str, Any]) -> None:
    texts = format_parameters_for_display(parameters)
    if not texts:
        axis.text(0.02, 0.5, "No parameters extracted from H5 files", fontsize=10, va="center")
        return
    axis.text(0.02, 0.98, "Computed Parameters:", fontsize=14, fontweight="bold", va="top")
    if len(texts) > 6:
        split = len(texts) // 2
        columns = ((0.02, texts[:split]), (0.52, texts[split:]))
        step = 0.055
    else:
        columns = ((0.02, texts),)
        step = 0.08
    for x, column in columns:
        for index, text in enumerate(column):
            y = 0.86 - index * step
            if y > 0.05:
                axis.text(x, y, text, fontsize=9, va="top")


def _try_load_vessel_image(report_images, mask_dir, hd_png_dir, folder_name, vessel_type):
    """Prefer the registered vessel map, then legacy external HD images."""
    paths: list[Path] = []
    exported = report_images.get(("vessel_map", vessel_type))
    if exported is not None:
        paths.append(Path(exported))
    if mask_dir:
        paths.extend((
            mask_dir / f"{folder_name}_M0_{vessel_type}.png",
            mask_dir / f"M0_{vessel_type}.png",
        ))
    if hd_png_dir:
        paths.extend((hd_png_dir / f"{folder_name}_M0.png", hd_png_dir / "M0.png"))
    return _load_or_placeholder(paths)


def _try_load_ri_image(report_images, vessel_type):
    path = report_images.get(("ri", vessel_type))
    return None if path is None else _load_or_placeholder([Path(path)])


def _try_load_systole_image(report_images, vessel_type):
    path = report_images.get(("systole", vessel_type))
    return None if path is None else _load_or_placeholder([Path(path)])


def _load_or_placeholder(paths: list[Path]) -> np.ndarray:
    for path in paths:
        if path.exists():
            try:
                with Image.open(path) as image:
                    return np.asarray(image.convert("RGB"))
            except Exception:
                continue
    return np.full((200, 200, 3), 255, dtype=np.uint8)
