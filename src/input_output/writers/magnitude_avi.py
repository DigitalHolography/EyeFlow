"""Display mapping for float displacement-magnitude AVI exports."""

from __future__ import annotations

from pathlib import Path
import warnings

import numpy as np

from .avi import MjpegAviWriter


def select_display_range(
    magnitudes,
    frame_count: int,
    valid_mask: np.ndarray,
    mode: str,
    low_percentile: float,
    high_percentile: float,
    fixed_maximum: float,
) -> tuple[float, float]:
    """Select the grayscale display range from the temporal median magnitude."""
    if mode == "fixed":
        return 0.0, float(fixed_maximum)
    if frame_count <= 0:
        return 0.0, 0.0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        median_image = np.nanmedian(
            np.asarray(magnitudes[:frame_count], dtype=np.float32), axis=0
        )
    values = np.asarray(
        median_image[np.asarray(valid_mask, dtype=bool) & np.isfinite(median_image)],
        dtype=np.float32,
    )
    if values.size == 0:
        return 0.0, 0.0
    if mode == "global-minmax":
        return float(np.min(values)), float(np.max(values))
    return (
        float(np.percentile(values, low_percentile)),
        float(np.percentile(values, high_percentile)),
    )


def write_magnitude_avi(
    path: str | Path,
    magnitudes,
    *,
    frame_count: int,
    valid_mask: np.ndarray,
    fps: float,
    normalization: str,
    low_percentile: float,
    high_percentile: float,
    fixed_maximum: float,
    gamma: float,
    visualization_sigma: float,
) -> tuple[float, float]:
    """Render float displacement magnitudes as a grayscale MJPEG AVI."""
    from tqdm import tqdm

    target = Path(path)
    if target.suffix.lower() != ".avi":
        raise ValueError("Magnitude AVI output path must end in '.avi'.")
    height, width = magnitudes.shape[-2:]
    mask = np.asarray(valid_mask, dtype=bool)
    if mask.shape != (height, width):
        raise ValueError("Magnitude AVI mask does not match the frame shape.")
    minimum, maximum = select_display_range(
        magnitudes,
        frame_count,
        mask,
        normalization,
        low_percentile,
        high_percentile,
        fixed_maximum,
    )
    if visualization_sigma > 0:
        import cv2

    with MjpegAviWriter(target, width=width, height=height, fps=fps) as writer:
        for index in tqdm(range(frame_count), desc="Écriture magnitude", unit="frame"):
            image = np.asarray(magnitudes[index], dtype=np.float32)
            if visualization_sigma > 0:
                image = cv2.GaussianBlur(
                    image,
                    (0, 0),
                    visualization_sigma,
                    borderType=cv2.BORDER_REFLECT101,
                )
            if maximum <= minimum + 1e-12:
                gray = np.zeros(image.shape, np.uint8)
            else:
                normalized = np.clip((image - minimum) / (maximum - minimum), 0.0, 1.0)
                normalized = np.power(normalized, gamma, dtype=np.float32)
                gray = np.round(normalized * 255.0).astype(np.uint8)
            gray[~mask] = 0
            writer.write_frame(gray)
    return minimum, maximum
