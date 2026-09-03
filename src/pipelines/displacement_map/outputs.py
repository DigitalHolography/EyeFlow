"""Temporary numeric caches and magnitude-video output helpers."""

from __future__ import annotations

from pathlib import Path
import warnings

import numpy as np

try:
    import cv2
except ImportError:
    cv2 = None

try:
    from tqdm import tqdm
except ImportError:
    tqdm = None


class OutputCaches:
    def __init__(
        self,
        output_dir: Path,
        frame_count: int,
        height: int,
        width: int,
        save_field: bool,
        valid_mask: np.ndarray,
    ) -> None:
        self.magnitude_path = output_dir / "displacement_magnitude.npy"
        self.magnitude = np.lib.format.open_memmap(
            self.magnitude_path,
            mode="w+",
            dtype=np.float32,
            shape=(frame_count, height, width),
        )

        self.field_path = output_dir / "displacement_field.npy"
        self.field: np.memmap | None = None
        if save_field:
            self.field = np.lib.format.open_memmap(
                self.field_path,
                mode="w+",
                dtype=np.float32,
                shape=(frame_count, height, width, 2),
            )

        self.valid_mask = valid_mask.astype(bool, copy=True)
        self.count = 0

    def append(
        self,
        field: np.ndarray,
        magnitude: np.ndarray,
    ) -> None:
        self.magnitude[self.count] = magnitude

        if self.field is not None:
            self.field[self.count] = field

        self.count += 1

    def flush(self) -> None:
        self.magnitude.flush()
        if self.field is not None:
            self.field.flush()

    def close(self) -> None:
        self.magnitude.flush()
        del self.magnitude

        if self.field is not None:
            self.field.flush()
            del self.field
            self.field = None


def select_display_range(
    cache: OutputCaches,
    mode: str,
    low_percentile: float,
    high_percentile: float,
    fixed_maximum: float,
) -> tuple[float, float]:
    if mode == "fixed":
        return 0.0, float(fixed_maximum)

    values = _median_magnitude_image_values(cache)
    if values.size == 0:
        return 0.0, 0.0

    if mode == "global-minmax":
        return float(np.min(values)), float(np.max(values))

    return (
        float(np.percentile(values, low_percentile)),
        float(np.percentile(values, high_percentile)),
    )


def _median_magnitude_image_values(cache: OutputCaches) -> np.ndarray:
    if cache.count <= 0:
        return np.empty(0, dtype=np.float32)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        median_image = np.nanmedian(
            np.asarray(cache.magnitude[: cache.count], dtype=np.float32),
            axis=0,
        )
    valid = cache.valid_mask & np.isfinite(median_image)
    return np.asarray(median_image[valid], dtype=np.float32)


def create_writer(
    path: Path,
    codec: str,
    fps: float,
    size: tuple[int, int],
) -> cv2.VideoWriter:
    writer = cv2.VideoWriter(
        str(path),
        cv2.VideoWriter_fourcc(*codec),
        fps,
        size,
        True,
    )
    if not writer.isOpened():
        raise RuntimeError(f"Impossible de créer : {path}")
    return writer


def write_magnitude_video(
    cache: OutputCaches,
    output_path: Path,
    codec: str,
    fps: float,
    size: tuple[int, int],
    minimum: float,
    maximum: float,
    gamma: float,
    visualization_sigma: float,
    valid_mask: np.ndarray,
) -> None:
    writer = create_writer(output_path, codec, fps, size)

    try:
        for index in tqdm(
            range(cache.count),
            desc="Écriture magnitude",
            unit="frame",
        ):
            image = np.asarray(
                cache.magnitude[index],
                dtype=np.float32,
            )

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
                normalized = np.clip(
                    (image - minimum) / (maximum - minimum),
                    0.0,
                    1.0,
                )
                normalized = np.power(
                    normalized,
                    gamma,
                    dtype=np.float32,
                )
                gray = np.round(normalized * 255.0).astype(np.uint8)

            gray[~valid_mask] = 0
            writer.write(cv2.cvtColor(gray, cv2.COLOR_GRAY2BGR))
    finally:
        writer.release()
