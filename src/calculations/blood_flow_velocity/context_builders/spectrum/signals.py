"""Build display-ready signal and histogram inputs for spectrum figures."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


FREQUENCY_DISPLAY_SCALE = np.float32(1.0 / 1000.0)
VELOCITY_DISPLAY_SCALE = np.float32(1.0 / 1000.0)


@dataclass(frozen=True)
class HistogramData:
    counts: np.ndarray
    vmin: float
    vmax: float
    count_max: float


def display_frequency(values) -> np.ndarray:
    return np.asarray(values, dtype=np.float32) * FREQUENCY_DISPLAY_SCALE


def display_velocity(values) -> np.ndarray:
    return np.asarray(values, dtype=np.float32) * VELOCITY_DISPLAY_SCALE


def section_mask(analysis: dict[str, object], shape: tuple[int, int]) -> np.ndarray:
    section = analysis.get("velocity_section_mask")
    if section is not None:
        return np.asarray(section, dtype=bool)
    return np.ones(shape, dtype=bool)


def array_or_none(values) -> np.ndarray | None:
    if values is None:
        return None
    return np.asarray(values, dtype=np.float32)


def masked_video_signal(video: np.ndarray, mask: np.ndarray) -> np.ndarray:
    data = np.asarray(video, dtype=np.float32)
    count = np.count_nonzero(mask)
    if count == 0:
        return np.full((data.shape[0],), np.nan, dtype=np.float32)
    return np.nanmean(np.where(mask[None, :, :], data, np.nan), axis=(1, 2)).astype(np.float32)


def histogram_matrix(video: np.ndarray, mask: np.ndarray, bins: int = 256) -> HistogramData:
    data = np.asarray(video, dtype=np.float32)
    selected = data[:, mask]
    selected = selected[np.isfinite(selected)]
    if selected.size == 0:
        return HistogramData(
            np.zeros((bins, data.shape[0]), dtype=np.float32),
            0.0,
            1.0,
            1.0,
        )
    vmin = float(np.nanmin(selected))
    vmax = float(np.nanmax(selected))
    if vmax <= vmin:
        vmax = vmin + 1.0
    counts = np.zeros((bins, data.shape[0]), dtype=np.float32)
    edges = np.linspace(vmin, vmax, bins + 1, dtype=np.float32)
    for frame_idx, frame in enumerate(data):
        counts[:, frame_idx] = np.histogram(frame[mask], bins=edges)[0]
    count_max = float(np.nanmax(counts))
    return HistogramData(counts, vmin, vmax, count_max if count_max > 0 else 1.0)
