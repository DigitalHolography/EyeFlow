"""Build display-ready signal and histogram inputs for pulse figures."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


FREQUENCY_DISPLAY_SCALE = np.float32(1.0 / 1000.0)
VELOCITY_DISPLAY_SCALE = np.float32(1.0)


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


def histogram_matrix(video: np.ndarray, mask: np.ndarray, bins: int = 256) -> HistogramData:
    vmin = np.inf
    vmax = -np.inf
    has_finite = False
    for _, chunk in video_chunks(video):
        selected = chunk[:, mask]
        finite = selected[np.isfinite(selected)]
        if finite.size:
            has_finite = True
            vmin = min(vmin, float(np.min(finite)))
            vmax = max(vmax, float(np.max(finite)))
    if not has_finite:
        return HistogramData(
            np.zeros((bins, video.shape[0]), dtype=np.float32),
            0.0,
            1.0,
            1.0,
        )
    if vmax <= vmin:
        vmax = vmin + 1.0
    counts = np.zeros((bins, video.shape[0]), dtype=np.float32)
    edges = np.linspace(vmin, vmax, bins + 1, dtype=np.float32)
    for start, chunk in video_chunks(video):
        for offset, frame in enumerate(chunk):
            counts[:, start + offset] = np.histogram(frame[mask], bins=edges)[0]
    count_max = float(np.nanmax(counts))
    return HistogramData(counts, vmin, vmax, count_max if count_max > 0 else 1.0)


def video_chunks(video, frame_chunk_size: int = 8):
    for start in range(0, int(video.shape[0]), frame_chunk_size):
        stop = min(start + frame_chunk_size, int(video.shape[0]))
        yield start, np.asarray(video[start:stop], dtype=np.float32)
