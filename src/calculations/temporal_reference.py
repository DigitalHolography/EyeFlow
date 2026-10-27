"""Compute a mean reference image from a frame sequence."""

from __future__ import annotations

from typing import Iterator, Protocol

import numpy as np

try:
    import cv2
except ImportError:
    cv2 = None

try:
    from tqdm import tqdm
except ImportError:
    tqdm = None


class FrameSource(Protocol):
    frame_count: int

    def iter_frames(self, max_frames: int | None = None) -> Iterator[np.ndarray]:
        ...


def compute_mean_reference(
    sequence: FrameSource,
    max_frames: int | None,
) -> tuple[np.ndarray, int]:
    accumulator: np.ndarray | None = None
    count = 0
    total = sequence.frame_count if max_frames is None else min(sequence.frame_count, max_frames)
    for frame in tqdm(
        sequence.iter_frames(max_frames), total=total, desc="Référence moyenne", unit="frame"
    ):
        gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY).astype(np.float64) / 255.0
        if accumulator is None:
            accumulator = np.zeros_like(gray, dtype=np.float64)
        accumulator += gray
        count += 1
    if accumulator is None or count == 0:
        raise RuntimeError("Aucune frame lisible.")
    return (accumulator / count).astype(np.float32), count
