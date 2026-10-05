"""Temporal reference computation for displacement maps."""

from __future__ import annotations

import numpy as np

from input_output.frame_sequences import FrameSequence

try:
    import cv2
except ImportError:
    cv2 = None

try:
    from tqdm import tqdm
except ImportError:
    tqdm = None


def compute_mean_reference(
    sequence: FrameSequence,
    max_frames: int | None,
) -> tuple[np.ndarray, int]:
    accumulator: np.ndarray | None = None
    count = 0
    total = sequence.frame_count if max_frames is None else min(sequence.frame_count, max_frames)
    for frame in tqdm(
        sequence.iter_frames(max_frames), total=total, desc="RÃ©fÃ©rence moyenne", unit="frame"
    ):
        gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY).astype(np.float64) / 255.0
        if accumulator is None:
            accumulator = np.zeros_like(gray, dtype=np.float64)
        accumulator += gray
        count += 1
    if accumulator is None or count == 0:
        raise RuntimeError("Aucune frame lisible.")
    return (accumulator / count).astype(np.float32), count
