"""Temporary numeric caches for displacement-map calculation."""

from __future__ import annotations

from pathlib import Path

import numpy as np


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
