"""Persisted and temporary displacement-map array storage."""

from __future__ import annotations

from pathlib import Path
from collections.abc import Mapping

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

def load_displacement_maps(
    field_paths_by_vessel: dict[str, Path], method: str,
) -> dict[str, dict[str, object]]:
    """Open each persisted field once, even when vessels share a file."""
    loaded_by_path: dict[str, object] = {}
    displacement_maps: dict[str, dict[str, object]] = {}
    for vessel, field_path in field_paths_by_vessel.items():
        normalized_path = str(field_path.resolve())
        displacement_map = loaded_by_path.get(normalized_path)
        if displacement_map is None:
            displacement_map = np.load(field_path, mmap_mode="r")
            loaded_by_path[normalized_path] = displacement_map
        displacement_maps[vessel] = {method: displacement_map}
    return displacement_maps


def release_displacement_maps(
    displacement_maps: Mapping[str, Mapping[str, object]],
) -> None:
    """Close shared memory maps exactly once."""
    closed: set[int] = set()
    for maps_for_vessel in displacement_maps.values():
        for displacement_map in maps_for_vessel.values():
            identity = id(displacement_map)
            if identity in closed:
                continue
            closed.add(identity)
            mmap = getattr(displacement_map, "_mmap", None)
            if mmap is not None:
                mmap.close()
