"""Output paths, destinations, and collision checks for EyeFlow runs."""

import os
import shutil
from collections.abc import Sequence
from dataclasses import dataclass, replace
from enum import Enum
from pathlib import Path

from .runs.layout import HoloRunLayout
from .writers import open_h5, write_png_file
from .writers.artifact_names import prefixed_artifact_path


class OutputType(Enum):
    H5 = "h5"
    PNG = "png"
    MP4 = "mp4"
    AVI = "avi"
    PDF = "pdf"
    EPS = "eps"


def reset_output_dir(path: str | Path) -> None:
    """Replace the EyeFlow-owned output directory for a fresh run."""
    output = Path(path)
    try:
        if output.is_dir():
            shutil.rmtree(output)
        elif output.exists():
            output.unlink()
        output.mkdir(parents=True, exist_ok=False)
    except OSError as exc:
        raise RuntimeError(
            "Could not replace the existing output directory. Close any File Explorer "
            f"window, terminal, or application using this folder, then retry:\n{output}"
        ) from exc


@dataclass(frozen=True)
class OutputManager:
    layout: HoloRunLayout
    artifact_folder: str | None = None

    def for_artifact_namespace(self, namespace: str) -> "OutputManager":
        """Return a manager whose sidecars live below one variant namespace."""

        return OutputManager(self.layout, str(namespace))

    def prepare(self, *, replace: bool = False) -> None:
        output_dir = self.layout.ef_dir
        if replace and output_dir.exists():
            reset_output_dir(output_dir)
        else:
            output_dir.mkdir(parents=True, exist_ok=True)

    def dir_for(self, output_type: OutputType) -> Path:
        directory = self.layout.ef_dir / output_type.value
        if self.artifact_folder is not None:
            directory /= self.artifact_folder
        return directory

    def path_for(
        self,
        output_type: OutputType,
        filename: str | None = None,
    ) -> Path:
        return self.dir_for(output_type) / self._filename_for(output_type, filename)

    def open_h5(self, filename: str | None = None, mode: str = "w"):
        path = self.path_for(OutputType.H5, filename)
        path.parent.mkdir(parents=True, exist_ok=True)
        return open_h5(path, mode)

    def write_png(self, output, filename: str | None = None) -> Path:
        return write_png_file(self.path_for(OutputType.PNG, filename), output)

    def _filename_for(self, output_type: OutputType, filename: str | None) -> str:
        if filename:
            return str(prefixed_artifact_path(filename, self.layout.stem))
        if output_type is OutputType.H5:
            return f"{self.layout.stem}_EF.h5"
        return f"{self.layout.stem}_output.{output_type.value}"


def output_manager_for_layout(
    layout: HoloRunLayout,
    *,
    output_root: Path | None,
    batch_root: Path,
) -> OutputManager:
    if output_root is None:
        return OutputManager(layout)
    relative_path = _relative_to_batch(layout.holo_path, batch_root)
    target_dir = output_root / relative_path.parent
    output_layout = replace(layout, _root_dir=target_dir / layout.stem)
    return OutputManager(output_layout)


def _relative_to_batch(holo_path: Path, batch_root: Path) -> Path:
    try:
        return holo_path.relative_to(batch_root)
    except ValueError:
        anchor = Path(holo_path.anchor)
        drive_token = holo_path.drive.rstrip(":\\/") or "root"
        tail = holo_path.relative_to(anchor) if anchor != holo_path else Path()
        return Path(drive_token) / tail


def reject_duplicate_destinations(requests: Sequence[tuple[HoloRunLayout, OutputManager]]) -> None:
    destinations: dict[str, Path] = {}
    for layout, manager in requests:
        destination = manager.layout.ef_dir.resolve(strict=False)
        key = os.path.normcase(str(destination))
        previous = destinations.get(key)
        if previous is not None:
            raise ValueError(
                "Multiple inputs resolve to the same EyeFlow output directory: "
                f"{previous} and {layout.input_path} -> {destination}"
            )
        destinations[key] = layout.input_path
