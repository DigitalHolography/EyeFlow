"""Discover HOLO markers and matching HD/DV run folders without opening HDF5."""

from __future__ import annotations

import os
from pathlib import Path

from .holo_run_layout import INPUT_LAYOUTS, HoloRunLayout


def is_run_directory(directory: Path) -> bool:
    """Require matching HD and DV directories; existing EF/AE output is optional."""
    directory = directory.expanduser().resolve()
    return directory.is_dir() and all(
        (directory / f"{directory.name}_{layout.companion_suffix}").is_dir()
        for layout in INPUT_LAYOUTS
    )


def run_input_key(path: Path) -> str:
    """Deduplicate a marker and its companion-data directory as one run."""
    identity = path
    if path.suffix.lower() == ".holo" and not path.is_dir():
        identity = HoloRunLayout.from_holo(path).root_dir
    return os.path.normcase(str(identity.resolve()))


def discover_run_inputs(directory: Path) -> list[Path]:
    """Search the selected directory itself and its descendants for input runs.

    Keep existing markers when available. Otherwise return the real run folder;
    the canonical sibling .holo path is used only as an internal naming anchor.
    Validation of HDF5 files belongs to normal run-spec resolution.
    """
    directory = directory.expanduser().resolve()
    if not directory.is_dir():
        raise NotADirectoryError(f"Input folder does not exist: {directory}")
    markers: list[Path] = []
    run_directories: list[Path] = []

    def fail_scan(error: OSError) -> None:
        raise error

    for root_name, folders, filenames in os.walk(directory, onerror=fail_scan):
        root = Path(root_name)
        required = {os.path.normcase(f"{root.name}_{layout.companion_suffix}")
                    for layout in INPUT_LAYOUTS}
        if required.issubset(os.path.normcase(name) for name in folders):
            run_directories.append(root)
        markers.extend(root / name for name in filenames
                       if Path(name).suffix.lower() == ".holo" and (root / name).is_file())

    found: dict[str, Path] = {}
    for marker in sorted(markers, key=lambda path: str(path).casefold()):
        found.setdefault(run_input_key(marker), marker)
    for root in run_directories:
        marker = root.parent / f"{root.name}.holo"
        found.setdefault(run_input_key(root), marker if marker.is_file() else root)
    return sorted(found.values(), key=lambda path: str(path).casefold())
