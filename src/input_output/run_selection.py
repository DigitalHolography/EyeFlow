"""Input expansion and output destinations for batch runs."""

from __future__ import annotations

import os
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path

from .archives import extracted_zip_tree
from .inputs import INPUT_LIST_SUFFIX, HoloRunLayout
from .output_manager import OutputManager

HOLO_SUFFIX = ".holo"


@dataclass(frozen=True)
class ExpandedRunInputs:
    """Input paths produced from one CLI source selection."""

    paths: tuple[Path, ...]
    batch_root: Path | None = None


@contextmanager
def expand_run_inputs(data_path: Path) -> Iterator[ExpandedRunInputs]:
    """Expand a CLI file, recursive folder, or ZIP source safely."""

    source = data_path.expanduser().resolve()
    if source.is_file() and source.suffix.lower() == ".zip":
        with extracted_zip_tree(source) as extracted_root:
            paths = _find_holo_inputs(extracted_root)
            if not paths:
                raise ValueError(f"No {HOLO_SUFFIX} files found in {source}")
            yield ExpandedRunInputs(tuple(paths), extracted_root)
        return

    paths = _find_holo_inputs(source)
    if not paths:
        raise ValueError(f"No {HOLO_SUFFIX} files found under {source}")
    batch_root = source if source.is_dir() else None
    yield ExpandedRunInputs(tuple(paths), batch_root)


def _find_holo_inputs(path: Path) -> list[Path]:
    if path.is_file():
        if path.suffix.lower() in {HOLO_SUFFIX, INPUT_LIST_SUFFIX}:
            return [path]
        raise ValueError(
            f"File is not a {HOLO_SUFFIX} or {INPUT_LIST_SUFFIX} file: {path}"
        )
    if path.is_dir():
        return sorted(
            candidate
            for candidate in path.rglob("*")
            if candidate.is_file() and candidate.suffix.lower() == HOLO_SUFFIX
        )
    raise FileNotFoundError(f"Input path does not exist: {path}")


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
    output_layout = HoloRunLayout.from_holo(
        layout.holo_path,
        output_root=target_dir,
    )
    return OutputManager(output_layout)


def batch_root(holo_paths: Sequence[Path]) -> Path:
    if not holo_paths:
        return Path.cwd()
    if len(holo_paths) == 1:
        return holo_paths[0].parent
    try:
        return Path(os.path.commonpath([str(path.parent) for path in holo_paths]))
    except ValueError:
        return Path.cwd()


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
                f"{previous} and {layout.holo_path} -> {destination}"
            )
        destinations[key] = layout.holo_path
