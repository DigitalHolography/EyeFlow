"""Resolve and expand selected HOLO markers, run folders, lists, and ZIP inputs."""

import os
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path

from .layout import HoloRunLayout
from .discovery import discover_run_inputs, is_run_directory
from ..archives import extracted_zip_tree

HOLO_SUFFIX = ".holo"
INPUT_LIST_SUFFIX = ".txt"


@dataclass(frozen=True)
class HoloInputStatus:
    hd: bool
    dv: bool


@dataclass(frozen=True)
class HoloInputList:
    path_stem_pairs: tuple[tuple[Path, str], ...]


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
                raise ValueError(f"No {HOLO_SUFFIX} files or HD/DV run folders found in {source}")
            batch_root = extracted_root.parent if is_run_directory(extracted_root) else extracted_root
            yield ExpandedRunInputs(tuple(paths), batch_root)
        return

    paths = _find_holo_inputs(source)
    if not paths:
        raise ValueError(f"No {HOLO_SUFFIX} files or HD/DV run folders found under {source}")
    batch_root = (source.parent if is_run_directory(source) else source) if source.is_dir() else None
    yield ExpandedRunInputs(tuple(paths), batch_root)


def _find_holo_inputs(path: Path) -> list[Path]:
    if path.is_file():
        if path.suffix.lower() in {HOLO_SUFFIX, INPUT_LIST_SUFFIX}:
            return [path]
        raise ValueError(
            f"File is not a {HOLO_SUFFIX} or {INPUT_LIST_SUFFIX} file: {path}"
        )
    if path.is_dir():
        return discover_run_inputs(path)
    raise FileNotFoundError(f"Input path does not exist: {path}")


def batch_root(holo_paths: Sequence[Path]) -> Path:
    if not holo_paths:
        return Path.cwd()
    if len(holo_paths) == 1:
        return holo_paths[0].parent
    try:
        return Path(os.path.commonpath([str(path.parent) for path in holo_paths]))
    except ValueError:
        return Path.cwd()


def resolve_holo_run_layout(
    holo_path: Path,
) -> HoloRunLayout:
    holo_path = _absolute(holo_path)
    _validate_holo_file(holo_path)
    run_layout = HoloRunLayout.from_holo(holo_path)
    run_layout.require_inputs()
    return run_layout


def resolve_stem_run_layout(stem: str, root_dir: Path) -> HoloRunLayout:
    root_dir = _absolute(root_dir)
    run_layout = HoloRunLayout(
        _holo_path=root_dir / f"{stem}{HOLO_SUFFIX}",
        _stem=stem,
        _root_dir=root_dir / stem,
    )
    run_layout.require_inputs()
    return run_layout


def resolve_selected_run_layouts(
    input_paths: Sequence[Path],
) -> list[HoloRunLayout]:
    normalized = [_absolute(path) for path in input_paths]
    if not normalized:
        raise ValueError(
            f"Select {HOLO_SUFFIX} files, matching HD/DV run folders, "
            f"or one {INPUT_LIST_SUFFIX} list."
        )
    if len(normalized) == 1 and normalized[0].suffix.lower() == INPUT_LIST_SUFFIX and (
        not normalized[0].is_dir()
    ):
        input_list_path = normalized[0]
        input_list = read_holo_input_list(input_list_path)
        return [
            resolve_stem_run_layout(stem, root_dir)
            for root_dir, stem in input_list.path_stem_pairs
        ]
    if any(path.suffix.lower() == INPUT_LIST_SUFFIX and not path.is_dir() for path in normalized):
        raise ValueError(
            f"Select either one {INPUT_LIST_SUFFIX} list or "
            f"{HOLO_SUFFIX} files and HD/DV run folders."
        )

    resolved: list[HoloRunLayout] = []
    errors: list[str] = []
    for holo_path in normalized:
        try:
            if holo_path.is_dir():
                if not is_run_directory(holo_path):
                    raise ValueError(
                        f"Run folder must contain {holo_path.name}_HD and "
                        f"{holo_path.name}_DV directories:\n{holo_path}"
                    )
                layout = HoloRunLayout.from_run_directory(holo_path)
                layout.require_inputs()
                resolved.append(layout)
            else:
                resolved.append(resolve_holo_run_layout(holo_path))
        except (FileNotFoundError, ValueError) as exc:
            errors.append(f"{holo_path}:\n{exc}")

    if errors:
        raise FileNotFoundError(
            "Missing required input data for one or more selected "
            "runs:\n\n"
            + "\n\n".join(errors)
        )
    return resolved


def read_holo_input_list(input_list_path: Path) -> HoloInputList:
    input_list_path = _absolute(input_list_path)
    entries: list[str] = []
    for line in input_list_path.read_text(encoding="utf-8-sig").splitlines():
        entry = line.strip()
        if entry[:1] in {"\"", "'"} and entry[-1:] == entry[:1]:
            entry = entry[1:-1]
        if entry:
            entries.append(entry)
    if not entries:
        raise ValueError(f"Input list is empty:\n{input_list_path}")
    return HoloInputList(
        path_stem_pairs=tuple(
            _parse_input_list_entries(entries, input_list_path.parent)
        )
    )


def holo_input_status(
    holo_path: Path,
) -> HoloInputStatus:
    holo_path = _absolute(holo_path)
    if holo_path.is_dir():
        layout = HoloRunLayout.from_run_directory(holo_path)
        return HoloInputStatus(hd=layout.has_hd_h5, dv=layout.has_dv_h5)
    try:
        _validate_holo_file(holo_path)
    except (FileNotFoundError, ValueError):
        return HoloInputStatus(hd=False, dv=False)

    run_layout = HoloRunLayout.from_holo(holo_path)
    return HoloInputStatus(
        hd=run_layout.has_hd_h5,
        dv=run_layout.has_dv_h5,
    )


def stem_input_status(stem: str, root_dir: Path) -> HoloInputStatus:
    root_dir = _absolute(root_dir)
    run_layout = HoloRunLayout(
        _holo_path=root_dir / f"{stem}{HOLO_SUFFIX}",
        _stem=stem,
        _root_dir=root_dir / stem,
    )
    return HoloInputStatus(
        hd=run_layout.has_hd_h5,
        dv=run_layout.has_dv_h5,
    )


def _parse_input_list_entries(
    entries: Sequence[str],
    default_root_dir: Path,
) -> list[tuple[Path, str]]:
    parsed: list[tuple[Path, str]] = []
    for entry in entries:
        path = Path(entry).expanduser()
        if path.suffix.lower() == HOLO_SUFFIX:
            holo_path = path if path.is_absolute() else default_root_dir / path
            holo_path = _absolute(holo_path)
            parsed.append((holo_path.parent, holo_path.stem))
        else:
            parsed.append((default_root_dir, entry))
    return parsed


def _absolute(path: str | Path) -> Path:
    resolved = Path(path).expanduser()
    return resolved if resolved.is_absolute() else Path.cwd() / resolved


def _validate_holo_file(holo_path: Path) -> None:
    if holo_path.suffix.lower() != HOLO_SUFFIX:
        raise ValueError(f"HOLO input must be a {HOLO_SUFFIX} file:\n{holo_path}")
    if not holo_path.exists():
        raise FileNotFoundError(f"HOLO input does not exist:\n{holo_path}")
    if not holo_path.is_file():
        raise ValueError(f"HOLO input must be a file:\n{holo_path}")
