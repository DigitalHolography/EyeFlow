"""Resolve selected HOLO markers, run folders, and input lists."""

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import h5py

from .holo_run_layout import HoloRunLayout
from .run_discovery import is_run_directory
from .schema import DOPPLER_VIEW_LAYOUT, HOLODOPPLER_LAYOUT, SourceFileLayout
from .schema.base import load_h5_sidecar_config as _load_h5_sidecar_config
from .schema.base import sidecar_dir_for_h5

HOLO_SUFFIX = ".holo"
INPUT_LIST_SUFFIX = ".txt"


@dataclass(frozen=True)
class HoloInputStatus:
    hd: bool
    dv: bool


@dataclass(frozen=True)
class HoloInputList:
    path_stem_pairs: tuple[tuple[Path, str], ...]


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


def load_h5_sidecar_config(
    h5file: h5py.File | None, *, source: str | SourceFileLayout
) -> dict[str, object]:
    """Resolve a named source and load its optional sidecar JSON."""
    if isinstance(source, SourceFileLayout):
        layout = source
    elif source == "hd":
        layout = HOLODOPPLER_LAYOUT
    elif source == "dv":
        layout = DOPPLER_VIEW_LAYOUT
    else:
        raise ValueError(f"Unknown sidecar config source: {source}")
    return _load_h5_sidecar_config(h5file, source=layout)
