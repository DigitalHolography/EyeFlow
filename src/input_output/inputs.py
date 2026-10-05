"""Resolve selected HOLO runs and their HD/DV input paths."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import h5py

from .schema import DOPPLER_VIEW_LAYOUT, HOLODOPPLER_LAYOUT, SourceFileLayout

HOLO_SUFFIX = ".holo"
INPUT_LIST_SUFFIX = ".txt"
HDF5_SUFFIXES = (".h5", ".hdf5")
INPUT_LAYOUTS = (HOLODOPPLER_LAYOUT, DOPPLER_VIEW_LAYOUT)


@dataclass(frozen=True)
class HoloRunLayout:
    """Path layout and input discovery for one selected HOLO run."""

    _holo_path: Path
    _stem: str
    _root_dir: Path

    @classmethod
    def from_holo(
        cls,
        holo_path: str | Path,
        *,
        output_root: str | Path | None = None,
    ) -> HoloRunLayout:
        path = _absolute(Path(holo_path).expanduser())
        root = _absolute(Path(output_root).expanduser()) if output_root else path.parent
        return cls(_holo_path=path, _stem=path.stem, _root_dir=root / path.stem)

    @property
    def holo_path(self) -> Path:
        return self._holo_path

    @property
    def stem(self) -> str:
        return self._stem

    @property
    def root_dir(self) -> Path:
        return self._root_dir

    @property
    def ef_dir(self) -> Path:
        return self._run_dir("EF")

    @property
    def hd_h5(self) -> Path:
        return self._input_h5(HOLODOPPLER_LAYOUT)

    @property
    def dv_h5(self) -> Path:
        return self._input_h5(DOPPLER_VIEW_LAYOUT)

    @property
    def has_hd_h5(self) -> bool:
        return self._has_h5(HOLODOPPLER_LAYOUT)

    @property
    def has_dv_h5(self) -> bool:
        return self._has_h5(DOPPLER_VIEW_LAYOUT)

    def require_inputs(self) -> None:
        if not self.root_dir.is_dir():
            raise FileNotFoundError(f"Could not find data folder:\n{self.root_dir}")
        errors: list[str] = []
        for schema in INPUT_LAYOUTS:
            try:
                self._input_h5(schema)
            except FileNotFoundError as exc:
                errors.append(str(exc))
        if errors:
            raise FileNotFoundError(
                "Missing required input data for the selected .holo file:\n\n"
                + "\n\n".join(errors)
            )

    def _run_dir(self, suffix: str) -> Path:
        return self._root_dir / f"{self._stem}_{suffix}"

    def _input_dir(self, schema: SourceFileLayout) -> Path:
        return self._run_dir(schema.companion_suffix)

    def _h5_dir(self, schema: SourceFileLayout) -> Path:
        return self._input_dir(schema) / schema.h5_folder_name

    def _preferred_h5(self, schema: SourceFileLayout) -> Path:
        folder_name = f"{self._stem}_{schema.companion_suffix}"
        filename = schema.h5_filename_template.format(
            stem=self._stem,
            folder=folder_name,
            companion=schema.companion_suffix,
        )
        return self._h5_dir(schema) / filename

    def _input_h5(self, schema: SourceFileLayout) -> Path:
        files = self._h5_files(schema)
        if not files:
            raise FileNotFoundError(
                f"{schema.label} HDF5 file missing in expected folder:\n"
                f"{self._h5_dir(schema)}"
            )
        preferred = self._preferred_h5(schema)
        if preferred in files:
            return preferred
        if len(files) == 1:
            return files[0]
        candidates = "\n".join(str(path) for path in files)
        raise FileNotFoundError(
            f"Multiple {schema.label} HDF5 files found in:\n{self._h5_dir(schema)}\n\n"
            f"Expected one file, preferably named:\n{preferred.name}\n\n"
            f"Candidates:\n{candidates}"
        )

    def _has_h5(self, schema: SourceFileLayout) -> bool:
        try:
            self._input_h5(schema)
        except FileNotFoundError:
            return False
        return True

    def _h5_files(self, schema: SourceFileLayout) -> list[Path]:
        folder = self._h5_dir(schema)
        if not folder.is_dir():
            return []
        return sorted(path for path in folder.iterdir() if _is_hdf5_file(path))


def _is_hdf5_file(path: Path) -> bool:
    return path.is_file() and path.suffix.lower() in HDF5_SUFFIXES and h5py.is_hdf5(path)


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
            f"Select one or more {HOLO_SUFFIX} files "
            f"or one {INPUT_LIST_SUFFIX} list."
        )
    if len(normalized) == 1 and normalized[0].suffix.lower() == INPUT_LIST_SUFFIX:
        input_list_path = normalized[0]
        input_list = read_holo_input_list(input_list_path)
        return [
            resolve_stem_run_layout(stem, root_dir)
            for root_dir, stem in input_list.path_stem_pairs
        ]
    if any(path.suffix.lower() == INPUT_LIST_SUFFIX for path in normalized):
        raise ValueError(
            f"Select either one {INPUT_LIST_SUFFIX} list or one or more "
            f"{HOLO_SUFFIX} files."
        )

    resolved: list[HoloRunLayout] = []
    errors: list[str] = []
    for holo_path in normalized:
        try:
            resolved.append(resolve_holo_run_layout(holo_path))
        except (FileNotFoundError, ValueError) as exc:
            errors.append(f"{holo_path}:\n{exc}")

    if errors:
        raise FileNotFoundError(
            "Missing required input data for one or more selected "
            f"{HOLO_SUFFIX} files:\n\n"
            + "\n\n".join(errors)
        )
    return resolved


def read_holo_input_list(input_list_path: Path) -> HoloInputList:
    input_list_path = _absolute(input_list_path)
    entries = [
        line.strip()
        for line in input_list_path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
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
