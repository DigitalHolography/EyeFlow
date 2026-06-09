"""Resolve HOLO selections and expose HD/DV/work HDF5 inputs to pipelines."""

import json
from collections.abc import Iterator, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import h5py

from .holo_run_layout import HoloRunLayout
from .schema import DOPPLER_VIEW_LAYOUT, HOLODOPPLER_LAYOUT, SourceFileLayout
from .writers.h5 import normalize_h5_path

HOLO_SUFFIX = ".holo"
INPUT_LIST_SUFFIX = ".txt"


@dataclass(frozen=True)
class HoloInputStatus:
    hd: bool
    dv: bool
    ef: bool


def resolve_holo_run_layout(
    holo_path: Path,
    *,
    require_holo_file: bool = True,
) -> HoloRunLayout:
    holo_path = _absolute(holo_path)
    _validate_holo_file(holo_path, require_file=require_holo_file)
    run_layout = HoloRunLayout.from_holo(holo_path)
    run_layout.require_inputs()
    return run_layout


def resolve_stem_run_layout(stem: str, root_dir: Path) -> HoloRunLayout:
    root_dir = _absolute(root_dir)
    run_layout = HoloRunLayout(
        _holo_path=root_dir / stem,
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
        return [
            resolve_stem_run_layout(stem, input_list_path.parent)
            for stem in read_stems_from_input_list(input_list_path)
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


def read_stems_from_input_list(input_list_path: Path) -> list[str]:
    stems = [
        line.strip()
        for line in _absolute(input_list_path).read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if not stems:
        raise ValueError(f"Input list is empty:\n{input_list_path}")
    return stems


def default_output_dir_for_input(input_path: Path) -> Path:
    if input_path.suffix.lower() == HOLO_SUFFIX:
        return HoloRunLayout.from_holo(input_path).ef_dir
    return input_path.parent if input_path.is_file() else input_path


def holo_input_status(
    holo_path: Path,
    *,
    require_holo_file: bool,
) -> HoloInputStatus:
    holo_path = _absolute(holo_path)
    try:
        _validate_holo_file(holo_path, require_file=require_holo_file)
    except (FileNotFoundError, ValueError):
        return HoloInputStatus(hd=False, dv=False, ef=False)

    run_layout = HoloRunLayout.from_holo(holo_path)
    return HoloInputStatus(
        hd=run_layout.has_hd_h5,
        dv=run_layout.has_dv_h5,
        ef=run_layout.has_ef_h5,
    )


def stem_input_status(stem: str, root_dir: Path) -> HoloInputStatus:
    root_dir = _absolute(root_dir)
    run_layout = HoloRunLayout(
        _holo_path=root_dir / stem,
        _stem=stem,
        _root_dir=root_dir / stem,
    )
    return HoloInputStatus(
        hd=run_layout.has_hd_h5,
        dv=run_layout.has_dv_h5,
        ef=run_layout.has_ef_h5,
    )


def _lookup_key(path: str) -> str:
    return normalize_h5_path(path)


def _absolute(path: str | Path) -> Path:
    resolved = Path(path).expanduser()
    return resolved if resolved.is_absolute() else Path.cwd() / resolved


def _validate_holo_file(holo_path: Path, *, require_file: bool) -> None:
    if holo_path.suffix.lower() != HOLO_SUFFIX:
        raise ValueError(f"HOLO input must be a {HOLO_SUFFIX} file:\n{holo_path}")
    if not require_file:
        return
    if not holo_path.exists():
        raise FileNotFoundError(f"HOLO input does not exist:\n{holo_path}")
    if not holo_path.is_file():
        raise ValueError(f"HOLO input must be a file:\n{holo_path}")


class MergedAttrs(Mapping[str, object]):
    def __init__(self, *sources: h5py.File | Mapping[str, object] | None) -> None:
        self._sources = [
            _attr_source(source) for source in sources if source is not None
        ]

    def __getitem__(self, key: str) -> object:
        sentinel = object()
        value = self.get(key, sentinel)
        if value is sentinel:
            raise KeyError(key)
        return value

    def __iter__(self) -> Iterator[str]:
        seen: set[str] = set()
        for source in self._sources:
            for key in source.keys():
                if key not in seen:
                    seen.add(key)
                    yield str(key)

    def __len__(self) -> int:
        return sum(1 for _ in self.__iter__())

    def get(self, key: str, default=None):
        for source in self._sources:
            if key in source:
                return source[key]
        return default


class EyeFlowView:
    def __init__(self, work_h5: h5py.File) -> None:
        self.work_h5 = work_h5

    def get(self, key: str, default=None):
        normalized_key = _lookup_key(key)
        if not normalized_key:
            return default
        return self.work_h5.get(normalized_key, default)

    def __getitem__(self, key: str):
        found = self.get(key)
        if found is None:
            raise KeyError(key)
        return found

    def __contains__(self, key: object) -> bool:
        return isinstance(key, str) and self.get(key) is not None


def _attr_source(source: h5py.File | Mapping[str, object]) -> Mapping[str, object]:
    return source.attrs if isinstance(source, h5py.File) else source


def _load_sidecar_config(
    h5file: h5py.File | None,
    *,
    source_schema: SourceFileLayout,
) -> dict[str, object]:
    if h5file is None or h5file.filename is None:
        return {}
    if not source_schema.config_dir_name or not source_schema.config_filename:
        return {}
    config_path = _sidecar_config_path(
        Path(h5file.filename),
        folder_name=source_schema.config_dir_name,
        preferred_name=source_schema.config_filename,
    )
    if config_path is None:
        return {}
    try:
        payload = json.loads(config_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    return _normalize_config_keys(payload)


def load_h5_sidecar_config(
    h5file: h5py.File | None,
    *,
    source: str | SourceFileLayout,
) -> dict[str, object]:
    if isinstance(source, SourceFileLayout):
        source_schema = source
    elif source == "hd":
        source_schema = HOLODOPPLER_LAYOUT
    elif source == "dv":
        source_schema = DOPPLER_VIEW_LAYOUT
    else:
        raise ValueError(f"Unknown sidecar config source: {source}")
    return _load_sidecar_config(h5file, source_schema=source_schema)


def _sidecar_config_path(
    h5_path: Path,
    *,
    folder_name: str,
    preferred_name: str,
) -> Path | None:
    config_dir = h5_path.parent.parent / folder_name
    if not config_dir.is_dir():
        return None
    preferred = config_dir / preferred_name
    if preferred.is_file():
        return preferred
    json_files = sorted(config_dir.glob("*.json"))
    return json_files[0] if json_files else None


def _normalize_config_keys(value):
    if isinstance(value, dict):
        return {
            str(key).replace(" ", ""): _normalize_config_keys(val)
            for key, val in value.items()
        }
    if isinstance(value, list):
        return [_normalize_config_keys(item) for item in value]
    return value
