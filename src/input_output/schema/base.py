"""Small source metadata and reader helpers for EyeFlow HDF5 adapters."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import h5py
import numpy as np

MISSING = object()


@dataclass(frozen=True)
class SourceFileLayout:
    """Companion application file layout metadata used for run discovery."""

    label: str
    companion_suffix: str
    h5_folder_name: str
    h5_filename_template: str
    config_dir_name: str | None = None
    config_filename: str | None = None


class TypedSource:
    """Thin typed facade over one raw HDF5 reader and its sidecar config."""

    def __init__(self, reader, config: dict[str, object] | None = None) -> None:
        self._reader = reader
        self._config = dict(config or {})

    @property
    def filename(self) -> str | None:
        return self._reader.filename

    def _array(
        self,
        path: str,
        *,
        dtype=None,
        default: Any = MISSING,
    ) -> Any:
        if default is MISSING:
            return self._reader.array(path, dtype=dtype)
        return self._reader.array(path, dtype=dtype, default=default)

    def _dataset(self, path: str):
        return self._reader.dataset(path)

    def _value(self, path: str, *, default: Any = MISSING):
        if default is MISSING:
            return self._reader.value(path)
        return self._reader.value(path, default=default)

    def _config_value(self, section: str, key: str, default):
        source = self._config.get(section, {})
        if not isinstance(source, dict):
            return default
        return source.get(key, default)


def scalar_from_value(value):
    """Return the first scalar from a scalar-like HDF5/JSON value."""

    if value is None:
        return None
    array = np.asarray(value).reshape(-1)
    if array.size == 0:
        return None
    scalar = array[0]
    if isinstance(scalar, bytes):
        return scalar.decode("utf-8")
    return scalar.item() if hasattr(scalar, "item") else scalar


def sidecar_dir_for_h5(h5_path: str | Path, folder_name: str) -> Path:
    """Return a sibling sidecar folder next to an exported HDF5 folder."""
    return Path(h5_path).parent.parent / folder_name


def load_h5_sidecar_config(
    h5file: h5py.File | None, *, source: SourceFileLayout
) -> dict[str, object]:
    """Load a source's JSON sidecar while retaining legacy filename fallbacks."""
    if h5file is None or h5file.filename is None:
        return {}
    if not source.config_dir_name or not source.config_filename:
        return {}
    config_path = _sidecar_config_path(Path(h5file.filename), source)
    if config_path is None:
        return {}
    try:
        payload = json.loads(config_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    return _normalize_config_keys(payload)


def _sidecar_config_path(h5_path: Path, source: SourceFileLayout) -> Path | None:
    folders = [source.config_dir_name]
    for fallback in ("json", "config"):
        if fallback not in folders:
            folders.append(fallback)
    for folder_name in folders:
        config_dir = sidecar_dir_for_h5(h5_path, folder_name)
        if not config_dir.is_dir():
            continue
        for name in (
            source.config_filename,
            "parameters_holodoppler.json",
            "parameters_holodoppler",
        ):
            candidate = config_dir / name
            if candidate.is_file():
                return candidate
        json_files = sorted(config_dir.glob("*.json"))
        if json_files:
            return json_files[0]
    return None


def _normalize_config_keys(value):
    if isinstance(value, dict):
        return {
            str(key).replace(" ", ""): _normalize_config_keys(val)
            for key, val in value.items()
        }
    if isinstance(value, list):
        return [_normalize_config_keys(item) for item in value]
    return value
