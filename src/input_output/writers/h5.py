"""Write EyeFlow runtime values into HDF5 files."""

import json
import re
from pathlib import Path

import h5py
import numpy as np

from app_settings import app_version
from ..schema.holodoppler import HD_OUTPUT_PASSTHROUGH_PATHS


def normalize_h5_path(path: object) -> str:
    return str(path).replace("\\", "/").strip("/")


def open_h5(path: Path | str, mode: str = "r") -> h5py.File:
    h5_path = Path(path)
    try:
        return h5py.File(h5_path, mode)
    except OSError as exc:
        raise OSError(f"Could not open HDF5 file '{h5_path}': {exc}") from exc


def resolve_dataset_target(root_group: h5py.Group, key: str) -> tuple[h5py.Group, str]:
    normalized_key = normalize_h5_path(key)
    parts = [part for part in normalized_key.split("/") if part]
    if not parts:
        raise ValueError("Dataset key cannot be empty.")

    parent = root_group
    for part in parts[:-1]:
        existing = parent.get(part)
        if existing is None:
            parent = parent.create_group(part)
            continue
        if isinstance(existing, h5py.Group):
            parent = existing
            continue
        raise ValueError(
            f"Cannot create subgroup '{part}' for key '{key}': "
            "a dataset already exists at that path."
        )

    return parent, parts[-1]


def set_attr_safe(h5obj: h5py.File | h5py.Group | h5py.Dataset, key: str, value) -> None:
    if isinstance(value, str):
        h5obj.attrs.create(key, value, dtype=h5py.string_dtype(encoding="utf-8"))
        return
    data = value
    if isinstance(value, (list, tuple)):
        if all(isinstance(item, str) for item in value):
            data = np.asarray(value, dtype=h5py.string_dtype(encoding="utf-8"))
        else:
            data = np.asarray(value)
    data = _downcast_numeric_payload(data)
    try:
        h5obj.attrs[key] = data
    except (TypeError, ValueError):
        h5obj.attrs[key] = str(value)


def write_value_dataset(group: h5py.Group, key: str, value) -> None:
    from pipeline_engine import DatasetValue

    ds_attrs = None
    data = value

    if hasattr(value, "data") and hasattr(value, "attrs"):
        data = value.data
        ds_attrs = value.attrs
    elif isinstance(value, DatasetValue):
        data = value.data
        ds_attrs = value.attrs
    elif isinstance(value, tuple) and len(value) == 2 and isinstance(value[1], dict):
        data, ds_attrs = value

    target_group, dataset_key = resolve_dataset_target(group, str(key))
    if dataset_key in target_group:
        del target_group[dataset_key]

    payload, ds_attrs = _normalize_dataset_payload(data, ds_attrs)
    dataset = _create_dataset(target_group, dataset_key, payload)

    if ds_attrs:
        for attr_key, attr_val in ds_attrs.items():
            set_attr_safe(dataset, attr_key, attr_val)

    if "nameID" not in (ds_attrs or {}):
        set_attr_safe(dataset, "nameID", str(key))


def initialize_output_h5(
    h5file: h5py.File,
    *,
    holodoppler_source_file: str | None = None,
    doppler_vision_source_file: str | None = None,
) -> None:
    set_attr_safe(h5file, "eyeflow_version", _project_version())
    if holodoppler_source_file:
        h5file.attrs["holodoppler_source_file"] = holodoppler_source_file

        # Pass registration and zernike data through directly
        with open_h5(holodoppler_source_file, "r") as source_h5:
            for path in HD_OUTPUT_PASSTHROUGH_PATHS:
                if source_h5.get(path) is None:
                    continue
                source_h5.copy(path, h5file, name=path)
    if doppler_vision_source_file:
        h5file.attrs["doppler_vision_source_file"] = doppler_vision_source_file
    _initialize_app_versions(h5file, doppler_vision_source_file)
    primary_source = holodoppler_source_file or doppler_vision_source_file
    if primary_source:
        h5file.attrs["source_file"] = primary_source


def _initialize_app_versions(
    h5file: h5py.File,
    doppler_vision_source_file: str | None,
) -> None:
    """Write application versions as one scalar JSON dataset."""

    versions: dict[str, object] = {}
    if doppler_vision_source_file:
        with open_h5(doppler_vision_source_file, "r") as source_h5:
            source_versions = source_h5.get("app_versions")
            versions = _read_app_versions(source_versions)

    versions["EF_version"] = _project_version()
    if "app_versions" in h5file:
        del h5file["app_versions"]
    h5file.create_dataset(
        "app_versions",
        data=json.dumps(versions),
        dtype=h5py.string_dtype(encoding="utf-8"),
    )


def _read_app_versions(source_versions) -> dict[str, object]:
    if isinstance(source_versions, h5py.Dataset):
        value = _scalar_text(source_versions[()])
        if value is None:
            return {}
        try:
            parsed = json.loads(value)
        except json.JSONDecodeError:
            return {}
        return parsed if isinstance(parsed, dict) else {}

    # Accept the group representation emitted by older EyeFlow builds.
    if isinstance(source_versions, h5py.Group):
        versions: dict[str, object] = {}
        for key, value in source_versions.items():
            if isinstance(value, h5py.Dataset):
                scalar = _scalar_text(value[()])
                if scalar is not None:
                    versions[str(key)] = scalar
        return versions
    return {}


def _scalar_text(value) -> str | None:
    if isinstance(value, np.ndarray):
        if value.size != 1:
            return None
        value = value.reshape(-1)[0]
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value) if value is not None else None


def _project_version() -> str:
    pyproject_path = Path(__file__).resolve().parents[3] / "pyproject.toml"
    try:
        pyproject_text = pyproject_path.read_text(encoding="utf-8")
    except OSError:
        pyproject_text = ""

    match = re.search(r'(?m)^version\s*=\s*"([^"]+)"\s*$', pyproject_text)
    return match.group(1) if match else (app_version() or "unknown")


def _normalize_dataset_payload(data, ds_attrs):
    original_class = None
    payload = data

    if isinstance(payload, bool) or (
        isinstance(payload, np.ndarray) and payload.dtype == np.bool_
    ):
        payload = np.asarray(payload, dtype=np.uint8)
        original_class = "bool"
    elif isinstance(payload, (list, tuple)):
        payload = np.asarray(payload)

    payload = _downcast_numeric_payload(payload)

    if original_class is not None:
        ds_attrs = {} if ds_attrs is None else dict(ds_attrs)
        ds_attrs.setdefault("original_class", original_class)

    return payload, ds_attrs


def _downcast_numeric_payload(payload):
    if isinstance(payload, bool):
        return payload
    if isinstance(payload, float):
        return np.float32(payload)
    if isinstance(payload, int):
        return np.int32(payload)
    if isinstance(payload, complex):
        return np.complex64(payload)
    if not isinstance(payload, np.ndarray):
        return payload
    if payload.dtype.kind == "f":
        return payload.astype(np.float32, copy=False)
    if payload.dtype.kind == "c":
        return payload.astype(np.complex64, copy=False)
    if payload.dtype.kind == "i":
        return payload.astype(np.int32, copy=False)
    if payload.dtype.kind == "u":
        return payload.astype(np.uint32, copy=False)
    return payload


def _create_dataset(group: h5py.Group, dataset_key: str, payload):
    if isinstance(payload, str):
        return group.create_dataset(
            dataset_key,
            data=payload,
            dtype=h5py.string_dtype(encoding="utf-8"),
        )
    if isinstance(payload, (list, tuple)) and all(isinstance(item, str) for item in payload):
        return group.create_dataset(
            dataset_key,
            data=np.asarray(payload, dtype=object),
            dtype=h5py.string_dtype(encoding="utf-8"),
        )
    try:
        return group.create_dataset(dataset_key, data=payload)
    except (TypeError, ValueError):
        return _create_fallback_dataset(group, dataset_key, payload)


def _create_fallback_dataset(group: h5py.Group, dataset_key: str, payload):
    if isinstance(payload, np.ndarray) and payload.dtype.kind in {"U", "O"}:
        return group.create_dataset(
            dataset_key,
            data=np.asarray(payload, dtype=object),
            dtype=h5py.string_dtype(encoding="utf-8"),
        )
    return group.create_dataset(
        dataset_key,
        data=str(payload),
        dtype=h5py.string_dtype(encoding="utf-8"),
    )
