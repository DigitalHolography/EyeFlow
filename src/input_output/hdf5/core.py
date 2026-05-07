"""Read and write EyeFlow runtime HDF5 files.

This file contains concrete HDF5 operations used by the runtime output manager.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import h5py
import numpy as np

if TYPE_CHECKING:
    from pipeline_engine import DatasetValue


def open_h5(path: Path | str, mode: str = "r") -> h5py.File:
    return h5py.File(Path(path), mode)


def normalize_h5_path(path: object) -> str:
    return str(path).replace("\\", "/").strip("/")


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
            f"Cannot create subgroup '{part}' for key '{key}': a dataset already exists at that path."
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


def _get_dataset_creation_kwargs(payload: np.ndarray) -> dict[str, object]:
    return {}


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

    if isinstance(payload, str):
        dataset = target_group.create_dataset(
            dataset_key, data=payload, dtype=h5py.string_dtype(encoding="utf-8")
        )
    elif isinstance(payload, (list, tuple)) and all(isinstance(item, str) for item in payload):
        dataset = target_group.create_dataset(
            dataset_key,
            data=np.asarray(payload, dtype=object),
            dtype=h5py.string_dtype(encoding="utf-8"),
        )
    else:
        try:
            create_kwargs = {}
            if isinstance(payload, np.ndarray):
                create_kwargs = _get_dataset_creation_kwargs(payload)
            dataset = target_group.create_dataset(dataset_key, data=payload, **create_kwargs)
        except (TypeError, ValueError):
            if isinstance(payload, np.ndarray) and payload.dtype.kind in {"U", "O"}:
                dataset = target_group.create_dataset(
                    dataset_key,
                    data=np.asarray(payload, dtype=object),
                    dtype=h5py.string_dtype(encoding="utf-8"),
                )
            else:
                dataset = target_group.create_dataset(
                    dataset_key, data=str(payload), dtype=h5py.string_dtype(encoding="utf-8")
                )

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
    if holodoppler_source_file:
        h5file.attrs["holodoppler_source_file"] = holodoppler_source_file
    if doppler_vision_source_file:
        h5file.attrs["doppler_vision_source_file"] = doppler_vision_source_file
    primary_source = holodoppler_source_file or doppler_vision_source_file
    if primary_source:
        h5file.attrs["source_file"] = primary_source

