"""Generic HDF5 file, dataset, and attribute write operations."""

from contextlib import contextmanager
from pathlib import Path

import h5py
import numpy as np

from ..schema.base import normalize_h5_path
from ..schema.output_payload import dataset_parts
from ..schema.provenance import inherited_provenance
from .h5_storage import (
    _array_from_sequence,
    _create_dataset,
    _downcast_numeric_payload,
    _normalize_dataset_payload,
)


def open_h5(path: Path | str, mode: str = "r") -> h5py.File:
    h5_path = Path(path)
    try:
        return h5py.File(h5_path, mode)
    except OSError as exc:
        raise OSError(f"Could not open HDF5 file '{h5_path}': {exc}") from exc


@contextmanager
def h5_write_checkpoint(h5file: h5py.File):
    """Undo a failed pipeline's writes without copying scientific arrays.

    Open object handles retain replaced datasets until the checkpoint ends.
    Pipeline writers replace datasets; in-place array edits are not supported.
    """
    objects = {"": h5file}

    def retain(name, obj):
        objects[name] = obj

    h5file.visititems(retain)
    attributes = {
        name: {key: (obj.attrs[key], obj.attrs.get_id(key).dtype) for key in obj.attrs}
        for name, obj in objects.items()
    }
    try:
        yield
    except Exception:
        current_paths = []
        h5file.visit(current_paths.append)
        for name in sorted(current_paths, key=lambda path: path.count("/"), reverse=True):
            if name not in objects:
                del h5file[name]
        for name, obj in objects.items():
            if name:
                current = h5file.get(name)
                if current is None or current.id != obj.id:
                    if current is not None:
                        del h5file[name]
                    h5file[name] = obj
            for key in list(obj.attrs):
                del obj.attrs[key]
            for key, (value, dtype) in attributes[name].items():
                obj.attrs.create(key, value, dtype=dtype)
        h5file.flush()
        raise


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
    data = _normalize_attribute_payload(value)
    try:
        h5obj.attrs[key] = data
    except (TypeError, ValueError) as exc:
        raise type(exc)(f"Could not write HDF5 attribute '{key}': {exc}") from exc


def write_root_attributes(h5file: h5py.File, values: dict[str, object]) -> None:
    """Apply engine bookkeeping to the root output object."""
    for key, value in values.items():
        h5file.attrs[key] = value


def delete_output_path(h5file: h5py.File, path: str) -> None:
    normalized = normalize_h5_path(path)
    if normalized in h5file:
        del h5file[normalized]


def _normalize_attribute_payload(value):
    data = value
    if isinstance(value, (list, tuple)):
        if all(isinstance(item, str) for item in value):
            data = np.asarray(value, dtype=h5py.string_dtype(encoding="utf-8"))
        else:
            data = _array_from_sequence(value)
    return _downcast_numeric_payload(data)


def _matches_inherited_attribute(value, inherited) -> bool:
    actual = np.asarray(_normalize_attribute_payload(value))
    expected = np.asarray(inherited)
    if actual.shape != expected.shape or actual.dtype != expected.dtype:
        return False
    if actual.dtype.kind in {"b", "i", "u", "f", "c"}:
        return actual.tobytes() == expected.tobytes()
    return np.array_equal(actual, expected)


def write_value_dataset(group: h5py.Group, key: str, value) -> None:
    data, ds_attrs, h5_options = dataset_parts(value)

    target_group, dataset_key = resolve_dataset_target(group, str(key))
    if dataset_key in target_group:
        del target_group[dataset_key]
    backing_name = f"__{dataset_key}_unique_frequency_storage"
    if backing_name in target_group:
        del target_group[backing_name]

    payload, ds_attrs = _normalize_dataset_payload(data, ds_attrs)
    dataset = _create_dataset(
        target_group,
        dataset_key,
        payload,
        h5_options=h5_options,
    )

    if ds_attrs:
        inherited_attrs = inherited_provenance(target_group)
        for attr_key, attr_val in ds_attrs.items():
            if attr_key in inherited_attrs:
                if _matches_inherited_attribute(attr_val, inherited_attrs[attr_key]):
                    continue
            set_attr_safe(dataset, attr_key, attr_val)
