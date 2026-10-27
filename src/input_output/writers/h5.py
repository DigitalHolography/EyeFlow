"""Write EyeFlow runtime values into HDF5 files."""

import json
import re
from contextlib import contextmanager
from pathlib import Path

import h5py
import numpy as np

from app_settings import app_version
from ..payloads import dataset_parts
from ..schema.holodoppler import HD_OUTPUT_PASSTHROUGH_PATHS
from ..schema.provenance import inherited_provenance


# Small arrays save less payload than the per-dataset chunk index costs.
_MIN_COMPRESSION_BYTES = 4096
# The virtual mapping itself costs several KiB; use it only when the unique
# frequencies can repay that metadata overhead.
_MIN_MIRRORED_FFT_SAVINGS_BYTES = 32 * 1024


def profile_h5_options(shape: tuple[int, ...]) -> dict[str, object]:
    """Compress profiles in chunks aligned to one segment."""
    options: dict[str, object] = {
        "compression": "gzip", "compression_opts": 9, "shuffle": True,
    }
    if len(shape) not in (4, 5) or not all(shape):
        return options
    sample_count, time_count = shape[:2]
    target_elements = (1024 * 1024) // np.dtype(np.float32).itemsize
    if len(shape) == 4:
        time_chunk = min(time_count, 16, max(target_elements // sample_count, 1))
        options["chunks"] = (sample_count, time_chunk, 1, 1)
    else:
        time_chunk = min(time_count, max(target_elements // sample_count, 1))
        options["chunks"] = (sample_count, time_chunk, 1, 1, 1)
    return options


def velocity_map_h5_options(shape: tuple[int, ...]) -> dict[str, object]:
    """Compress per-beat maps while retaining segment-oriented chunks."""
    options: dict[str, object] = {
        "compression": "gzip", "compression_opts": 9, "shuffle": True,
    }
    if len(shape) in (6, 7) and all(shape):
        target_elements = (1024 * 1024) // np.dtype(np.float32).itemsize
        tail_count = shape[-1] if len(shape) == 7 else 1
        x_chunk = min(shape[0], 512)
        y_chunk = min(shape[1], max(target_elements // (x_chunk * tail_count), 1))
        time_chunk = min(
            shape[2], 4, max(target_elements // (x_chunk * y_chunk * tail_count), 1)
        )
        options["chunks"] = (
            (x_chunk, y_chunk, time_chunk, 1, 1, 1)
            if len(shape) == 6
            else (x_chunk, y_chunk, time_chunk, 1, 1, 1, shape[-1])
        )
    return options


def segment_mask_h5_options(shape: tuple[int, ...]) -> dict[str, object]:
    options: dict[str, object] = {
        "dtype": np.bool_, "compression": "gzip", "compression_opts": 9,
        "shuffle": True,
    }
    if len(shape) == 4 and all(shape):
        options["chunks"] = (shape[0], shape[1], 1, 1)
    return options


def normalize_h5_path(path: object) -> str:
    return str(path).replace("\\", "/").strip("/")


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


def initialize_output_h5(
    h5file: h5py.File,
    *,
    holodoppler_source_file: str | None = None,
    doppler_vision_source_file: str | None = None,
    holodoppler_source_h5: h5py.File | None = None,
    doppler_vision_source_h5: h5py.File | None = None,
    copy_hd_passthrough: bool = True,
) -> None:
    set_attr_safe(h5file, "eyeflow_version", _project_version())
    if holodoppler_source_file:
        h5file.attrs["holodoppler_source_file"] = holodoppler_source_file

        # Pass registration and zernike data through to their output schema paths.
        if copy_hd_passthrough:
            with _source_handle(holodoppler_source_h5, holodoppler_source_file) as source_h5:
                for source_path, output_path in HD_OUTPUT_PASSTHROUGH_PATHS:
                    if source_h5.get(source_path) is None:
                        continue
                    output_group, output_name = resolve_dataset_target(
                        h5file,
                        output_path,
                    )
                    if output_name in output_group:
                        del output_group[output_name]
                    _copy_source_object(source_h5[source_path], output_group, output_name)
    if doppler_vision_source_file:
        h5file.attrs["doppler_vision_source_file"] = doppler_vision_source_file
    _initialize_app_versions(h5file, doppler_vision_source_file, doppler_vision_source_h5)
    primary_source = holodoppler_source_file or doppler_vision_source_file
    if primary_source:
        h5file.attrs["source_file"] = primary_source


def _initialize_app_versions(
    h5file: h5py.File,
    doppler_vision_source_file: str | None,
    doppler_vision_source_h5: h5py.File | None = None,
) -> None:
    """Write application versions as one scalar JSON dataset."""

    versions: dict[str, object] = {}
    if doppler_vision_source_file:
        with _source_handle(doppler_vision_source_h5, doppler_vision_source_file) as source_h5:
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


@contextmanager
def _source_handle(existing: h5py.File | None, path: str):
    if existing is not None:
        yield existing
    else:
        with open_h5(path, "r") as source_h5:
            yield source_h5


def _copy_source_object(source, parent: h5py.Group, name: str) -> None:
    """Copy external metadata with lossless compression where HDF5 permits it."""
    if isinstance(source, h5py.Group):
        # Retain link semantics for exports that use aliases or external links.
        links = [source.get(key, getlink=True) for key in source.keys()]
        if any(not isinstance(link, h5py.HardLink) for link in links):
            source.file.copy(source, parent, name=name)
            return
        target = parent.create_group(name)
        for key, value in source.attrs.items():
            target.attrs[key] = value
        for child_name, child in source.items():
            _copy_source_object(child, target, child_name)
        return
    if not isinstance(source, h5py.Dataset):
        source.file.copy(source, parent, name=name)
        return
    if (
        source.ndim == 0 or source.dtype.kind not in "biufc"
        or source.size * source.dtype.itemsize < _MIN_COMPRESSION_BYTES
        or source.is_virtual
        or any(
            h5py.check_dtype(ref=source.attrs.get_id(key).dtype) is not None
            for key in source.attrs
        )
    ):
        source.file.copy(source, parent, name=name)
        return
    target = parent.create_dataset(
        name,
        shape=source.shape,
        dtype=source.dtype,
        maxshape=source.maxshape,
        chunks=source.chunks or True,
        compression="gzip",
        compression_opts=9,
        shuffle=True,
        fillvalue=source.fillvalue,
    )
    for selection in target.iter_chunks():
        target[selection] = source[selection]
    for key, value in source.attrs.items():
        target.attrs[key] = value


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
        payload = _array_from_sequence(payload)

    payload = _downcast_numeric_payload(payload)

    if original_class is not None:
        ds_attrs = {} if ds_attrs is None else dict(ds_attrs)
        ds_attrs.setdefault("original_class", original_class)

    return payload, ds_attrs


def _array_from_sequence(values: list | tuple) -> np.ndarray:
    # NumPy can silently coerce mixed signed/large-unsigned integers to float.
    objects = np.asarray(values, dtype=object)
    if objects.size and all(
        isinstance(item, (int, np.integer)) and not isinstance(item, (bool, np.bool_))
        for item in objects.flat
    ):
        minimum = min(int(item) for item in objects.flat)
        maximum = max(int(item) for item in objects.flat)
        if minimum < np.iinfo(np.int64).min or maximum > np.iinfo(np.uint64).max:
            raise OverflowError("Integer sequence is outside the HDF5 64-bit range.")
        if minimum < 0:
            if maximum > np.iinfo(np.int64).max:
                raise OverflowError(
                    "Integer sequence cannot fit in a single signed HDF5 64-bit dtype."
                )
            return np.asarray(values, dtype=np.int64)
        dtype = np.uint64 if maximum > np.iinfo(np.int64).max else np.int64
        return np.asarray(values, dtype=dtype)
    return np.asarray(values)


def _downcast_numeric_payload(payload):
    if isinstance(payload, bool):
        return payload
    if isinstance(payload, float):
        return np.float32(payload)
    if isinstance(payload, int):
        if np.iinfo(np.int64).min <= payload <= np.iinfo(np.int64).max:
            return np.int64(payload)
        if 0 <= payload <= np.iinfo(np.uint64).max:
            return np.uint64(payload)
        raise OverflowError(f"Integer value is outside the HDF5 64-bit range: {payload}")
    if isinstance(payload, complex):
        return np.complex64(payload)
    if not isinstance(payload, np.ndarray):
        return payload
    if payload.dtype.kind == "f":
        return payload.astype(np.float32, copy=False)
    if payload.dtype.kind == "c":
        return payload.astype(np.complex64, copy=False)
    if payload.dtype.kind in {"i", "u"}:
        # Preserve integer arrays rather than narrowing them without a range check.
        return payload
    return payload


def _create_dataset(
    group: h5py.Group,
    dataset_key: str,
    payload,
    *,
    h5_options: dict[str, object] | None = None,
):
    if isinstance(payload, str):
        return group.create_dataset(
            dataset_key,
            data=payload,
            dtype=h5py.string_dtype(encoding="utf-8"),
        )
    options = dict(h5_options or {})
    sparse_nan_chunks = bool(options.pop("sparse_nan_chunks", False))
    mirror_fft_axis = options.pop("mirror_fft_axis", None)
    # HDF5 filters require chunked arrays; scalars and strings stay unfiltered.
    if (
        isinstance(payload, np.ndarray)
        and payload.ndim > 0
        and payload.dtype.kind in {"b", "i", "u", "f", "c"}
    ):
        options.setdefault(
            "compression",
            "gzip" if payload.nbytes >= _MIN_COMPRESSION_BYTES else None,
        )
        if options["compression"] == "gzip":
            options.setdefault("compression_opts", 9)
        if options["compression"] is not None:
            options.setdefault("shuffle", True)
    if isinstance(payload, np.ndarray) and payload.dtype.kind in {"U", "O"}:
        if all(isinstance(item, (str, bytes)) for item in payload.flat):
            payload = np.asarray(payload, dtype=object)
            options.setdefault("dtype", h5py.string_dtype(encoding="utf-8"))
    try:
        if mirror_fft_axis is not None:
            return _create_mirrored_fft_dataset(
                group, dataset_key, payload, int(mirror_fft_axis), options,
                sparse_nan_chunks=sparse_nan_chunks,
            )
        if sparse_nan_chunks:
            return _create_sparse_nan_dataset(group, dataset_key, payload, options)
        return group.create_dataset(dataset_key, data=payload, **options)
    except (TypeError, ValueError) as exc:
        raise type(exc)(f"Could not write HDF5 dataset '{dataset_key}': {exc}") from exc


def _create_sparse_nan_dataset(group, name, payload, options):
    """Leave exact fill-value chunks unallocated while preserving dense reads."""
    if not isinstance(payload, np.ndarray) or payload.dtype != np.dtype(np.float32):
        raise TypeError("Sparse NaN storage requires a float32 array.")
    if not payload.size:
        return group.create_dataset(name, data=payload, **options)
    options = dict(options)
    options.setdefault("chunks", True)
    options["fillvalue"] = np.float32(np.nan)
    dataset = group.create_dataset(name, shape=payload.shape, dtype=payload.dtype, **options)
    fill_bits = np.asarray(dataset.fillvalue, dtype=np.float32).view(np.uint32).item()
    for selection in dataset.iter_chunks():
        chunk = payload[selection]
        if np.any(chunk.view(np.uint32) != fill_bits):
            dataset[selection] = chunk
    return dataset


def _create_mirrored_fft_dataset(
    group, name, payload, axis, options, *, sparse_nan_chunks,
):
    """Store unique real-signal FFT magnitudes behind the full public shape."""
    if not isinstance(payload, np.ndarray) or payload.dtype != np.dtype(np.float32):
        raise TypeError("Mirrored FFT storage requires a float32 array.")
    if axis < 0 or axis >= payload.ndim:
        raise ValueError("Mirrored FFT axis is out of range.")
    backing_name = f"__{name}_unique_frequency_storage"
    count = payload.shape[axis]
    unique_count = count // 2 + 1
    savings_bytes = payload.nbytes * (count - unique_count) // max(count, 1)
    if (
        count < 2 or savings_bytes < _MIN_MIRRORED_FFT_SAVINGS_BYTES
        or not _is_exactly_mirrored(payload, axis)
    ):
        if sparse_nan_chunks:
            return _create_sparse_nan_dataset(group, name, payload, options)
        return group.create_dataset(name, data=payload, **options)

    selection = [slice(None)] * payload.ndim
    selection[axis] = slice(0, unique_count)
    unique = payload[tuple(selection)]
    backing_options = dict(options)
    if backing_options.get("chunks") not in (None, True):
        chunks = list(backing_options["chunks"])
        chunks[axis] = min(chunks[axis], unique_count)
        backing_options["chunks"] = tuple(chunks)
    backing = (
        _create_sparse_nan_dataset(group, backing_name, unique, backing_options)
        if sparse_nan_chunks else group.create_dataset(backing_name, data=unique, **backing_options)
    )
    layout = h5py.VirtualLayout(shape=payload.shape, dtype=payload.dtype)
    source = h5py.VirtualSource(".", backing.name, shape=backing.shape)
    layout[tuple(selection)] = source
    for frequency in range(unique_count, count):
        destination = [slice(None)] * payload.ndim
        origin = [slice(None)] * payload.ndim
        destination[axis] = frequency
        origin[axis] = count - frequency
        layout[tuple(destination)] = source[tuple(origin)]
    return group.create_virtual_dataset(name, layout, fillvalue=np.float32(np.nan))


def _is_exactly_mirrored(values: np.ndarray, axis: int) -> bool:
    count = values.shape[axis]
    for frequency in range(count // 2 + 1, count):
        first = [slice(None)] * values.ndim
        second = [slice(None)] * values.ndim
        first[axis] = frequency
        second[axis] = count - frequency
        if not np.array_equal(
            values[tuple(first)].view(np.uint32),
            values[tuple(second)].view(np.uint32),
        ):
            return False
    return True
