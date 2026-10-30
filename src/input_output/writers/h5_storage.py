"""HDF5 payload conversion, compression, and dataset storage policy."""

import h5py
import numpy as np


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
    if hasattr(group, "store_dataset"):
        return group.store_dataset(dataset_key, payload, h5_options)
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
