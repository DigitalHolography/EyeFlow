"""Initialize EyeFlow HDF5 files and copy selected source metadata."""

import json
import re
from contextlib import contextmanager
from pathlib import Path

import h5py
import numpy as np

from app_settings import app_version
from ..schema.holodoppler import HD_OUTPUT_PASSTHROUGH_PATHS
from .h5 import open_h5, resolve_dataset_target, set_attr_safe
from .h5_storage import _MIN_COMPRESSION_BYTES


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


