"""Tests for HDF5 output writer metadata."""

from __future__ import annotations

import json
import re
import sys
import tempfile
import unittest
from pathlib import Path

import h5py

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from input_output.writers.h5 import initialize_output_h5  # noqa: E402


class H5WriterTests(unittest.TestCase):
    def test_initialize_output_h5_writes_eyeflow_version(self) -> None:
        version = _pyproject_version()
        with tempfile.TemporaryDirectory() as tmp_dir:
            output_path = Path(tmp_dir) / "output.h5"
            with h5py.File(output_path, "w") as h5file:
                initialize_output_h5(h5file)

            with h5py.File(output_path, "r") as h5file:
                self.assertEqual(version, h5file.attrs["eyeflow_version"])
                self.assertEqual(
                    {"EF_version": version},
                    json.loads(_read_scalar_string(h5file["app_versions"])),
                )

    def test_initialize_output_h5_copies_dopplerview_app_versions(self) -> None:
        version = _pyproject_version()
        with tempfile.TemporaryDirectory() as tmp_dir:
            source_path = Path(tmp_dir) / "input_DV.h5"
            output_path = Path(tmp_dir) / "output.h5"
            with h5py.File(source_path, "w") as source_h5:
                source_h5.create_dataset(
                    "app_versions",
                    data=json.dumps(
                        {
                            "HD_version": "v0.5.0",
                            "DV_version": "v1.17.2",
                        }
                    ),
                    dtype=h5py.string_dtype(encoding="utf-8"),
                )

            with h5py.File(output_path, "w") as output_h5:
                initialize_output_h5(
                    output_h5,
                    doppler_vision_source_file=str(source_path),
                )

            with h5py.File(output_path, "r") as output_h5:
                self.assertEqual(
                    {
                        "HD_version": "v0.5.0",
                        "DV_version": "v1.17.2",
                        "EF_version": version,
                    },
                    json.loads(_read_scalar_string(output_h5["app_versions"])),
                )
                self.assertEqual((), output_h5["app_versions"].shape)


def _pyproject_version() -> str:
    pyproject_path = Path(__file__).resolve().parents[1] / "pyproject.toml"
    text = pyproject_path.read_text(encoding="utf-8")
    match = re.search(r'(?m)^version\s*=\s*"([^"]+)"\s*$', text)
    if match is None:
        raise AssertionError(f"Missing project version in {pyproject_path}")
    return match.group(1)


def _read_scalar_string(dataset: h5py.Dataset) -> str:
    value = dataset[()]
    return value.decode("utf-8") if isinstance(value, bytes) else str(value)


if __name__ == "__main__":
    unittest.main()
