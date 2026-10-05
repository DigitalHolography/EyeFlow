"""Tests for HDF5 output writer metadata."""

from __future__ import annotations

import json
import re
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import h5py
import numpy as np

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from input_output.writers.h5 import (  # noqa: E402
    initialize_output_h5,
    set_attr_safe,
    write_value_dataset,
)
from app_settings import app_version  # noqa: E402
from pipeline_engine.base import DatasetValue  # noqa: E402


class H5WriterTests(unittest.TestCase):
    def test_output_version_uses_real_package_version_not_environment_override(self) -> None:
        with patch.dict("os.environ", {"EYEFLOW_VERSION": "9.9-test"}):
            with h5py.File("version_writer_test.h5", "w", driver="core", backing_store=False) as h5file:
                initialize_output_h5(h5file, eyeflow_version=app_version() or "unknown")
                self.assertEqual(_pyproject_version(), app_version())
                self.assertEqual(_pyproject_version(), h5file.attrs["eyeflow_version"])
                self.assertEqual(
                    _pyproject_version(),
                    json.loads(_read_scalar_string(h5file["app_versions"]))["EF_version"],
                )

    def test_integer_values_are_not_narrowed_out_of_range(self) -> None:
        with h5py.File("integer_writer_test.h5", "w", driver="core", backing_store=False) as h5file:
            signed = np.array([2**40, -(2**40)], dtype=np.int64)
            unsigned = np.array([2**63 + 1], dtype=np.uint64)
            write_value_dataset(h5file, "signed", signed)
            write_value_dataset(h5file, "unsigned", unsigned)
            write_value_dataset(h5file, "large_scalar", 2**63 + 1)
            write_value_dataset(h5file, "unsigned_list", [0, 2**63 + 1])
            set_attr_safe(h5file, "large_integer", 2**40)

            np.testing.assert_array_equal(h5file["signed"][()], signed)
            np.testing.assert_array_equal(h5file["unsigned"][()], unsigned)
            self.assertEqual(2**63 + 1, h5file["large_scalar"][()])
            np.testing.assert_array_equal(
                h5file["unsigned_list"][()], np.array([0, 2**63 + 1], dtype=np.uint64)
            )
            self.assertEqual(2**40, h5file.attrs["large_integer"])

    def test_integer_outside_hdf5_range_raises(self) -> None:
        with h5py.File("integer_writer_test.h5", "w", driver="core", backing_store=False) as h5file:
            with self.assertRaises(OverflowError):
                write_value_dataset(h5file, "too_large", 2**64)
            with self.assertRaises(OverflowError):
                write_value_dataset(h5file, "mixed_range", [-1, 2**63 + 1])
            self.assertNotIn("too_large", h5file)
            self.assertNotIn("mixed_range", h5file)

    def test_invalid_h5_options_do_not_stringify_numeric_data(self) -> None:
        with h5py.File("options_writer_test.h5", "w", driver="core", backing_store=False) as h5file:
            value = DatasetValue(np.array([1, 2]), h5_options={"chunks": (3,)})
            with self.assertRaisesRegex(ValueError, "Could not write HDF5 dataset 'values'"):
                write_value_dataset(h5file, "values", value)
            self.assertNotIn("values", h5file)

    def test_invalid_attribute_does_not_become_a_string(self) -> None:
        with h5py.File("attribute_writer_test.h5", "w", driver="core", backing_store=False) as h5file:
            with self.assertRaisesRegex(TypeError, "Could not write HDF5 attribute 'bad'"):
                set_attr_safe(h5file, "bad", {"unexpected": "mapping"})
            self.assertNotIn("bad", h5file.attrs)

    def test_unicode_array_is_written_as_string_array(self) -> None:
        with h5py.File("strings_writer_test.h5", "w", driver="core", backing_store=False) as h5file:
            write_value_dataset(h5file, "labels", np.array(["artery", "vein"]))
            self.assertEqual([b"artery", b"vein"], list(h5file["labels"][()]))

    def test_initialize_output_h5_writes_eyeflow_version(self) -> None:
        version = _pyproject_version()
        with tempfile.TemporaryDirectory() as tmp_dir:
            output_path = Path(tmp_dir) / "output.h5"
            with h5py.File(output_path, "w") as h5file:
                initialize_output_h5(h5file, eyeflow_version=app_version() or "unknown")

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
                    eyeflow_version=app_version() or "unknown",
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

    def test_initialize_output_h5_places_registration_under_meta(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            source_path = Path(tmp_dir) / "input_HD.h5"
            output_path = Path(tmp_dir) / "output.h5"
            registration = np.asarray(
                [[1.0, 2.0], [3.0, 4.0]],
                dtype=np.float32,
            )
            with h5py.File(source_path, "w") as source_h5:
                dataset = source_h5.create_dataset(
                    "registration",
                    data=registration,
                )
                dataset.attrs["unit"] = "pixel"
                source_h5.create_dataset(
                    "zernike_coefs_radians",
                    data=np.asarray([0.1, 0.2], dtype=np.float32),
                )

            with h5py.File(output_path, "w") as output_h5:
                initialize_output_h5(
                    output_h5,
                    eyeflow_version=app_version() or "unknown",
                    holodoppler_source_file=str(source_path),
                )

            with h5py.File(output_path, "r") as output_h5:
                self.assertNotIn("registration", output_h5)
                np.testing.assert_array_equal(
                    output_h5["Meta/registration"][...],
                    registration,
                )
                self.assertEqual(
                    "pixel",
                    output_h5["Meta/registration"].attrs["unit"],
                )
                self.assertIn("zernike_coefs_radians", output_h5)


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
