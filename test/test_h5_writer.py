"""Tests for HDF5 output writer metadata."""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from unittest.mock import patch
from pathlib import Path

import h5py
import numpy as np

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from input_output.writers.h5 import (  # noqa: E402
    set_attr_safe,
    write_value_dataset,
)
from input_output.readers.h5 import PipelineInputSource, RawH5SourceReader  # noqa: E402
from input_output.writers.h5_initialization import initialize_output_h5  # noqa: E402
from input_output.writers.h5_output import PipelineH5Output  # noqa: E402
from input_output.schema import EyeFlowOutputPaths  # noqa: E402
from pipeline_engine.base import DatasetValue  # noqa: E402


class H5WriterTests(unittest.TestCase):
    def test_spatial_gradient_paths_and_references_follow_workflow(self) -> None:
        schema = EyeFlowOutputPaths.active()
        profile = f"{schema.spatial_gradient_profiles_root}/Artery/Transverse/value"
        metric = f"{schema.spatial_gradient_metrics_root}/Artery/Transverse/value"
        flow = schema.blood_volume_rate.artery.dynamic_edges
        with h5py.File("gradient_paths.h5", "w", driver="core", backing_store=False) as h5file:
            for root, branch, value in (
                ("ProcessingMoments", "Moments", 1),
                ("ProcessingBandRatio", "BandRatio", 2),
            ):
                output = PipelineH5Output(h5file, processing_root=root)
                output.write_many({
                    profile: DatasetValue(np.asarray(value)),
                    metric: DatasetValue(
                        np.asarray(value),
                        attrs={"source_profile": f"/{profile}"},
                    ),
                    flow: DatasetValue(
                        np.asarray(value),
                        attrs={"source_left_edge_index": f"/{metric}"},
                    ),
                })
                self.assertEqual(value, output.read(profile))
                self.assertEqual(value, output.read(metric))
                published = f"ProcessingWidth/SpatialGradient/{branch}"
                expected_profile_reference = (
                    f"/{profile}"
                    if branch == "Moments"
                    else f"{published}/Profiles/Artery/Transverse/value"
                )
                self.assertEqual(
                    expected_profile_reference,
                    h5file[f"{published}/Metrics/Artery/Transverse/value"].attrs["source_profile"],
                )
                expected_metric_reference = (
                    f"/{metric}"
                    if branch == "Moments"
                    else f"{published}/Metrics/Artery/Transverse/value"
                )
                self.assertEqual(
                    expected_metric_reference,
                    h5file[f"{root}/BloodVolumeRate/Artery/dynamicEdges/value"].attrs[
                        "source_left_edge_index"
                    ],
                )
            self.assertNotIn("ProcessingMoments/SpatialGradientProfiles", h5file)
            self.assertNotIn("ProcessingBandRatio/SpatialGradientProfiles", h5file)

    def test_runtime_typed_sources_read_small_fields_once(self) -> None:
        with h5py.File("cached_inputs.h5", "w", driver="core", backing_store=False) as h5file:
            h5file.create_dataset("HD_parameters", data=json.dumps({"sampling_freq": 1000}))
            h5file.create_dataset("segmentation/Retina/artery_mask", data=np.ones((2, 3)))
            reader = RawH5SourceReader(h5file=h5file, label="source")
            inputs = PipelineInputSource(reader, {})
            with patch.object(reader, "value", wraps=reader.value) as read_value:
                hd = inputs.as_holodoppler()
                self.assertIs(hd, inputs.as_holodoppler())
                self.assertEqual(1000, hd.sampling_frequency_hz())
                self.assertEqual(1000, hd.sampling_frequency_hz())
                self.assertEqual(1, read_value.call_count)
            with patch.object(reader, "array", wraps=reader.array) as read_array:
                dv = inputs.as_dopplerview()
                self.assertIs(dv, inputs.as_dopplerview())
                dv.retinal_mask("segmentation/Retina/artery_mask")
                dv.retinal_artery_mask()
                self.assertEqual(1, read_array.call_count)

    def test_default_compression_preserves_float32_bits(self) -> None:
        values = np.tile(
            np.array([0.0, -0.0, np.nan, np.inf, -np.inf, 1.234567], dtype=np.float32),
            256,
        )
        with tempfile.TemporaryDirectory() as tmp_dir:
            path = Path(tmp_dir) / "compressed.h5"
            with h5py.File(path, "w") as h5file:
                write_value_dataset(h5file, "values", values)

            with h5py.File(path, "r") as h5file:
                dataset = h5file["values"]
                self.assertEqual("gzip", dataset.compression)
                self.assertEqual(9, dataset.compression_opts)
                self.assertTrue(dataset.shuffle)
                self.assertEqual(values.dtype, dataset.dtype)
                self.assertEqual(values.shape, dataset.shape)
                self.assertEqual(values.tobytes(), dataset[()].tobytes())

    def test_sparse_nan_chunks_preserve_dense_values_without_allocating_empty_chunks(self) -> None:
        values = np.full((16, 16, 4, 2), np.nan, dtype=np.float32)
        values[:, :, :2, 0] = np.arange(16 * 16 * 2, dtype=np.float32).reshape(16, 16, 2)
        with h5py.File("sparse_nan.h5", "w", driver="core", backing_store=False) as h5file:
            write_value_dataset(h5file, "values", DatasetValue(
                values,
                h5_options={
                    "compression": "gzip", "shuffle": True,
                    "chunks": (16, 16, 2, 1), "sparse_nan_chunks": True,
                },
            ))
            dataset = h5file["values"]
            self.assertEqual(1, dataset.id.get_num_chunks())
            self.assertEqual(values.tobytes(), dataset[()].tobytes())

    def test_mirrored_fft_virtual_dataset_reads_full_exact_spectrum(self) -> None:
        values = np.abs(np.fft.fft(
            np.arange(181 * 128, dtype=np.float32).reshape(181, 128), axis=1,
        )).astype(np.float32)
        with tempfile.TemporaryDirectory() as tmp_dir:
            path = Path(tmp_dir) / "mirrored_fft.h5"
            with h5py.File(path, "w") as h5file:
                write_value_dataset(h5file, "profiles/fft", DatasetValue(
                    values,
                    attrs={"unit": "a.u."},
                    h5_options={
                        "compression": "gzip", "shuffle": True,
                        "chunks": (181, 128), "mirror_fft_axis": 1,
                        "sparse_nan_chunks": True,
                    },
                ))
            relocated = path.with_name("relocated.h5")
            path.rename(relocated)
            with h5py.File(relocated, "r") as h5file:
                public = h5file["profiles/fft"]
                backing = h5file["profiles/__fft_unique_frequency_storage"]
                self.assertTrue(public.is_virtual)
                self.assertEqual(values.shape, public.shape)
                self.assertEqual(values.tobytes(), public[()].tobytes())
                self.assertEqual(values.shape[:1] + (65,), backing.shape)
                self.assertEqual("gzip", backing.compression)
                self.assertEqual("a.u.", public.attrs["unit"])

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

    def test_initialize_uses_existing_handles_and_compresses_large_passthrough(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            hd_path = Path(tmp_dir) / "input_HD.h5"
            dv_path = Path(tmp_dir) / "input_DV.h5"
            output_path = Path(tmp_dir) / "output.h5"
            values = np.zeros((64, 64), dtype=np.float32)
            with h5py.File(hd_path, "w") as hd:
                hd.create_dataset("registration", data=values).attrs["unit"] = "pixel"
            with h5py.File(dv_path, "w") as dv:
                dv.create_dataset("app_versions", data=json.dumps({"DV_version": "test"}))
                laterality = dv.create_group("laterality")
                laterality.create_dataset("eye_laterality", data="right")
                confidence = laterality.create_dataset(
                    "eye_laterality_confidence",
                    data=np.float32(0.875),
                )
                confidence.attrs["unit"] = "1"
            with h5py.File(hd_path, "r") as hd, h5py.File(dv_path, "r") as dv:
                with h5py.File(output_path, "w") as output_h5:
                    with patch(
                        "input_output.writers.h5_initialization.open_h5",
                        side_effect=AssertionError("input was reopened"),
                    ):
                        initialize_output_h5(
                            output_h5,
                            holodoppler_source_file=str(hd_path),
                            doppler_vision_source_file=str(dv_path),
                            holodoppler_source_h5=hd,
                            doppler_vision_source_h5=dv,
                        )
            with h5py.File(output_path, "r") as output_h5:
                registration = output_h5["Meta/registration"]
                self.assertEqual("gzip", registration.compression)
                self.assertEqual(9, registration.compression_opts)
                self.assertEqual("pixel", registration.attrs["unit"])
                np.testing.assert_array_equal(values, registration[()])
                self.assertEqual(
                    "test", json.loads(_read_scalar_string(output_h5["app_versions"]))["DV_version"],
                )
                self.assertEqual(
                    "right",
                    _read_scalar_string(
                        output_h5["Segmentation/Laterality/EyeLaterality"]
                    ),
                )
                self.assertAlmostEqual(
                    0.875,
                    float(
                        output_h5[
                            "Segmentation/Laterality/EyeLateralityConfidence"
                        ][()]
                    ),
                )
                self.assertEqual(
                    "1",
                    output_h5[
                        "Segmentation/Laterality/EyeLateralityConfidence"
                    ].attrs["unit"],
                )

def _read_scalar_string(dataset: h5py.Dataset) -> str:
    value = dataset[()]
    return value.decode("utf-8") if isinstance(value, bytes) else str(value)


if __name__ == "__main__":
    unittest.main()
