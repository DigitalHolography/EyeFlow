"""Tests for pipeline runtime context helpers."""

from __future__ import annotations

import unittest

import h5py
import numpy as np

from pipeline_engine import ExecutionVariant, PipelineContext
from pipeline_engine.context import RawH5SourceReader
from input_output.h5_access import PipelineH5Output, PipelineInputSource
from input_output.schema.base import TypedSource
from input_output.schema.provenance import (
    INHERITED_PROVENANCE_KEYS_ATTR,
    read_dataset_attrs,
)
from input_output.writers.h5 import h5_write_checkpoint, set_attr_safe
from utils.logger import Logger


class PipelineContextTests(unittest.TestCase):
    def test_missing_arrays_return_none_for_explicit_none_default(self) -> None:
        with h5py.File("context_missing_test.h5", "w", driver="core", backing_store=False) as h5file:
            source = RawH5SourceReader(h5file=h5file, label="HD")
            wrapped = PipelineInputSource(source, {})
            output = PipelineH5Output(h5file)
            for reader in (source, wrapped, output):
                self.assertIsNone(reader.array("missing", dtype=np.float32, default=None))
                with self.assertRaises(KeyError):
                    reader.array("missing")

            typed = TypedSource(source)
            with self.assertRaises(KeyError):
                typed._array("missing")
            self.assertIsNone(typed._array("missing", default=None))

    def test_output_facade_writes_dataset_attributes(self) -> None:
        with h5py.File("context_output_test.h5", "w", driver="core", backing_store=False) as h5file:
            output = PipelineH5Output(h5file)
            output.write("result", np.array([1, 2]), unit="pixel")

            np.testing.assert_array_equal(h5file["result"][()], [1, 2])
            self.assertEqual("pixel", h5file["result"].attrs["unit"])

    def test_workflow_provenance_is_shared_with_local_overrides(self) -> None:
        with h5py.File("provenance.h5", "w", driver="core", backing_store=False) as h5file:
            h5file.attrs["band_ratio_frequency_scale_hz"] = 999.0
            for root, method, scale in [
                ("ProcessingMoments", "doppler_moments", 1.0),
                ("ProcessingBandRatio", "frequency_bands", 3.0),
            ]:
                output = PipelineH5Output(h5file, processing_root=root, provenance={
                    "velocity_estimation_method": method,
                    "velocity_unit": "mm/s",
                    "band_ratio_frequency_scale_hz": scale,
                })
                output.write("ProcessingMoments/signal", np.arange(4),
                             velocity_estimation_method=method,
                             band_ratio_frequency_scale_hz=scale,
                             unit="Hz", dimDesc=["time"], definition="Keep this definition",
                             source_path="ProcessingMoments/source")
                dataset = h5file[f"{root}/signal"]
                self.assertNotIn("velocity_estimation_method", dataset.attrs)
                self.assertNotIn("band_ratio_frequency_scale_hz", dataset.attrs)
                effective = output.read_attrs("ProcessingMoments/signal")
                self.assertEqual(method, effective["velocity_estimation_method"])
                self.assertEqual(scale, effective["band_ratio_frequency_scale_hz"])
                self.assertEqual("mm/s", effective["velocity_unit"])
                self.assertEqual("Hz", dataset.attrs["unit"])
                self.assertEqual(["time"], list(dataset.attrs["dimDesc"]))
                self.assertEqual("Keep this definition", dataset.attrs["definition"])
                self.assertEqual(f"{root}/source", dataset.attrs["source_path"])

                # Explicit dataset overrides remain local, including dtype changes.
                output.write("ProcessingMoments/override", [1],
                             band_ratio_frequency_scale_hz=scale + 1)
                self.assertEqual(scale + 1, output.read_attrs("ProcessingMoments/override")[
                    "band_ratio_frequency_scale_hz"])
                output.write("ProcessingMoments/typed", [1],
                             band_ratio_frequency_scale_hz=np.int32(scale))
                self.assertEqual(np.dtype("int32"),
                                 h5file[f"{root}/typed"].attrs.get_id(
                                     "band_ratio_frequency_scale_hz").dtype)

                # Another facade can resolve and reuse the persisted declaration.
                reopened = PipelineH5Output(h5file, processing_root=root)
                reopened.write("ProcessingMoments/later", [1],
                               velocity_estimation_method=method)
                self.assertNotIn("velocity_estimation_method", h5file[f"{root}/later"].attrs)
                self.assertEqual(method, reopened.read_attrs("ProcessingMoments/later")[
                    "velocity_estimation_method"])

                output.write("Segmentation/shared", [1], unit="pixel")
                self.assertNotIn("velocity_estimation_method",
                                 read_dataset_attrs(h5file["Segmentation/shared"]))
            self.assertEqual("doppler_moments", read_dataset_attrs(
                h5file["ProcessingMoments/signal"])["velocity_estimation_method"])
            self.assertEqual("frequency_bands", read_dataset_attrs(
                h5file["ProcessingBandRatio/signal"])["velocity_estimation_method"])

    def test_inherited_attributes_preserve_legacy_and_nearest_group_values(self) -> None:
        with h5py.File("legacy_provenance.h5", "w", driver="core", backing_store=False) as f:
            root = f.create_group("ProcessingMoments")
            root.attrs["velocity_estimation_method"] = "doppler_moments"
            legacy = root.create_dataset("legacy", data=[1])
            legacy.attrs["velocity_estimation_method"] = "frequency_bands"
            self.assertEqual(dict(legacy.attrs), read_dataset_attrs(legacy))
            missing = root.create_dataset("missing", data=[1])
            self.assertEqual({}, read_dataset_attrs(missing))
            set_attr_safe(root, INHERITED_PROVENANCE_KEYS_ATTR, ["velocity_estimation_method"])
            child = root.create_group("child")
            child.attrs["velocity_estimation_method"] = "frequency_bands"
            child.attrs["unit"] = "do not inherit"
            set_attr_safe(child, INHERITED_PROVENANCE_KEYS_ATTR, ["velocity_estimation_method"])
            dataset = child.create_dataset("value", data=[1])
            self.assertEqual({"velocity_estimation_method": "frequency_bands"},
                             read_dataset_attrs(dataset))
            dataset.attrs["velocity_estimation_method"] = "doppler_moments"
            self.assertEqual("doppler_moments",
                             read_dataset_attrs(dataset)["velocity_estimation_method"])
            output = PipelineH5Output(f, processing_root="ProcessingMoments")
            output.write("ProcessingMoments/child/written", [1],
                         velocity_estimation_method="doppler_moments")
            self.assertEqual("doppler_moments", output.read_attrs(
                "ProcessingMoments/child/written")["velocity_estimation_method"])
            self.assertIn("velocity_estimation_method", f[
                "ProcessingMoments/child/written"].attrs)
            with self.assertRaises(KeyError):
                PipelineH5Output(f).read_attrs("absent")

    def test_checkpoint_restores_shared_provenance_and_dataset_overrides(self) -> None:
        with h5py.File("checkpoint_provenance.h5", "w", driver="core", backing_store=False) as f:
            output = PipelineH5Output(f, processing_root="ProcessingMoments",
                                      provenance={"velocity_unit": "mm/s"})
            output.write("ProcessingMoments/signal", [1], velocity_unit="local")
            with self.assertRaisesRegex(ValueError, "injected"):
                with h5_write_checkpoint(f):
                    PipelineH5Output(f, processing_root="ProcessingMoments", provenance={
                        "velocity_unit": "invalid", "new_provenance": "remove on rollback",
                    }).write("ProcessingMoments/signal", [2])
                    raise ValueError("injected")
            self.assertEqual("mm/s", f["ProcessingMoments"].attrs["velocity_unit"])
            self.assertEqual(["velocity_unit"], list(
                f["ProcessingMoments"].attrs[INHERITED_PROVENANCE_KEYS_ATTR]))
            self.assertNotIn("new_provenance", output.read_attrs("ProcessingMoments/signal"))
            self.assertEqual("local", output.read_attrs("ProcessingMoments/signal")["velocity_unit"])
            np.testing.assert_array_equal([1], output.read("ProcessingMoments/signal"))

    def test_log_emits_to_callback(self) -> None:
        messages: list[str] = []
        Logger.configure(on_log=messages.append)
        with h5py.File(
            "context_log_test.h5",
            "w",
            driver="core",
            backing_store=False,
        ) as h5file:
            ctx = PipelineContext(
                work_h5=h5file,
                holodoppler_h5=None,
                doppler_vision_h5=None,
            )

            ctx.log("Starting test pipeline...")

        Logger.reset_current()
        self.assertEqual(["Starting test pipeline..."], messages)

    def test_log_without_callback_is_noop(self) -> None:
        Logger.reset_current()
        with h5py.File(
            "context_log_test.h5",
            "w",
            driver="core",
            backing_store=False,
        ) as h5file:
            ctx = PipelineContext(
                work_h5=h5file,
                holodoppler_h5=None,
                doppler_vision_h5=None,
            )

            ctx.log("No listener")

    def test_pipeline_options_and_schedule_are_available_to_runners(self) -> None:
        with h5py.File(
            "context_options_test.h5",
            "w",
            driver="core",
            backing_store=False,
        ) as h5file:
            ctx = PipelineContext(
                work_h5=h5file,
                holodoppler_h5=None,
                doppler_vision_h5=None,
                pipeline_name="velocity_analysis",
                pipeline_options={
                    "velocity_analysis": ("segments", "quadrants"),
                    "waveform_shape_metrics": (),
                },
                pipeline_order=(
                    "velocity",
                    "velocity_analysis",
                ),
                pipeline_targets=("velocity_analysis",),
            )

            self.assertTrue(ctx.option_enabled("segments"))
            self.assertTrue(
                ctx.option_enabled("quadrants", pipeline="velocity_analysis")
            )
            self.assertFalse(
                ctx.option_enabled("quadrants", pipeline="waveform_shape_metrics")
            )
            self.assertEqual(
                frozenset({"segments", "quadrants"}),
                ctx.options_for("velocity_analysis"),
            )
            self.assertTrue(ctx.pipeline_scheduled("velocity"))
            self.assertFalse(ctx.pipeline_scheduled("pdf_report"))
            self.assertTrue(ctx.pipeline_targeted("velocity_analysis"))
            self.assertFalse(ctx.pipeline_targeted("velocity"))

    def test_velocity_estimation_method_is_available_to_runners(self) -> None:
        with h5py.File(
            "context_velocity_method_test.h5",
            "w",
            driver="core",
            backing_store=False,
        ) as h5file:
            default_ctx = PipelineContext(
                work_h5=h5file,
                holodoppler_h5=None,
                doppler_vision_h5=None,
            )
            frequency_band_ctx = PipelineContext(
                work_h5=h5file,
                holodoppler_h5=None,
                doppler_vision_h5=None,
                velocity_estimation_method="frequency_bands",
                band_ratio_frequency_scale_hz=2.0,
            )

            self.assertEqual(
                "doppler_moments",
                default_ctx.velocity_estimation_method,
            )
            self.assertEqual(
                "frequency_bands",
                frequency_band_ctx.velocity_estimation_method,
            )
            self.assertEqual(
                2.0,
                frequency_band_ctx.band_ratio_frequency_scale_hz,
            )

            with self.assertRaisesRegex(ValueError, "velocity_estimation_method"):
                PipelineContext(
                    work_h5=h5file,
                    holodoppler_h5=None,
                    doppler_vision_h5=None,
                    velocity_estimation_method="unknown",
                )
            with self.assertRaisesRegex(
                ValueError,
                "band_ratio_frequency_scale_hz",
            ):
                PipelineContext(
                    work_h5=h5file,
                    holodoppler_h5=None,
                    doppler_vision_h5=None,
                    band_ratio_frequency_scale_hz=np.nan,
                )

    def test_execution_variants_are_validated_and_published(self) -> None:
        with h5py.File(
            "context_variants_test.h5",
            "w",
            driver="core",
            backing_store=False,
        ) as h5file:
            ctx = PipelineContext(
                work_h5=h5file,
                holodoppler_h5=None,
                doppler_vision_h5=None,
            )
            variant = ExecutionVariant(
                name="alternate",
                state={"result": object()},
                output_namespace="ProcessingBandRatio",
                artifact_namespace="alternate",
                provenance={"method": "alternate"},
            )

            ctx.publish_execution_variants(
                [variant],
                failures={"unavailable": "missing input"},
            )

            self.assertEqual((variant,), ctx.execution_variants)
            self.assertEqual(
                {"unavailable": "missing input"},
                ctx.execution_variant_failures,
            )
            with self.assertRaisesRegex(ValueError, "Duplicate execution variant"):
                ctx.publish_execution_variants([variant, variant])
            with self.assertRaisesRegex(ValueError, "invalid artifact namespace"):
                ctx.publish_execution_variants(
                    [
                        ExecutionVariant(
                            name="unsafe",
                            state={},
                            output_namespace="ProcessingMoments",
                            artifact_namespace="../outside",
                            provenance={},
                        )
                    ]
                )

    def test_source_array_casts_during_numeric_hdf5_read(self) -> None:
        with h5py.File(
            "context_array_test.h5",
            "w",
            driver="core",
            backing_store=False,
        ) as h5file:
            h5file.create_dataset(
                "values",
                data=np.arange(6, dtype=np.float64).reshape(2, 3),
            )
            reader = RawH5SourceReader(h5file=h5file, label="HD")

            values = reader.array("values", dtype=np.float32)

        self.assertEqual(np.float32, values.dtype)
        np.testing.assert_array_equal(values, np.arange(6, dtype=np.float32).reshape(2, 3))


if __name__ == "__main__":
    unittest.main()
