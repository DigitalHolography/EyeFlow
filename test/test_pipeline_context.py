"""Tests for pipeline runtime context helpers."""

from __future__ import annotations

import unittest

import h5py
import numpy as np

from pipeline_engine import PipelineContext
from pipeline_engine.context import RawH5SourceReader
from input_output.h5_access import PipelineH5Output, PipelineInputSource
from input_output.schema.base import TypedSource
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
                pipeline_name="waveform_velocity",
                pipeline_options={
                    "waveform_velocity": ("segments", "quadrants"),
                    "waveform_shape_metrics": (),
                },
                pipeline_order=(
                    "retinal_velocity",
                    "waveform_velocity",
                ),
                pipeline_targets=("waveform_velocity",),
            )

            self.assertTrue(ctx.option_enabled("segments"))
            self.assertTrue(
                ctx.option_enabled("quadrants", pipeline="waveform_velocity")
            )
            self.assertFalse(
                ctx.option_enabled("quadrants", pipeline="waveform_shape_metrics")
            )
            self.assertEqual(
                frozenset({"segments", "quadrants"}),
                ctx.options_for("waveform_velocity"),
            )
            self.assertTrue(ctx.pipeline_scheduled("retinal_velocity"))
            self.assertFalse(ctx.pipeline_scheduled("pdf_report"))
            self.assertTrue(ctx.pipeline_targeted("waveform_velocity"))
            self.assertFalse(ctx.pipeline_targeted("retinal_velocity"))

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
