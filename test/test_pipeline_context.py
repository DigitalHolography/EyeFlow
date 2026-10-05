"""Tests for pipeline runtime context helpers."""

from __future__ import annotations

import unittest

import h5py
import numpy as np

from pipeline_engine import PipelineContext, ProcessResult
from pipeline_engine.context import PipelineH5Output, apply_pipeline_result
from input_output.h5_access import H5Output, PipelineInputSource
from input_output.schema.base import TypedSource
from utils.logger import Logger


class PipelineContextTests(unittest.TestCase):
    def test_missing_arrays_return_none_for_explicit_none_default(self) -> None:
        with h5py.File("context_missing_test.h5", "w", driver="core", backing_store=False) as h5file:
            source = PipelineInputSource(h5file=h5file, label="HD")
            output = H5Output(h5file)
            for reader in (source, output):
                self.assertIsNone(reader.array("missing", dtype=np.float32, default=None))
                with self.assertRaises(KeyError):
                    reader.array("missing")

            typed = TypedSource(source)
            with self.assertRaises(KeyError):
                typed._array("missing")
            self.assertIsNone(typed._array("missing", default=None))

    def test_pipeline_attribute_policy_belongs_to_engine(self) -> None:
        with h5py.File(
            "context_attrs_test.h5", "w", driver="core", backing_store=False
        ) as h5file:
            generic_output = H5Output(h5file)
            generic_output.set_attr("pipeline", "generic")
            self.assertEqual("generic", h5file.attrs["pipeline"])
            del h5file.attrs["pipeline"]

            ctx = PipelineContext(
                work_h5=h5file,
                holodoppler_h5=None,
                doppler_vision_h5=None,
            )
            self.assertIsInstance(ctx.output.h5, PipelineH5Output)
            ctx.output.h5.set_attr("pipeline", "direct")
            ctx.output.h5.set_attrs({"pipeline": "direct-batch", "source": "direct"})
            apply_pipeline_result(
                ctx,
                ProcessResult(metrics={}, attrs={"pipeline": "result", "unit": "pixel"}),
            )

            self.assertNotIn("pipeline", h5file.attrs)
            self.assertEqual("direct", h5file.attrs["source"])
            self.assertEqual("pixel", h5file.attrs["unit"])

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
                    "waveform_velocity": ("per_beat", "quadrants"),
                    "waveform_shape_metrics": (),
                },
                pipeline_order=(
                    "waveform_velocity_core",
                    "waveform_velocity",
                ),
            )

            self.assertTrue(ctx.option_enabled("per_beat"))
            self.assertTrue(
                ctx.option_enabled("quadrants", pipeline="waveform_velocity")
            )
            self.assertFalse(
                ctx.option_enabled("quadrants", pipeline="waveform_shape_metrics")
            )
            self.assertEqual(
                frozenset({"per_beat", "quadrants"}),
                ctx.options_for("waveform_velocity"),
            )
            self.assertTrue(ctx.pipeline_scheduled("waveform_velocity_core"))
            self.assertFalse(ctx.pipeline_scheduled("pdf_report"))

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
            reader = PipelineInputSource(h5file=h5file, label="HD")

            values = reader.array("values", dtype=np.float32)

        self.assertEqual(np.float32, values.dtype)
        np.testing.assert_array_equal(values, np.arange(6, dtype=np.float32).reshape(2, 3))


if __name__ == "__main__":
    unittest.main()
