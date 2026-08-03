"""Tests for pipeline runtime context helpers."""

from __future__ import annotations

import unittest

import h5py
import numpy as np

from pipeline_engine import PipelineContext
from pipeline_engine.context import RawH5SourceReader


class PipelineContextTests(unittest.TestCase):
    def test_log_emits_to_callback(self) -> None:
        messages: list[str] = []
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
                on_log=messages.append,
            )

            ctx.log("Starting test pipeline...")

        self.assertEqual(["Starting test pipeline..."], messages)

    def test_log_without_callback_is_noop(self) -> None:
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
                    "waveform_velocity": ("per_beat", "hemifield"),
                    "waveform_shape_metrics": (),
                },
                pipeline_order=(
                    "waveform_velocity_core",
                    "waveform_velocity",
                ),
            )

            self.assertTrue(ctx.option_enabled("per_beat"))
            self.assertTrue(
                ctx.option_enabled("hemifield", pipeline="waveform_velocity")
            )
            self.assertFalse(
                ctx.option_enabled("hemifield", pipeline="waveform_shape_metrics")
            )
            self.assertEqual(
                frozenset({"per_beat", "hemifield"}),
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
            reader = RawH5SourceReader(h5file=h5file, label="HD")

            values = reader.array("values", dtype=np.float32)

        self.assertEqual(np.float32, values.dtype)
        np.testing.assert_array_equal(values, np.arange(6, dtype=np.float32).reshape(2, 3))


if __name__ == "__main__":
    unittest.main()
