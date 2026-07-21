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
