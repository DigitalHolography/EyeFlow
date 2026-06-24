"""Tests for pipeline runtime context helpers."""

from __future__ import annotations

import unittest

import h5py

from pipeline_engine import PipelineContext


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


if __name__ == "__main__":
    unittest.main()
