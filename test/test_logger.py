"""Tests for run-log file ownership."""

from __future__ import annotations

import tempfile
import unittest
import multiprocessing
import os
import sys
from pathlib import Path
from unittest.mock import patch

from app_settings import LAST_RUN_LOG_FILENAME
from utils.logger import (
    Logger,
    format_log_message,
)


def _spawn_log_probe(output_queue):
    Logger.reset_current()
    messages = []
    Logger.configure(on_log=messages.append)
    Logger.log_debug("worker debug")
    Logger.log_perf("worker performance")
    Logger.log_info("worker operational")
    output_queue.put(messages)


class LoggerTests(unittest.TestCase):
    def setUp(self) -> None:
        self._logging_env = {
            key: os.environ.pop(key, None)
            for key in ("EYEFLOW_DIAGNOSTICS", "EYEFLOW_PERF_LOGGING", "EYEFLOW_LOG_LEVEL")
        }

    def tearDown(self) -> None:
        Logger.reset_current()
        for key, value in self._logging_env.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value

    def test_writes_versioned_settings_sibling_log_file(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            settings_path = Path(temp_dir) / "settings.json"
            logger = Logger.configure(settings_path)

            log_path = logger.write_snapshot("hello")

            self.assertEqual(settings_path.with_name(LAST_RUN_LOG_FILENAME), log_path)
            self.assertEqual(log_path, logger.last_saved_path)
            self.assertEqual("hello", log_path.read_text(encoding="utf-8"))

    def test_ensure_file_creates_empty_log(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            logger = Logger.configure(Path(temp_dir) / "settings.json")

            log_path = logger.ensure_file()

            self.assertTrue(log_path.exists())
            self.assertEqual("", log_path.read_text(encoding="utf-8"))

    def test_current_requires_configuration(self) -> None:
        with self.assertRaises(RuntimeError):
            Logger.current()

    def test_configure_sets_current_logger(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            settings_path = Path(temp_dir) / "settings.json"

            logger = Logger.configure(settings_path)

            self.assertIs(logger, Logger.current())
            self.assertEqual(settings_path.with_name(LAST_RUN_LOG_FILENAME), logger.path)

    def test_reconfigure_updates_the_existing_logger(self) -> None:
        first = Logger.configure()
        second = Logger.configure(on_log=lambda _message: None)

        self.assertIs(first, second)
        self.assertIs(first, Logger.current())

    def test_log_level_is_selected_without_embedded_markers(self) -> None:
        messages: list[str] = []
        Logger.configure(on_log=messages.append)

        Logger.log_error("Something went wrong")
        Logger.log_warning("The vessel crop is small")
        Logger.log_debug("Timing detail")

        self.assertEqual(
            [
                "[ERROR] Something went wrong",
                "[WARNING] The vessel crop is small",
                "[DEBUG] Timing detail",
            ],
            messages,
        )

    def test_development_mode_shows_debug_and_performance(self) -> None:
        messages: list[str] = []
        with patch.dict(os.environ, {
            "EYEFLOW_DIAGNOSTICS": "1", "EYEFLOW_PERF_LOGGING": "1",
            "EYEFLOW_LOG_LEVEL": "DEBUG",
        }):
            Logger.configure(on_log=messages.append)
            Logger.log_debug("debug detail")
            Logger.log_perf("timing detail")
            Logger.log_info("normal operation")

        self.assertEqual(
            ["[DEBUG] debug detail", "[PERF] timing detail", "[INFO] normal operation"],
            messages,
        )

    def test_packaged_mode_hides_diagnostics_but_keeps_operational_levels(self) -> None:
        messages: list[str] = []
        with (
            patch.object(sys, "frozen", True, create=True),
            patch.dict(os.environ, {
                "EYEFLOW_DIAGNOSTICS": "0", "EYEFLOW_PERF_LOGGING": "0",
                "EYEFLOW_LOG_LEVEL": "INFO",
            }),
        ):
            Logger.configure(on_log=messages.append)
            Logger.log_perf("hidden performance")
            Logger.log("[PERF] legacy performance")
            Logger.log_debug("hidden detail")
            Logger.log_info("started")
            Logger.log_warning("warning")
            Logger.log_error("failure")

        self.assertEqual(
            ["[INFO] started", "[WARNING] warning", "[ERROR] failure"],
            messages,
        )

    def test_configured_severity_threshold(self) -> None:
        messages: list[str] = []
        with patch.dict(os.environ, {
            "EYEFLOW_DIAGNOSTICS": "0", "EYEFLOW_PERF_LOGGING": "0",
            "EYEFLOW_LOG_LEVEL": "INFO",
        }):
            Logger.configure(on_log=messages.append, level="WARNING")
            Logger.log("[INPUT] source.holo")
            Logger.log_info("below threshold")
            Logger.log_warning("at threshold")
            Logger.log_error("above threshold")

        self.assertEqual(["[WARNING] at threshold", "[ERROR] above threshold"], messages)

    def test_diagnostic_override_is_inherited_by_spawned_workers(self) -> None:
        context = multiprocessing.get_context("spawn")
        for enabled, level, perf, expected in (
            ("0", "INFO", "0", ["[INFO] worker operational"]),
            (
                "1",
                "DEBUG",
                "1",
                [
                    "[DEBUG] worker debug",
                    "[PERF] worker performance",
                    "[INFO] worker operational",
                ],
            ),
        ):
            with patch.dict(
                os.environ,
                {
                    "EYEFLOW_DIAGNOSTICS": enabled,
                    "EYEFLOW_PERF_LOGGING": perf,
                    "EYEFLOW_LOG_LEVEL": level,
                },
            ):
                output = context.Queue()
                worker = context.Process(target=_spawn_log_probe, args=(output,))
                worker.start()
                worker.join(timeout=20)
                if worker.is_alive():
                    worker.terminate()
                    worker.join()
                self.assertEqual(0, worker.exitcode)
                self.assertEqual(expected, output.get(timeout=2))
                output.close()
                output.join_thread()

    def test_explicit_diagnostic_override_enables_production_logging(self) -> None:
        messages: list[str] = []
        with (
            patch.object(sys, "frozen", True, create=True),
            patch.dict(os.environ, {}, clear=False),
        ):
            os.environ.pop("EYEFLOW_DIAGNOSTICS", None)
            os.environ.pop("EYEFLOW_PERF_LOGGING", None)
            os.environ.pop("EYEFLOW_LOG_LEVEL", None)
            Logger.configure(on_log=messages.append, diagnostics=True)
            Logger.log_debug("diagnostic detail")
            Logger.log_perf("diagnostic timing")

        self.assertEqual(
            ["[DEBUG] diagnostic detail", "[PERF] diagnostic timing"],
            messages,
        )

    def test_legacy_level_markers_are_normalized(self) -> None:
        self.assertEqual("[ERROR] Failure details", format_log_message("[FAIL] Failure details"))
        self.assertEqual("[PERF] Timing details", format_log_message("[PERF] Timing details"))
        self.assertEqual(
            "[WARNING] Warning details",
            format_log_message("[WARNING] Warning details"),
        )
        self.assertEqual("[INPUT] source.h5", format_log_message("[INPUT] source.h5"))

    def test_performance_category_is_independently_filterable(self) -> None:
        messages: list[str] = []
        with patch.dict(os.environ, {
            "EYEFLOW_DIAGNOSTICS": "0", "EYEFLOW_PERF_LOGGING": "1",
            "EYEFLOW_LOG_LEVEL": "INFO",
        }):
            Logger.configure(on_log=messages.append)
            Logger.log_debug("hidden detail")
            Logger.log_perf("visible timing")

        self.assertEqual(["[PERF] visible timing"], messages)


if __name__ == "__main__":
    unittest.main()
