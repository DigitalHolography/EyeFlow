"""Tests for run-log file ownership."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from app_settings import LAST_RUN_LOG_FILENAME
from utils.logger import (
    Logger,
    format_log_message,
)


class LoggerTests(unittest.TestCase):
    def tearDown(self) -> None:
        Logger.reset_current()

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

    def test_legacy_level_markers_are_normalized(self) -> None:
        self.assertEqual("[ERROR] Failure details", format_log_message("[FAIL] Failure details"))
        self.assertEqual(
            "[WARNING] Warning details",
            format_log_message("[WARNING] Warning details"),
        )
        self.assertEqual("[INPUT] source.h5", format_log_message("[INPUT] source.h5"))


if __name__ == "__main__":
    unittest.main()
