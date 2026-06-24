"""Tests for run-log file ownership."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from app_settings import LAST_RUN_LOG_FILENAME
from utils.logger import RunLogger, configure_logger, current_logger


class RunLoggerTests(unittest.TestCase):
    def tearDown(self) -> None:
        RunLogger.reset_current()

    def test_writes_versioned_settings_sibling_log_file(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            settings_path = Path(temp_dir) / "settings.json"
            logger = RunLogger(settings_path)

            log_path = logger.write_snapshot("hello")

            self.assertEqual(settings_path.with_name(LAST_RUN_LOG_FILENAME), log_path)
            self.assertEqual(log_path, logger.last_saved_path)
            self.assertEqual("hello", log_path.read_text(encoding="utf-8"))

    def test_ensure_file_creates_empty_log(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            logger = RunLogger(Path(temp_dir) / "settings.json")

            log_path = logger.ensure_file()

            self.assertTrue(log_path.exists())
            self.assertEqual("", log_path.read_text(encoding="utf-8"))

    def test_current_requires_configuration(self) -> None:
        with self.assertRaises(RuntimeError):
            RunLogger.current()

    def test_configure_sets_current_logger(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            settings_path = Path(temp_dir) / "settings.json"

            logger = configure_logger(settings_path)

            self.assertIs(logger, current_logger())
            self.assertEqual(settings_path.with_name(LAST_RUN_LOG_FILENAME), logger.path)


if __name__ == "__main__":
    unittest.main()
