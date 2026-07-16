"""Tests for GUI/CLI launcher dispatch."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path
from unittest.mock import patch

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

import launcher  # noqa: E402


class LauncherTests(unittest.TestCase):
    def test_eyeflow_without_arguments_launches_gui(self) -> None:
        with patch("launcher._call_entry", return_value="gui") as call_entry:
            result = launcher.main([])

        self.assertEqual("gui", result)
        call_entry.assert_called_once_with("eye_flow", "main")

    def test_eyeflow_with_any_argument_dispatches_cli(self) -> None:
        argv = ["--data", "scan.holo"]
        with patch("launcher._call_entry", return_value="cli") as call_entry:
            result = launcher.main(argv)

        self.assertEqual("cli", result)
        call_entry.assert_called_once_with("cli", "main", argv)

    def test_dedicated_cli_with_empty_arguments_never_launches_gui(self) -> None:
        with patch("launcher._call_entry", return_value="cli") as call_entry:
            result = launcher.cli_main([])

        self.assertEqual("cli", result)
        call_entry.assert_called_once_with("cli", "main", [])


if __name__ == "__main__":
    unittest.main()
