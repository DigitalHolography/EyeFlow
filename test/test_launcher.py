"""Tests for GUI/CLI launcher dispatch."""

from __future__ import annotations

import runpy
import sys
import unittest
from contextlib import redirect_stderr, redirect_stdout
from io import StringIO
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

import launcher  # noqa: E402

FROZEN_ENTRY = SRC_DIR.parent / "installer/entry_point.py"


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


class FrozenLauncherTests(unittest.TestCase):
    def _run(self, arguments):
        with (
            patch.object(sys, "argv", ["EyeFlow.exe", *arguments]),
            patch("multiprocessing.freeze_support") as freeze_support,
            self.assertRaises(SystemExit) as exited,
        ):
            runpy.run_path(str(FROZEN_ENTRY), run_name="__main__")
        freeze_support.assert_called_once_with()
        return exited.exception.code

    def test_frozen_cli_help_prints_and_exits_successfully(self):
        stdout = StringIO()
        with redirect_stdout(stdout):
            status = self._run(["--help"])
        self.assertEqual(0, status)
        self.assertIn("--data", stdout.getvalue())

    def test_frozen_cli_argument_errors_keep_status_two(self):
        stderr = StringIO()
        with redirect_stderr(stderr):
            status = self._run(["--data"])
        self.assertEqual(2, status)
        self.assertIn("expected one argument", stderr.getvalue())

    def test_frozen_cli_runtime_failure_prints_and_keeps_status_one(self):
        stderr = StringIO()
        with (
            patch("cli.run_cli", side_effect=ValueError("invalid input")),
            redirect_stderr(stderr),
        ):
            status = self._run(["--data", "missing.holo"])
        self.assertEqual(1, status)
        self.assertIn("Error: invalid input", stderr.getvalue())

    def test_frozen_cli_success_keeps_status_zero(self):
        with patch("cli.run_cli", return_value=0):
            status = self._run(["--data", "scan.holo"])
        self.assertEqual(0, status)

    def test_frozen_gui_without_arguments_exits_normally(self):
        with patch("launcher._call_entry", return_value=None) as call_entry:
            status = self._run([])
        self.assertIsNone(status)
        call_entry.assert_called_once_with("eye_flow", "main")

    def test_packaging_check_preserves_its_exit_code(self):
        check_main = unittest.mock.Mock(return_value=7)
        with patch.dict(sys.modules, {"frozen_check": SimpleNamespace(main=check_main)}):
            status = self._run(["--packaging-check", "report.json"])
        self.assertEqual(7, status)
        check_main.assert_called_once_with(["report.json"])


if __name__ == "__main__":
    unittest.main()
