import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

import launcher  # noqa: E402


class LauncherTests(unittest.TestCase):
    def test_call_entry_imports_module_and_calls_function(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            module_dir = Path(tmp_dir)
            (module_dir / "fake_entry.py").write_text(
                "def main(value=None):\n"
                "    return ('imported', value)\n",
                encoding="utf-8",
            )
            sys.path.insert(0, str(module_dir))
            sys.modules.pop("fake_entry", None)
            try:
                result = launcher._call_entry("fake_entry", "main", "value")
            finally:
                sys.path.remove(str(module_dir))
                sys.modules.pop("fake_entry", None)

            self.assertEqual(("imported", "value"), result)

    def test_main_configures_threads_and_calls_eye_flow_entry(self) -> None:
        module = types.SimpleNamespace(main=mock.Mock(return_value="ok"))
        with (
            mock.patch.dict(sys.modules, {"eye_flow": module}),
            mock.patch("launcher.configure_numeric_threads") as configure,
        ):
            result = launcher.main()

        self.assertEqual("ok", result)
        configure.assert_called_once_with()
        module.main.assert_called_once_with()

    def test_cli_main_configures_threads_and_passes_argv(self) -> None:
        module = types.SimpleNamespace(main=mock.Mock(return_value=0))
        with (
            mock.patch.dict(sys.modules, {"cli": module}),
            mock.patch("launcher.configure_numeric_threads") as configure,
        ):
            result = launcher.cli_main(["--help"])

        self.assertEqual(0, result)
        configure.assert_called_once_with()
        module.main.assert_called_once_with(["--help"])


if __name__ == "__main__":
    unittest.main()
