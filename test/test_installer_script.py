"""Static regressions for the generated Inno Setup installer script."""

from __future__ import annotations

import re
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


class InstallerScriptTests(unittest.TestCase):
    def test_desktop_shortcut_task_is_checked_by_default(self) -> None:
        script = (REPO_ROOT / "installer/EyeFlow.iss.in").read_text(encoding="utf-8")
        task = re.search(
            r'^Name: "desktopicon";[^\r\n]+$',
            script,
            flags=re.MULTILINE,
        )

        self.assertIsNotNone(task)
        self.assertNotIn("Flags:", task.group(0))
        self.assertNotIn("Flags: checked", script)

    def test_matplotlib_postscript_backend_is_bundled(self) -> None:
        spec = (REPO_ROOT / "installer/eyeflow.spec").read_text(encoding="utf-8")

        self.assertIn(
            '"matplotlib.backends.backend_ps"',
            spec,
        )


if __name__ == "__main__":
    unittest.main()
