"""Static regressions for the generated Inno Setup installer script."""

from __future__ import annotations

import re
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


class InstallerScriptTests(unittest.TestCase):
    def test_desktop_shortcut_task_is_checked_by_default(self) -> None:
        script = (REPO_ROOT / "build_installer.ps1").read_text(encoding="utf-8")
        task = re.search(
            r'^Name: "desktopicon";[^\r\n]+$',
            script,
            flags=re.MULTILINE,
        )

        self.assertIsNotNone(task)
        self.assertNotIn("Flags:", task.group(0))
        self.assertNotIn("Flags: checked", script)


if __name__ == "__main__":
    unittest.main()
