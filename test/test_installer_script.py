"""Static regressions for the generated Inno Setup installer script."""

from __future__ import annotations

import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


class InstallerScriptTests(unittest.TestCase):
    def test_matplotlib_postscript_backend_is_bundled(self) -> None:
        spec = (REPO_ROOT / "installer/eyeflow.spec").read_text(encoding="utf-8")

        self.assertIn(
            '"matplotlib.backends.backend_ps"',
            spec,
        )


if __name__ == "__main__":
    unittest.main()
