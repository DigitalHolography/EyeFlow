"""Keep shipped defaults synchronized with the visible pipeline catalog."""

from __future__ import annotations

import json
import unittest
from pathlib import Path

from pipelines import load_pipeline_catalog


ROOT = Path(__file__).resolve().parents[1]


class DefaultSettingsTests(unittest.TestCase):
    def test_defaults_cover_every_visible_pipeline_and_option(self) -> None:
        settings = json.loads(
            (ROOT / "default_settings.json").read_text(encoding="utf-8")
        )
        available, missing = load_pipeline_catalog()
        visible = {
            descriptor.name: descriptor
            for descriptor in (*available, *missing)
            if descriptor.visibility != "hidden"
        }

        self.assertEqual(set(visible), set(settings["pipeline_visibility"]))

        expected_options = {
            name: {option.name for option in descriptor.options}
            for name, descriptor in visible.items()
            if descriptor.options
        }
        configured_options = {
            name: set(options)
            for name, options in settings["pipeline_options"].items()
        }
        self.assertEqual(expected_options, configured_options)


if __name__ == "__main__":
    unittest.main()
