"""Keep shipped defaults synchronized with the visible pipeline catalog."""

from __future__ import annotations

import json
import unittest
from pathlib import Path

from app_settings import normalize_pipeline_visibility
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
        self.assertTrue(settings["pipeline_visibility"]["blood_volume_rate"])
        self.assertEqual(
            {"gradient_edges": True, "masked_edges": True},
            settings["pipeline_options"]["blood_volume_rate"],
        )

    def test_new_default_selected_pipeline_is_enabled_in_existing_settings(self) -> None:
        visibility, changed = normalize_pipeline_visibility(
            ("waveform_velocity", "blood_volume_rate"),
            {"waveform_velocity": False},
            missing_defaults={"blood_volume_rate": True},
        )

        self.assertTrue(changed)
        self.assertEqual(
            {"waveform_velocity": False, "blood_volume_rate": True},
            visibility,
        )


if __name__ == "__main__":
    unittest.main()
