"""Keep shipped defaults synchronized with the visible pipeline catalog."""

from __future__ import annotations

import json
import unittest
from pathlib import Path

from app_settings import normalize_pipeline_options
from pipeline_engine import PipelineDAG
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
            {"gradient_edges": False, "masked_edges": True},
            settings["pipeline_options"]["blood_volume_rate"],
        )
        self.assertFalse(
            settings["pipeline_options"]["velocity_analysis"][
                "velocity_profiles"
            ]
        )
        self.assertFalse(
            settings["pipeline_options"]["velocity_analysis"][
                "velocity_profile_analysis"
            ]
        )
        self.assertEqual(
            {"segment_maps": False},
            settings["pipeline_options"]["displacement_profiles"],
        )
        self.assertEqual("both", settings["velocity_estimation_method"])
        self.assertEqual(1.0, settings["band_ratio_frequency_scale_hz"])
        self.assertEqual(12, settings["band_limited_harmonics"])
        self.assertFalse(settings["copy_hd_passthrough"])

    def test_renamed_pipeline_options_are_migrated(self) -> None:
        options, changed = normalize_pipeline_options(
            {"velocity_analysis": ("segments", "quadrants")},
            {
                "waveform_velocity": {
                    "segments": False,
                    "quadrants": True,
                }
            },
        )

        self.assertTrue(changed)
        self.assertEqual(
            {"velocity_analysis": {"segments": False, "quadrants": True}},
            options,
        )

    def test_release_defaults_exclude_gradient_and_velocity_profile_outputs(
        self,
    ) -> None:
        settings = json.loads(
            (ROOT / "default_settings.json").read_text(encoding="utf-8")
        )
        available, missing = load_pipeline_catalog()
        targets = tuple(
            name
            for name, enabled in settings["pipeline_visibility"].items()
            if enabled
        )
        options = {
            name: tuple(
                option
                for option, enabled in values.items()
                if enabled
            )
            for name, values in settings["pipeline_options"].items()
        }

        plan = PipelineDAG((*available, *missing)).resolve_targets(
            targets,
            pipeline_options=options,
        )

        self.assertNotIn("spatial_gradient_moment0", plan.names)
        self.assertNotIn(
            "gradient_edges",
            options["blood_volume_rate"],
        )
        self.assertNotIn(
            "velocity_profiles",
            options["velocity_analysis"],
        )
        self.assertNotIn(
            "velocity_profile_analysis",
            options["velocity_analysis"],
        )


if __name__ == "__main__":
    unittest.main()
