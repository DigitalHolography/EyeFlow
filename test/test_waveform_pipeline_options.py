"""Tests for selectable waveform pipeline product groups."""

from __future__ import annotations

import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from pipelines.waveform_shape_metrics import runner as metric_runner
from pipelines.waveform_velocity import runner as velocity_runner
from pipelines.waveform_velocity_core import runner as core_runner


class _State:
    def __init__(self, values=None):
        self.values = dict(values or {})

    def get(self, key, default=None):
        return self.values.get(key, default)

    def set(self, key, value):
        self.values[key] = value


def _context(options, state_values=None):
    return SimpleNamespace(
        state=_State(state_values),
        options_for=lambda pipeline: frozenset(options.get(pipeline, ())),
        pipeline_scheduled=lambda pipeline: pipeline in {
            "waveform_velocity_core",
            "waveform_velocity",
            "waveform_shape_metrics",
        },
        option_enabled=lambda name, pipeline=None: name
        in options.get(pipeline or "waveform_velocity", ()),
    )


class WaveformPipelineOptionTests(unittest.TestCase):
    def test_pipeline_implementation_ownership_is_cleanly_split(self) -> None:
        pipeline_root = Path(__file__).resolve().parents[1] / "src" / "pipelines"
        metrics_root = pipeline_root / "waveform_shape_metrics"
        velocity_root = pipeline_root / "waveform_velocity"
        core_root = pipeline_root / "waveform_velocity_core"

        self.assertFalse((metrics_root / "velocity").exists())
        core_source = "\n".join(
            path.read_text(encoding="utf-8") for path in core_root.rglob("*.py")
        )
        velocity_source = "\n".join(
            path.read_text(encoding="utf-8") for path in velocity_root.rglob("*.py")
        )
        self.assertNotIn("pipelines.waveform_velocity.", core_source)
        self.assertNotIn("pipelines.waveform_shape_metrics", core_source)
        self.assertNotIn("pipelines.waveform_shape_metrics", velocity_source)

    def test_velocity_parent_always_publishes_base_velocity_only(self) -> None:
        context = SimpleNamespace(dopplerview_analysis={})
        ctx = _context(
            {"waveform_velocity": ()},
            {core_runner.WAVEFORM_CONTEXT_STATE: context},
        )

        with (
            patch.object(
                velocity_runner,
                "pack_continuous_velocity_outputs",
                return_value={"base": 1},
            ),
            patch.object(
                velocity_runner,
                "run_velocity_per_beat_metrics",
            ) as per_beat,
            patch.object(
                velocity_runner,
                "pack_cross_section_profile_outputs",
            ) as profiles,
            patch.object(
                velocity_runner,
                "pack_hemifield_velocity_outputs",
            ) as hemifield,
        ):
            outputs = velocity_runner.run_waveform_velocity(ctx)

        self.assertEqual({"base": 1}, outputs)
        per_beat.assert_not_called()
        profiles.assert_not_called()
        hemifield.assert_not_called()

    def test_velocity_children_publish_their_selected_products(self) -> None:
        per_beat_result = SimpleNamespace(cycle_boundary_indexes=(0, 5, 10))
        velocity_outputs = {"per_beat": 2}
        context = SimpleNamespace(
            dopplerview_analysis={},
            artery_segment_result="artery",
            vein_segment_result="vein",
            per_beat_analysis=SimpleNamespace(cycle_boundary_indexes=(1, 6, 11)),
            source_data=SimpleNamespace(provenance={"beat_index_base": 1}),
        )
        ctx = _context(
            {
                "waveform_velocity": (
                    "velocity_profiles",
                    "per_beat",
                    "hemifield",
                )
            },
            {
                core_runner.WAVEFORM_CONTEXT_STATE: context,
                core_runner.VELOCITY_PER_BEAT_RESULT_STATE: per_beat_result,
                core_runner.VELOCITY_PER_BEAT_OUTPUTS_STATE: velocity_outputs,
            },
        )

        with (
            patch.object(
                velocity_runner,
                "pack_continuous_velocity_outputs",
                return_value={"base": 1},
            ),
            patch.object(
                velocity_runner,
                "pack_cross_section_profile_outputs",
                return_value={"profile": 3},
            ) as profiles,
            patch.object(
                velocity_runner,
                "pack_hemifield_velocity_outputs",
                return_value={"hemifield": 4},
            ) as hemifield,
        ):
            outputs = velocity_runner.run_waveform_velocity(ctx)

        self.assertEqual(
            {"base": 1, "per_beat": 2, "profile": 3, "hemifield": 4},
            outputs,
        )
        profiles.assert_called_once_with(
            "artery",
            "vein",
            (0, 5, 10),
            index_base=0,
        )
        hemifield.assert_called_once_with(
            velocity_outputs,
            context.source_data,
            "artery",
            "vein",
        )

    def test_no_metric_children_skips_per_beat_metric_work(self) -> None:
        ctx = _context({"waveform_shape_metrics": ()})

        self.assertEqual({}, metric_runner.run_waveform_shape_metrics(ctx))
        self.assertFalse(core_runner._per_beat_required(ctx))


if __name__ == "__main__":
    unittest.main()
