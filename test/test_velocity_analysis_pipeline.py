"""Tests for selectable velocity-analysis product groups."""

from __future__ import annotations

import unittest
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pytest

from input_output.schema import EyeFlowOutputPaths
from pipeline_engine.base import PIPELINE_REGISTRY
from pipelines import load_pipeline_catalog
from pipelines.velocity_analysis import builder as analysis_builder
from pipelines.velocity_analysis import runner as analysis_runner


class _State:
    def __init__(self, values=None):
        self.values = dict(values or {})

    def get(self, key, default=None):
        return self.values.get(key, default)

    def set(self, key, value):
        self.values[key] = value


def _context(options, state_values=None, scheduled=None, targeted=None):
    scheduled = set(
        scheduled
        or {
            "velocity",
            "velocity_analysis",
            "waveform_shape_metrics",
        }
    )
    return SimpleNamespace(
        state=_State(state_values),
        options_for=lambda pipeline: frozenset(options.get(pipeline, ())),
        pipeline_scheduled=lambda pipeline: pipeline in scheduled,
        pipeline_targeted=lambda pipeline: pipeline in (
            scheduled if targeted is None else set(targeted)
        ),
        option_enabled=lambda name, pipeline=None: name
        in options.get(pipeline or "velocity_analysis", ()),
    )


@pytest.mark.parametrize(
    ("options", "scheduled", "expected"),
    [
        (("velocity_profiles",), True, True),
    ],
)
def test_longitudinal_velocity_profiles_require_a_profile_consumer(
    options, expected, scheduled,
):
    source = SimpleNamespace(
        source=SimpleNamespace(
            segmentation=SimpleNamespace(
                vessels=SimpleNamespace(artery="artery_mask", vein="vein_mask"),
                optic_disc="optic_disc",
            ),
        ),
        profile_settings="settings",
        provenance={"beat_index_base": 0},
    )
    ctx = _context(
        {"velocity_analysis": options},
        scheduled={"velocity_analysis"} if scheduled else {"waveform_shape_metrics"},
    )
    ctx.output = SimpleNamespace(available=False)
    with patch.object(
        analysis_builder,
        "analyze_velocity_segment_profiles",
        return_value={"artery": "artery", "vein": "vein"},
    ) as analyze:
        analysis_builder._segment_velocity_inputs(
            "velocity_map", source, "rings", ctx, cycle_boundary_indexes=(0, 5, 10),
        )
    assert analyze.call_args.kwargs["compute_longitudinal_profiles"] is (
        expected and scheduled
    )
    assert analyze.call_args.kwargs["compute_transverse_profiles"] is (
        expected and scheduled
    )


@pytest.mark.parametrize(
    ("blood_volume_rate_options", "expected"),
    [
        (("gradient_edges",), True),
    ],
)
def test_only_gradient_edge_flow_requires_transverse_velocity_profiles(
    blood_volume_rate_options, expected,
):
    source = SimpleNamespace(
        source=SimpleNamespace(
            segmentation=SimpleNamespace(
                vessels=SimpleNamespace(artery="artery_mask", vein="vein_mask"),
                optic_disc="optic_disc",
            ),
        ),
        profile_settings="settings",
        provenance={"beat_index_base": 0},
    )
    ctx = _context(
        {
            "velocity_analysis": ("segments",),
            "blood_volume_rate": blood_volume_rate_options,
        },
        scheduled={"velocity_analysis", "blood_volume_rate"},
    )
    ctx.output = SimpleNamespace(available=False)
    with patch.object(
        analysis_builder,
        "analyze_velocity_segment_profiles",
        return_value={"artery": "artery", "vein": "vein"},
    ) as analyze:
        analysis_builder._segment_velocity_inputs(
            "velocity_map", source, "rings", ctx, cycle_boundary_indexes=(0, 5, 10),
        )
    assert analyze.call_args.kwargs["compute_transverse_profiles"] is expected
    assert analyze.call_args.kwargs["compute_longitudinal_profiles"] is False


class VelocityAnalysisPipelineTests(unittest.TestCase):
    def test_profile_options_are_disabled_by_default(self) -> None:
        load_pipeline_catalog()
        options = {
            option.name: option
            for option in PIPELINE_REGISTRY["velocity_analysis"].options
        }
        profiles = options["velocity_profiles"]
        self.assertEqual(("segments",), options["segment_velocity_maps"].requires)
        self.assertEqual(("segments",), options["velocity_profiles"].requires)
        self.assertEqual(("segments",), options["quadrants"].requires)
        analysis = options["velocity_profile_analysis"]
        fft = options["velocity_profile_fft"]
        self.assertFalse(profiles.default_enabled)
        self.assertFalse(analysis.default_enabled)
        self.assertFalse(fft.default_enabled)
        self.assertEqual(("velocity_profiles",), analysis.requires)
        self.assertEqual(("velocity_profiles",), fft.requires)

    def test_velocity_analysis_always_publishes_base_and_per_beat_velocity(self) -> None:
        context = SimpleNamespace(
            velocity={},
            artery_segments=None,
            vein_segments=None,
            per_beat_result="per-beat",
        )
        ctx = _context(
            {"velocity_analysis": ()},
            {analysis_builder.VELOCITY_ANALYSIS_STATE: context},
        )

        with (
            patch.object(
                analysis_runner,
                "pack_continuous_velocity_outputs",
                return_value={"base": 1},
            ),
            patch.object(
                analysis_runner,
                "pack_cross_section_profile_outputs",
            ) as profiles,
            patch.object(
                analysis_runner,
                "pack_velocity_per_beat_outputs",
                return_value={"per_beat": 2},
            ) as per_beat,
            patch.object(
                analysis_runner,
                "pack_quadrant_velocity_outputs",
            ) as quadrants,
        ):
            outputs = analysis_runner.run_velocity_analysis(ctx)

        self.assertEqual({"base": 1, "per_beat": 2}, outputs)
        per_beat.assert_called_once_with("per-beat", velocity={})
        profiles.assert_not_called()
        quadrants.assert_not_called()

    def test_velocity_children_publish_their_selected_products(self) -> None:
        per_beat_result = SimpleNamespace(cycle_boundary_indexes=(0, 5, 10))
        schema = EyeFlowOutputPaths.active()
        artery_segments = SimpleNamespace(
            topology=SimpleNamespace(optic_disc_center_xy=(12.0, 13.0))
        )
        velocity_outputs = {
            "per_beat": 2,
            schema.artery_per_beat.segment_velocity_signal: 5,
        }
        context = SimpleNamespace(
            velocity={},
            artery_segments=artery_segments,
            vein_segments="vein",
            cycle_boundary_indexes=(0, 5, 10),
            per_beat_result=per_beat_result,
            source_data=SimpleNamespace(
                provenance={"beat_index_base": 1},
                profile_settings=SimpleNamespace(pixel_size_mm=0.01),
            ),
        )
        ctx = _context(
            {
                "velocity_analysis": (
                    "segments",
                    "velocity_profiles",
                    "velocity_profile_fft",
                    "segment_velocity_maps",
                    "quadrants",
                )
            },
            {
                analysis_builder.VELOCITY_ANALYSIS_STATE: context,
            },
        )

        with (
            patch.object(
                analysis_runner,
                "pack_continuous_velocity_outputs",
                return_value={"base": 1},
            ),
            patch.object(
                analysis_runner,
                "pack_velocity_per_beat_outputs",
                return_value=velocity_outputs,
            ),
            patch.object(
                analysis_runner,
                "pack_segment_velocity_outputs",
                return_value={"segment_signals": 6},
            ),
            patch.object(
                analysis_runner,
                "pack_cross_section_profile_outputs",
                return_value={"profile": 3},
            ) as profiles,
            patch.object(
                analysis_runner,
                "prepare_segment_velocity_maps_per_beat",
                return_value=("artery_maps", "vein_maps"),
            ) as prepare_maps,
            patch.object(
                analysis_runner,
                "pack_segment_map_outputs",
                return_value={"maps": 8},
            ) as maps,
            patch.object(
                analysis_runner,
                "pack_velocity_profile_fft_outputs",
                return_value={"fft_profile": 7},
            ) as fft_profiles,
            patch.object(
                analysis_runner,
                "pack_quadrant_velocity_outputs",
                return_value={"quadrants": 4},
            ) as quadrants,
        ):
            outputs = analysis_runner.run_velocity_analysis(ctx)

        self.assertEqual(
            {
                "base": 1,
                **velocity_outputs,
                "segment_signals": 6,
                "profile": 3,
                "fft_profile": 7,
                "maps": 8,
                "quadrants": 4,
            },
            outputs,
        )
        profiles.assert_called_once_with(
            artery_segments,
            "vein",
            (0, 5, 10),
            index_base=0,
            velocity={},
        )
        fft_profiles.assert_called_once_with(
            artery_segments,
            "vein",
        )
        prepare_maps.assert_called_once_with(
            artery_segments,
            "vein",
            (0, 5, 10),
            index_base=0,
        )
        maps.assert_called_once_with(
            artery_segments,
            "vein",
            "artery_maps",
            "vein_maps",
            velocity={},
        )
        quadrants.assert_called_once_with(
            velocity_outputs,
            context.source_data,
            artery_segments,
            "vein",
            velocity={},
        )

    def test_core_segment_requirement_uses_synchronized_segment_selection(self) -> None:
        ctx = _context(
            {
                "velocity_analysis": ("segments",),
                "waveform_shape_metrics": ("segments",),
            }
        )

        self.assertTrue(analysis_builder._segments_required(ctx))

        ctx = _context(
            {
                "velocity_analysis": ("segment_velocity_maps",),
                "waveform_shape_metrics": (),
            }
        )

        self.assertTrue(analysis_builder._segments_required(ctx))

        ctx = _context(
            {
                "velocity_analysis": (),
                "waveform_shape_metrics": ("segments",),
            }
        )

        self.assertTrue(analysis_builder._segments_required(ctx))

        ctx = _context(
            {
                "velocity_analysis": (),
                "waveform_shape_metrics": ("quadrants",),
            }
        )

        self.assertTrue(analysis_builder._segments_required(ctx))

        ctx = _context(
            {
                "velocity_analysis": (),
                "waveform_shape_metrics": (),
            }
        )

        self.assertFalse(analysis_builder._segments_required(ctx))

        ctx = _context(
            {
                "velocity_analysis": (),
                "waveform_shape_metrics": ("segments",),
            },
            scheduled={"velocity", "velocity_analysis"},
        )

        self.assertFalse(analysis_builder._segments_required(ctx))

    def test_core_explicit_fft_option_streams_fft_without_retaining_maps(self) -> None:
        source = SimpleNamespace(
            source=SimpleNamespace(
                segmentation=SimpleNamespace(
                    vessels=SimpleNamespace(
                        artery="artery_mask",
                        vein="vein_mask",
                    ),
                    optic_disc="optic_disc",
                ),
            ),
            profile_settings="settings",
            provenance={"beat_index_base": 0},
        )
        ctx = SimpleNamespace(
            pipeline_scheduled=lambda name: name == "velocity_analysis",
            option_enabled=lambda name, pipeline=None: name
            in {"velocity_profiles", "velocity_profile_fft"},
            inputs=SimpleNamespace(
                hd=SimpleNamespace(filename="hd.h5"),
                dv=SimpleNamespace(filename="dv.h5"),
            ),
            state=SimpleNamespace(raw={}),
            output=SimpleNamespace(available=False),
        )
        with patch.object(
            analysis_builder,
            "analyze_velocity_segment_profiles",
            return_value={"artery": "artery", "vein": "vein"},
        ) as analyze:
            analysis_builder._segment_velocity_inputs(
                "velocity_map",
                source,
                "rings",
                ctx,
                cycle_boundary_indexes=(0, 5, 10),
            )
        self.assertTrue(analyze.call_args.kwargs["velocity_profile_fft"])
        self.assertFalse(analyze.call_args.kwargs["retain_velocity_maps"])

    def test_analysis_option_generates_and_analyzes_both_source_profiles(self) -> None:
        context = SimpleNamespace(
            velocity={},
            artery_segments="artery",
            vein_segments="vein",
            cycle_boundary_indexes=(0, 5, 10),
            per_beat_result="per-beat",
            source_data=SimpleNamespace(provenance={"beat_index_base": 0}),
        )
        ctx = _context(
            {
                "velocity_analysis": (
                    "segments",
                    "velocity_profiles",
                    "velocity_profile_analysis",
                )
            },
            {analysis_builder.VELOCITY_ANALYSIS_STATE: context},
        )
        with (
            patch.object(
                analysis_runner,
                "pack_continuous_velocity_outputs",
                return_value={"base": 1},
            ),
            patch.object(
                analysis_runner,
                "pack_cross_section_profile_outputs",
                return_value={"both_profiles": 2},
            ) as profiles,
            patch.object(
                analysis_runner,
                "pack_segment_velocity_outputs",
                return_value={"signals": 3},
            ),
            patch.object(
                analysis_runner,
                "pack_velocity_per_beat_outputs",
                return_value={"per_beat": 4},
            ),
            patch.object(
                analysis_runner,
                "pack_velocity_profile_fft_outputs",
            ) as fft,
            patch.object(
                analysis_runner,
                "pack_velocity_profile_analysis_outputs",
                return_value={"analysis": 5},
            ) as analyze,
        ):
            outputs = analysis_runner.run_velocity_analysis(ctx)
        self.assertEqual(
            {
                "base": 1,
                "signals": 3,
                "per_beat": 4,
                "both_profiles": 2,
                "analysis": 5,
            },
            outputs,
        )
        profiles.assert_called_once_with(
            "artery",
            "vein",
            (0, 5, 10),
            index_base=0,
            velocity={},
        )
        analyze.assert_called_once_with({"both_profiles": 2})
        fft.assert_not_called()
        self.assertTrue(analysis_builder._segments_required(ctx))

if __name__ == "__main__":
    unittest.main()
