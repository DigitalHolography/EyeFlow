"""Tests for selectable velocity-analysis product groups."""

from __future__ import annotations

import unittest
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pytest

from pipeline_engine.planning.base import PIPELINE_REGISTRY
from pipelines import load_pipeline_catalog
from pipelines.velocity_analysis import runner as analysis_builder


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


@pytest.mark.parametrize("scheduled", [False, True])
@pytest.mark.parametrize(
    ("options", "expected"),
    [
        ((), False),
        (("segments", "quadrants"), False),
        (("segment_velocity_maps",), False),
        (("velocity_profiles",), True),
        (("velocity_profile_analysis",), True),
        (("velocity_profile_fft",), True),
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
        ((), False),
        (("masked_edges",), False),
        (("gradient_edges",), True),
        (("gradient_edges", "masked_edges"), True),
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
        descriptor = PIPELINE_REGISTRY["velocity_analysis"]
        for selected, present, absent in (
            ((), (), ("segments", "maps", "profiles", "profile_fits", "profile_fft", "quadrants")),
            (("segments",), ("segments",), ("maps", "profiles", "profile_fft")),
            (("velocity_profiles",), ("segments", "profiles"), ("maps", "profile_fits", "profile_fft")),
            (("velocity_profile_analysis",), ("segments", "profiles", "profile_fits"), ("maps", "profile_fft")),
            (("velocity_profile_fft",), ("segments", "profiles", "profile_fft"), ("maps", "profile_fits")),
            (("segment_velocity_maps",), ("segments", "maps", "map_videos"), ("profiles", "profile_fft")),
        ):
            stages = descriptor.build_graph(selected).plan().order
            self.assertTrue(set(present).issubset(stages))
            self.assertTrue(set(absent).isdisjoint(stages))


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


    def test_core_skips_segment_extraction_when_not_required(self) -> None:
        analysis = SimpleNamespace(
            cycle_boundary_indexes=np.asarray([0, 1], dtype=np.int32),
            cardiac_cycle=SimpleNamespace(spectral="cardiac_cycle"),
            continuous=lambda vessel, raw=False: np.asarray(
                [1.0, 2.0], dtype=np.float32
            ),
        )
        source = SimpleNamespace(
            timing=SimpleNamespace(dt_seconds=0.1),
            provenance={"beat_index_base": 0},
        )
        ctx = SimpleNamespace()

        with patch.object(analysis_builder, "_segment_velocity_inputs") as extract:
            artery, vein = analysis_builder._build_segments(
                ctx,
                analysis,
                source,
                segments_required=False,
            )

        extract.assert_not_called()
        self.assertIsNone(artery)
        self.assertIsNone(vein)

    def test_core_plain_profiles_do_not_stream_fft_or_retain_velocity_maps(self) -> None:
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
            provenance={"beat_index_base": 1},
        )
        ctx = SimpleNamespace(
            pipeline_scheduled=lambda name: name == "velocity_analysis",
            option_enabled=lambda name, pipeline=None: name == "velocity_profiles",
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
            artery, vein = analysis_builder._segment_velocity_inputs(
                "velocity_map",
                source,
                "rings",
                ctx,
                cycle_boundary_indexes=(1, 6, 11),
            )

        self.assertEqual(("artery", "vein"), (artery, vein))
        self.assertFalse(analyze.call_args.kwargs["retain_velocity_maps"])
        self.assertFalse(analyze.call_args.kwargs["velocity_profile_fft"])
        self.assertEqual(
            (1, 6, 11),
            analyze.call_args.kwargs["cycle_boundary_indexes"],
        )
        self.assertEqual(1, analyze.call_args.kwargs["index_base"])

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


    def test_pdf_report_requires_pulse_pngs(self) -> None:
        ctx = _context(
            {"velocity_analysis": (), "waveform_shape_metrics": ()},
            scheduled={
                "velocity",
                "velocity_analysis",
                "waveform_shape_metrics",
                "pdf_report",
            },
        )

        self.assertTrue(analysis_builder._pulse_pngs_required(ctx))

    def test_downstream_dependency_does_not_require_pulse_pngs(self) -> None:
        ctx = _context(
            {"velocity_analysis": ()},
            scheduled={"velocity", "velocity_analysis", "blood_volume_rate"},
            targeted={"blood_volume_rate"},
        )

        self.assertFalse(analysis_builder._pulse_pngs_required(ctx))


if __name__ == "__main__":
    unittest.main()
