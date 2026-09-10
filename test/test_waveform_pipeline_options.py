"""Tests for selectable waveform pipeline product groups."""

from __future__ import annotations

import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from input_output.schema import EyeFlowOutputPaths
from pipeline_engine.base import PIPELINE_REGISTRY
from pipelines import load_pipeline_catalog
from pipelines.lowrank_waveform_decomposition import runner as lowrank_runner
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


def _context(options, state_values=None, scheduled=None):
    scheduled = set(
        scheduled
        or {
            "waveform_velocity_core",
            "waveform_velocity",
            "waveform_shape_metrics",
        }
    )
    return SimpleNamespace(
        state=_State(state_values),
        options_for=lambda pipeline: frozenset(options.get(pipeline, ())),
        pipeline_scheduled=lambda pipeline: pipeline in scheduled,
        option_enabled=lambda name, pipeline=None: name
        in options.get(pipeline or "waveform_velocity", ()),
    )


class WaveformPipelineOptionTests(unittest.TestCase):
    def test_fft_option_is_disabled_by_default_and_requires_profiles(self) -> None:
        load_pipeline_catalog()
        options = {
            option.name: option
            for option in PIPELINE_REGISTRY["waveform_velocity"].options
        }
        fft = options["velocity_profile_fft"]
        self.assertFalse(fft.default_enabled)
        self.assertEqual(("velocity_profiles",), fft.requires)

    def test_lowrank_pipeline_includes_veins_and_selected_quadrants(self) -> None:
        velocity_outputs = {"per_beat": 1}
        context = SimpleNamespace(
            source_data="source",
            artery_segment_result="artery",
            vein_segment_result="vein",
        )
        ctx = SimpleNamespace(
            state=_State(
                {
                    core_runner.VELOCITY_PER_BEAT_OUTPUTS_STATE: velocity_outputs,
                    core_runner.WAVEFORM_CONTEXT_STATE: context,
                }
            ),
            options_for=lambda _pipeline: frozenset({"quadrants"}),
        )

        with patch.object(
            lowrank_runner,
            "pack_lowrank_waveform_decomposition_outputs",
            return_value={"lowrank": 2},
        ) as pack:
            outputs = lowrank_runner.run_lowrank_waveform_decomposition(ctx)

        self.assertEqual({"lowrank": 2}, outputs)
        pack.assert_called_once_with(
            velocity_outputs,
            vein_flag=True,
            include_quadrants=True,
            source_data="source",
            artery_segments="artery",
            vein_segments="vein",
        )

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
        context = SimpleNamespace(velocity_analysis={})
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
                "pack_quadrant_velocity_outputs",
            ) as quadrants,
        ):
            outputs = velocity_runner.run_waveform_velocity(ctx)

        self.assertEqual({"base": 1}, outputs)
        per_beat.assert_not_called()
        profiles.assert_not_called()
        quadrants.assert_not_called()

    def test_velocity_children_publish_their_selected_products(self) -> None:
        per_beat_result = SimpleNamespace(cycle_boundary_indexes=(0, 5, 10))
        schema = EyeFlowOutputPaths.active()
        velocity_outputs = {
            "per_beat": 2,
            schema.artery_per_beat.segment_velocity_signal: 5,
        }
        context = SimpleNamespace(
            velocity_analysis={},
            artery_segment_result="artery",
            vein_segment_result="vein",
            per_beat_analysis=SimpleNamespace(cycle_boundary_indexes=(1, 6, 11)),
            source_data=SimpleNamespace(provenance={"beat_index_base": 1}),
        )
        ctx = _context(
            {
                "waveform_velocity": (
                    "velocity_profiles",
                    "velocity_profile_fft",
                    "segment_velocity_maps",
                    "per_beat",
                    "quadrants",
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
                "prepare_segment_velocity_maps_per_beat",
                return_value=("artery_maps", "vein_maps"),
            ) as prepare_maps,
            patch.object(
                velocity_runner,
                "pack_segment_map_outputs",
                return_value={"maps": 8},
            ) as maps,
            patch.object(
                velocity_runner,
                "pack_displacement_segment_map_outputs",
                return_value={"displacement_maps": 9},
            ) as displacement_maps,
            patch.object(
                velocity_runner,
                "pack_velocity_profile_fft_outputs",
                return_value={"fft_profile": 7},
            ) as fft_profiles,
            patch.object(
                velocity_runner,
                "pack_quadrant_velocity_outputs",
                return_value={"quadrants": 4},
            ) as quadrants,
        ):
            outputs = velocity_runner.run_waveform_velocity(ctx)

        self.assertEqual(
            {
                "base": 1,
                "per_beat": 2,
                "profile": 3,
                "fft_profile": 7,
                "maps": 8,
                "displacement_maps": 9,
                "quadrants": 4,
            },
            outputs,
        )
        profiles.assert_called_once_with(
            "artery",
            "vein",
            (0, 5, 10),
            index_base=0,
        )
        fft_profiles.assert_called_once_with(
            "artery",
            "vein",
        )
        prepare_maps.assert_called_once_with(
            "artery",
            "vein",
            (1, 6, 11),
            index_base=1,
        )
        maps.assert_called_once_with(
            "artery",
            "vein",
            "artery_maps",
            "vein_maps",
        )
        displacement_maps.assert_called_once_with(
            "artery",
            "vein",
            (1, 6, 11),
            index_base=1,
        )
        quadrants.assert_called_once_with(
            velocity_outputs,
            context.source_data,
            "artery",
            "vein",
        )

    def test_segments_option_does_not_build_velocity_maps(self) -> None:
        context = SimpleNamespace(
            velocity_analysis={},
            artery_segment_result="artery",
            vein_segment_result="vein",
            per_beat_analysis=SimpleNamespace(cycle_boundary_indexes=(1, 6, 11)),
            source_data=SimpleNamespace(provenance={"beat_index_base": 1}),
        )
        ctx = _context(
            {"waveform_velocity": ("segments",)},
            {core_runner.WAVEFORM_CONTEXT_STATE: context},
        )
        ctx.output = SimpleNamespace(available=True)

        with (
            patch.object(
                velocity_runner,
                "pack_continuous_velocity_outputs",
                return_value={"base": 1},
            ),
            patch.object(
                velocity_runner,
                "pack_segment_velocity_outputs",
                return_value={"signals": 2},
            ),
            patch.object(
                velocity_runner,
                "pack_segment_map_outputs",
                return_value={"maps": 3},
            ) as maps,
            patch.object(
                velocity_runner,
                "export_segment_velocity_map_avis",
                return_value=["artery.avi", "vein.avi"],
            ) as avis,
        ):
            outputs = velocity_runner.run_waveform_velocity(ctx)

        self.assertEqual({"base": 1, "signals": 2}, outputs)
        maps.assert_not_called()
        avis.assert_not_called()

    def test_velocity_profiles_do_not_build_per_beat_velocity_maps(self) -> None:
        context = SimpleNamespace(
            velocity_analysis={},
            artery_segment_result="artery",
            vein_segment_result="vein",
            per_beat_analysis=SimpleNamespace(cycle_boundary_indexes=(1, 6, 11)),
            source_data=SimpleNamespace(provenance={"beat_index_base": 1}),
        )
        ctx = _context(
            {"waveform_velocity": ("velocity_profiles",)},
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
                "prepare_segment_velocity_maps_per_beat",
            ) as prepare_maps,
            patch.object(
                velocity_runner,
                "pack_cross_section_profile_outputs",
                return_value={"profiles": 2},
            ),
            patch.object(
                velocity_runner,
                "pack_velocity_profile_fft_outputs",
                return_value={"fft": 3},
            ) as fft_profiles,
        ):
            outputs = velocity_runner.run_waveform_velocity(ctx)

        self.assertEqual({"base": 1, "profiles": 2}, outputs)
        prepare_maps.assert_not_called()
        fft_profiles.assert_not_called()

    def test_segment_velocity_maps_option_publishes_maps_and_avis(self) -> None:
        context = SimpleNamespace(
            velocity_analysis={},
            artery_segment_result="artery",
            vein_segment_result="vein",
            per_beat_analysis=SimpleNamespace(cycle_boundary_indexes=(1, 6, 11)),
            source_data=SimpleNamespace(provenance={"beat_index_base": 1}),
        )
        ctx = _context(
            {"waveform_velocity": ("segment_velocity_maps",)},
            {core_runner.WAVEFORM_CONTEXT_STATE: context},
        )
        ctx.output = SimpleNamespace(available=True)

        with (
            patch.object(
                velocity_runner,
                "pack_continuous_velocity_outputs",
                return_value={"base": 1},
            ),
            patch.object(
                velocity_runner,
                "pack_segment_velocity_outputs",
            ) as segment_outputs,
            patch.object(
                velocity_runner,
                "pack_segment_map_outputs",
                return_value={"maps": 3},
            ) as maps,
            patch.object(
                velocity_runner,
                "prepare_segment_velocity_maps_per_beat",
                return_value=("artery_maps", "vein_maps"),
            ) as prepare_maps,
            patch.object(
                velocity_runner,
                "export_segment_velocity_map_avis",
                return_value=["artery.avi", "vein.avi"],
            ) as avis,
        ):
            outputs = velocity_runner.run_waveform_velocity(ctx)

        self.assertEqual({"base": 1, "maps": 3}, outputs)
        segment_outputs.assert_not_called()
        maps.assert_called_once_with(
            "artery",
            "vein",
            "artery_maps",
            "vein_maps",
        )
        prepare_maps.assert_called_once_with(
            "artery",
            "vein",
            (1, 6, 11),
            index_base=1,
        )
        avis.assert_called_once_with(
            ctx.output,
            "artery",
            "vein",
            {"maps": 3},
        )

    def test_no_metric_children_skips_per_beat_metric_work(self) -> None:
        ctx = _context({"waveform_shape_metrics": ()})

        self.assertEqual({}, metric_runner.run_waveform_shape_metrics(ctx))
        self.assertFalse(core_runner._per_beat_required(ctx))

    def test_core_segment_requirement_uses_synchronized_segment_selection(self) -> None:
        ctx = _context(
            {
                "waveform_velocity": ("per_beat", "segments"),
                "waveform_shape_metrics": ("per_beat", "segments"),
            }
        )

        self.assertTrue(core_runner._segments_required(ctx))

        ctx = _context(
            {
                "waveform_velocity": ("segment_velocity_maps",),
                "waveform_shape_metrics": (),
            }
        )

        self.assertTrue(core_runner._segments_required(ctx))

        ctx = _context(
            {
                "waveform_velocity": (),
                "waveform_shape_metrics": ("per_beat", "segments"),
            }
        )

        self.assertFalse(core_runner._segments_required(ctx))
        self.assertFalse(core_runner._per_beat_required(ctx))

        ctx = _context(
            {
                "waveform_velocity": (),
                "waveform_shape_metrics": ("per_beat", "quadrants"),
            }
        )

        self.assertTrue(core_runner._segments_required(ctx))
        self.assertTrue(core_runner._per_beat_required(ctx))

        ctx = _context(
            {
                "waveform_velocity": ("per_beat",),
                "waveform_shape_metrics": ("per_beat",),
            }
        )

        self.assertFalse(core_runner._segments_required(ctx))

    def test_global_shape_metrics_can_run_without_core_segments(self) -> None:
        context = SimpleNamespace(
            source_data="source",
            artery_segment_result=None,
            vein_segment_result=None,
        )
        ctx = _context(
            {
                "waveform_shape_metrics": ("per_beat",),
            },
            {
                core_runner.WAVEFORM_CONTEXT_STATE: context,
                core_runner.VELOCITY_PER_BEAT_OUTPUTS_STATE: {"global": 1},
            },
        )

        with patch.object(
            metric_runner,
            "pack_waveform_shape_outputs",
            return_value={"shape": 1},
        ) as pack:
            outputs = metric_runner.run_waveform_shape_metrics(ctx)

        self.assertEqual({"shape": 1}, outputs)
        self.assertFalse(pack.call_args.kwargs["include_segments"])

    def test_core_skips_segment_extraction_when_not_required(self) -> None:
        analysis = {
            "retinal_artery_velocity_signal": [1.0, 2.0],
            "retinal_vein_velocity_signal": [1.0, 2.0],
            "beat_indices": [0, 1],
        }
        source = SimpleNamespace(
            timing=SimpleNamespace(dt_seconds=0.1),
            provenance={"beat_index_base": 0},
            optic_disc_width=None,
            optic_disc_height=None,
        )
        ctx = SimpleNamespace()

        with (
            patch.object(core_runner, "_segment_velocity_inputs") as extract,
            patch.object(
                core_runner,
                "spectral_heartbeat_analysis",
                return_value="heartbeat",
            ),
        ):
            _, artery, vein = core_runner._per_beat_input_from_analysis(
                analysis,
                source,
                source.timing,
                4,
                ctx,
                segments_required=False,
            )

        extract.assert_not_called()
        self.assertIsNone(artery)
        self.assertIsNone(vein)

    def test_core_plain_profiles_do_not_stream_fft_or_retain_velocity_maps(self) -> None:
        source = SimpleNamespace(
            retinal_artery_mask="artery_mask",
            retinal_vein_mask="vein_mask",
            optic_disc_center=(10, 10),
            optic_disc_mask="disc_mask",
            cross_section_settings="settings",
            provenance={"beat_index_base": 1},
        )
        ctx = SimpleNamespace(
            pipeline_scheduled=lambda name: name == "waveform_velocity",
            option_enabled=lambda name, pipeline=None: name == "velocity_profiles",
            inputs=SimpleNamespace(
                hd=SimpleNamespace(filename="hd.h5"),
                dv=SimpleNamespace(filename="dv.h5"),
            ),
            state=SimpleNamespace(raw={}),
            output=SimpleNamespace(available=False),
        )

        with patch.object(
            core_runner,
            "analyze_velocity_segments",
            return_value={"artery": "artery", "vein": "vein"},
        ) as analyze:
            artery, vein = core_runner._segment_velocity_inputs(
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
            retinal_artery_mask="artery_mask",
            retinal_vein_mask="vein_mask",
            optic_disc_center=(10, 10),
            optic_disc_mask="disc_mask",
            cross_section_settings="settings",
            provenance={"beat_index_base": 0},
        )
        ctx = SimpleNamespace(
            pipeline_scheduled=lambda name: name == "waveform_velocity",
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
            core_runner,
            "analyze_velocity_segments",
            return_value={"artery": "artery", "vein": "vein"},
        ) as analyze:
            core_runner._segment_velocity_inputs(
                "velocity_map",
                source,
                "rings",
                ctx,
                cycle_boundary_indexes=(0, 5, 10),
            )
        self.assertTrue(analyze.call_args.kwargs["velocity_profile_fft"])
        self.assertFalse(analyze.call_args.kwargs["retain_velocity_maps"])

    def test_analysis_schedule_generates_both_source_profiles_automatically(self) -> None:
        context = SimpleNamespace(
            velocity_analysis={},
            artery_segment_result="artery",
            vein_segment_result="vein",
            per_beat_analysis=SimpleNamespace(cycle_boundary_indexes=(0, 5, 10)),
            source_data=SimpleNamespace(provenance={"beat_index_base": 0}),
        )
        ctx = _context(
            {"waveform_velocity": ()},
            {core_runner.WAVEFORM_CONTEXT_STATE: context},
            scheduled={"waveform_velocity", "velocity_profile_analysis"},
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
                return_value={"both_profiles": 2},
            ) as profiles,
            patch.object(
                velocity_runner,
                "pack_velocity_profile_fft_outputs",
            ) as fft,
        ):
            outputs = velocity_runner.run_waveform_velocity(ctx)
        self.assertEqual({"base": 1, "both_profiles": 2}, outputs)
        profiles.assert_called_once_with("artery", "vein", (0, 5, 10), index_base=0)
        fft.assert_not_called()
        self.assertTrue(core_runner._segments_required(ctx))

    def test_pdf_report_requires_shared_per_beat_products(self) -> None:
        ctx = _context(
            {"waveform_velocity": (), "waveform_shape_metrics": ()},
            scheduled={
                "waveform_velocity_core",
                "waveform_velocity",
                "waveform_shape_metrics",
                "pdf_report",
            },
        )

        self.assertTrue(core_runner._per_beat_required(ctx))
        self.assertTrue(core_runner._pulse_pngs_required(ctx))

    def test_pdf_report_publishes_velocity_per_beat_outputs(self) -> None:
        context = SimpleNamespace(velocity_analysis={})
        result = SimpleNamespace(cycle_boundary_indexes=(0, 2))
        ctx = _context(
            {"waveform_velocity": ()},
            {
                core_runner.WAVEFORM_CONTEXT_STATE: context,
                core_runner.VELOCITY_PER_BEAT_RESULT_STATE: result,
                core_runner.VELOCITY_PER_BEAT_OUTPUTS_STATE: {"per_beat": 1},
            },
            scheduled={"waveform_velocity", "pdf_report"},
        )

        with patch.object(
            velocity_runner,
            "pack_continuous_velocity_outputs",
            return_value={"base": 1},
        ):
            outputs = velocity_runner.run_waveform_velocity(ctx)

        self.assertEqual({"base": 1, "per_beat": 1}, outputs)

if __name__ == "__main__":
    unittest.main()
