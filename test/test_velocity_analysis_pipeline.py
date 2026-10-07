"""Tests for selectable velocity-analysis product groups."""

from __future__ import annotations

import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

from input_output.schema import EyeFlowOutputPaths
from pipeline_engine.base import PIPELINE_REGISTRY
from pipelines import load_pipeline_catalog
from pipelines.lowrank_waveform_decomposition import runner as lowrank_runner
from pipelines.velocity_analysis import builder as analysis_builder
from pipelines.velocity_analysis import runner as analysis_runner
from pipelines.waveform_shape_metrics import runner as metric_runner


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


class VelocityAnalysisPipelineTests(unittest.TestCase):
    def test_profile_options_are_disabled_by_default(self) -> None:
        load_pipeline_catalog()
        options = {
            option.name: option
            for option in PIPELINE_REGISTRY["velocity_analysis"].options
        }
        profiles = options["velocity_profiles"]
        self.assertNotIn("per_beat", options)
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

    def test_lowrank_pipeline_includes_veins_and_selected_quadrants(self) -> None:
        velocity_outputs = {"per_beat": 1}
        context = SimpleNamespace(
            velocity={},
            source_data="source",
            artery_segments="artery",
            vein_segments="vein",
            per_beat_result="per-beat",
        )
        ctx = SimpleNamespace(
            state=_State(
                {
                    analysis_builder.VELOCITY_ANALYSIS_STATE: context,
                }
            ),
            options_for=lambda _pipeline: frozenset({"quadrants"}),
        )

        with (
            patch.object(
                lowrank_runner,
                "velocity_analysis",
                return_value=context,
            ),
            patch.object(
                lowrank_runner,
                "pack_velocity_per_beat_inputs",
                return_value=velocity_outputs,
            ),
            patch.object(
                lowrank_runner,
                "pack_lowrank_waveform_decomposition_outputs",
                return_value={"lowrank": 2},
            ) as pack,
        ):
            outputs = lowrank_runner.run_lowrank_waveform_decomposition(ctx)

        self.assertEqual({"lowrank": 2}, outputs)
        pack.assert_called_once_with(
            velocity_outputs,
            vein_flag=True,
            include_quadrants=True,
            artery_segments="artery",
            vein_segments="vein",
        )

    def test_pipeline_implementation_ownership_is_cleanly_split(self) -> None:
        pipeline_root = Path(__file__).resolve().parents[1] / "src" / "pipelines"
        metrics_root = pipeline_root / "waveform_shape_metrics"
        velocity_root = pipeline_root / "velocity_analysis"
        gradient_root = pipeline_root / "spatial_gradient_moment0"

        self.assertFalse((metrics_root / "velocity").exists())
        self.assertFalse((pipeline_root / "velocity_analysis_core").exists())
        self.assertTrue((velocity_root / "builder.py").is_file())
        for package in ("analysis", "artifacts", "outputs"):
            self.assertTrue((velocity_root / package / "__init__.py").is_file())
        for obsolete in (
            "workflow.py",
            "per_beat.py",
            "per_beat_outputs.py",
            "segment_maps.py",
            "segment_velocity_map_avi.py",
        ):
            self.assertFalse((velocity_root / obsolete).exists())
        velocity_source = "\n".join(
            path.read_text(encoding="utf-8") for path in velocity_root.rglob("*.py")
        )
        self.assertNotIn("pipelines.waveform_shape_metrics", velocity_source)
        self.assertNotIn("spatial_gradient", velocity_source)
        self.assertNotIn("displacement", (velocity_root / "runner.py").read_text())
        self.assertTrue((gradient_root / "profiles.py").is_file())

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

    def test_velocity_figures_export_whenever_safe_per_beat_data_is_available(
        self,
    ) -> None:
        schema = EyeFlowOutputPaths.active()
        artery = (np.ones((3, 2, 1, 1), dtype=np.float32), {"unit": "mm/s"})
        vein = (np.full((3, 2, 1, 1), 2.0, dtype=np.float32), {"unit": "mm/s"})
        velocity_outputs = {
            schema.artery_per_beat_safe.velocity_signal: artery,
            schema.vein_per_beat_safe.velocity_signal: vein,
        }
        context = SimpleNamespace(
            velocity={},
            artery_segments=None,
            vein_segments=None,
            per_beat_result="per-beat",
        )
        ctx = _context(
            {"velocity_analysis": ()},
            {
                analysis_builder.VELOCITY_ANALYSIS_STATE: context,
            },
        )
        ctx.output = SimpleNamespace(available=True)

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
            patch.object(analysis_runner, "export_velocity_signals") as export,
        ):
            outputs = analysis_runner.run_velocity_analysis(ctx)

        self.assertEqual({"base": 1, **velocity_outputs}, outputs)
        export.assert_called_once_with(ctx.output, artery, vein)

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

    def test_segments_option_does_not_build_velocity_maps(self) -> None:
        context = SimpleNamespace(
            velocity={},
            artery_segments="artery",
            vein_segments="vein",
            cycle_boundary_indexes=(0, 5, 10),
            per_beat_result="per-beat",
            source_data=SimpleNamespace(provenance={"beat_index_base": 1}),
        )
        ctx = _context(
            {"velocity_analysis": ("segments",)},
            {analysis_builder.VELOCITY_ANALYSIS_STATE: context},
        )
        ctx.output = SimpleNamespace(available=True)

        with (
            patch.object(
                analysis_runner,
                "pack_continuous_velocity_outputs",
                return_value={"base": 1},
            ),
            patch.object(
                analysis_runner,
                "pack_segment_velocity_outputs",
                return_value={"signals": 2},
            ),
            patch.object(
                analysis_runner,
                "pack_velocity_per_beat_outputs",
                return_value={"per_beat": 4},
            ),
            patch.object(
                analysis_runner,
                "pack_segment_map_outputs",
                return_value={"maps": 3},
            ) as maps,
            patch.object(
                analysis_runner,
                "export_segment_velocity_map_avis",
                return_value=["artery.avi", "vein.avi"],
            ) as avis,
        ):
            outputs = analysis_runner.run_velocity_analysis(ctx)

        self.assertEqual(
            {"base": 1, "signals": 2, "per_beat": 4},
            outputs,
        )
        maps.assert_not_called()
        avis.assert_not_called()

    def test_velocity_profiles_do_not_build_per_beat_velocity_maps(self) -> None:
        context = SimpleNamespace(
            velocity={},
            artery_segments="artery",
            vein_segments="vein",
            cycle_boundary_indexes=(0, 5, 10),
            per_beat_result="per-beat",
            source_data=SimpleNamespace(provenance={"beat_index_base": 1}),
        )
        ctx = _context(
            {"velocity_analysis": ("velocity_profiles",)},
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
                "prepare_segment_velocity_maps_per_beat",
            ) as prepare_maps,
            patch.object(
                analysis_runner,
                "pack_segment_velocity_outputs",
                return_value={"signals": 4},
            ),
            patch.object(
                analysis_runner,
                "pack_velocity_per_beat_outputs",
                return_value={"per_beat": 5},
            ),
            patch.object(
                analysis_runner,
                "pack_cross_section_profile_outputs",
                return_value={"profiles": 2},
            ),
            patch.object(
                analysis_runner,
                "pack_velocity_profile_fft_outputs",
                return_value={"fft": 3},
            ) as fft_profiles,
        ):
            outputs = analysis_runner.run_velocity_analysis(ctx)

        self.assertEqual(
            {"base": 1, "signals": 4, "per_beat": 5, "profiles": 2},
            outputs,
        )
        prepare_maps.assert_not_called()
        fft_profiles.assert_not_called()

    def test_segment_velocity_maps_option_publishes_maps_and_avis(self) -> None:
        context = SimpleNamespace(
            velocity={},
            artery_segments="artery",
            vein_segments="vein",
            cycle_boundary_indexes=(0, 5, 10),
            per_beat_result="per-beat",
            source_data=SimpleNamespace(provenance={"beat_index_base": 1}),
        )
        ctx = _context(
            {"velocity_analysis": ("segment_velocity_maps",)},
            {analysis_builder.VELOCITY_ANALYSIS_STATE: context},
        )
        ctx.output = SimpleNamespace(available=True)

        with (
            patch.object(
                analysis_runner,
                "pack_continuous_velocity_outputs",
                return_value={"base": 1},
            ),
            patch.object(
                analysis_runner,
                "pack_segment_velocity_outputs",
                return_value={"signals": 2},
            ) as segment_outputs,
            patch.object(
                analysis_runner,
                "pack_velocity_per_beat_outputs",
                return_value={"per_beat": 4},
            ),
            patch.object(
                analysis_runner,
                "pack_segment_map_outputs",
                return_value={"maps": 3},
            ) as maps,
            patch.object(
                analysis_runner,
                "prepare_segment_velocity_maps_per_beat",
                return_value=("artery_maps", "vein_maps"),
            ) as prepare_maps,
            patch.object(
                analysis_runner,
                "export_segment_velocity_map_avis",
                return_value=["artery.avi", "vein.avi"],
            ) as avis,
        ):
            outputs = analysis_runner.run_velocity_analysis(ctx)

        self.assertEqual(
            {"base": 1, "signals": 2, "per_beat": 4, "maps": 3},
            outputs,
        )
        segment_outputs.assert_called_once()
        maps.assert_called_once_with(
            "artery",
            "vein",
            "artery_maps",
            "vein_maps",
            velocity={},
        )
        prepare_maps.assert_called_once_with(
            "artery",
            "vein",
            (0, 5, 10),
            index_base=0,
        )
        avis.assert_called_once_with(
            ctx.output,
            "artery",
            "vein",
            {"maps": 3},
        )

    def test_shape_metrics_include_global_per_beat_outputs_by_default(self) -> None:
        context = SimpleNamespace(
            velocity={},
            source_data="source",
            artery_segments=None,
            vein_segments=None,
            per_beat_result="per-beat",
        )
        ctx = _context({"waveform_shape_metrics": ()})

        with (
            patch.object(
                metric_runner,
                "velocity_analysis",
                return_value=context,
            ),
            patch.object(
                metric_runner,
                "pack_velocity_per_beat_inputs",
                return_value={"global": 1},
            ),
            patch.object(
                metric_runner,
                "pack_waveform_shape_outputs",
                return_value={"shape": 1},
            ) as pack,
        ):
            outputs = metric_runner.run_waveform_shape_metrics(ctx)

        self.assertEqual({"shape": 1}, outputs)
        self.assertTrue(pack.call_args.kwargs["include_per_beat"])
        self.assertFalse(pack.call_args.kwargs["include_segments"])
        self.assertFalse(pack.call_args.kwargs["include_quadrants"])

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

    def test_global_shape_metrics_can_run_without_core_segments(self) -> None:
        context = SimpleNamespace(
            velocity={},
            source_data="source",
            artery_segments=None,
            vein_segments=None,
            per_beat_result="per-beat",
        )
        ctx = _context(
            {
                "waveform_shape_metrics": (),
            },
            {
                analysis_builder.VELOCITY_ANALYSIS_STATE: context,
            },
        )

        with (
            patch.object(
                metric_runner,
                "velocity_analysis",
                return_value=context,
            ),
            patch.object(
                metric_runner,
                "pack_velocity_per_beat_inputs",
                return_value={"global": 1},
            ),
            patch.object(
                metric_runner,
                "pack_waveform_shape_outputs",
                return_value={"shape": 1},
            ) as pack,
        ):
            outputs = metric_runner.run_waveform_shape_metrics(ctx)

        self.assertEqual({"shape": 1}, outputs)
        self.assertFalse(pack.call_args.kwargs["include_segments"])

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
                "run_velocity_profile_analysis",
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

    def test_pdf_report_publishes_velocity_per_beat_outputs(self) -> None:
        result = SimpleNamespace(cycle_boundary_indexes=(0, 2))
        context = SimpleNamespace(
            velocity={},
            artery_segments=None,
            vein_segments=None,
            per_beat_result=result,
        )
        ctx = _context(
            {"velocity_analysis": ()},
            {
                analysis_builder.VELOCITY_ANALYSIS_STATE: context,
            },
            scheduled={"velocity_analysis", "pdf_report"},
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
                return_value={"per_beat": 1},
            ),
        ):
            outputs = analysis_runner.run_velocity_analysis(ctx)

        self.assertEqual({"base": 1, "per_beat": 1}, outputs)

if __name__ == "__main__":
    unittest.main()
