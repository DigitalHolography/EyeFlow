"""Tests for dependency-aware pipeline-library selection."""

from __future__ import annotations

import unittest
import tempfile
from pathlib import Path

from app_settings import UserSettings
from types import SimpleNamespace
from unittest.mock import Mock

from pipeline_engine import PipelineDAG, PipelineDescriptor, PipelineOption
from ui.controllers.pipeline_library import (
    PipelineLibraryController,
    option_selection_names,
    pipeline_ui_sort_key,
    status_text_wraplength,
)


def _descriptor(
    name: str,
    *,
    requires=(),
    produces=(),
    options=(),
    visibility="visible",
) -> PipelineDescriptor:
    return PipelineDescriptor(
        name=name,
        description=name,
        available=True,
        dag_requires=tuple(requires),
        dag_produces=tuple(produces),
        options=tuple(options),
        visibility=visibility,
    )


class PipelineLibraryDependencyTests(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.settings_path = Path(directory.name) / "settings.json"

    def make_app(self, **values):
        settings = UserSettings(path=self.settings_path)
        settings.pipeline_visibility.replace(values.pop("pipeline_visibility", {}))
        settings.pipeline_options.replace(values.pop("pipeline_options", {}))
        return SimpleNamespace(settings=settings, settings_controller=Mock(), **values)

    def test_status_text_is_limited_to_half_the_library_width(self) -> None:
        self.assertEqual(
            399,
            status_text_wraplength(
                900,
                horizontal_padding=48,
                divider_width=6,
                fallback=320,
            ),
        )
        self.assertEqual(
            320,
            status_text_wraplength(
                1,
                horizontal_padding=48,
                divider_width=6,
                fallback=320,
            ),
        )

    def test_pipeline_ui_order_is_explicit_and_extensible(self) -> None:
        rows = [
            _descriptor("pdf_report"),
            _descriptor("another_pipeline"),
            _descriptor("waveform_shape_metrics"),
            _descriptor("velocity_analysis"),
        ]

        self.assertEqual(
            [
                "velocity_analysis",
                "waveform_shape_metrics",
                "pdf_report",
                "another_pipeline",
            ],
            [item.name for item in sorted(rows, key=pipeline_ui_sort_key)],
        )

    def test_waveform_status_reports_both_velocity_estimators(self) -> None:
        estimator_var = SimpleNamespace(get=lambda: "frequency_bands")
        controller = PipelineLibraryController(
            SimpleNamespace(velocity_estimation_method_var=estimator_var)
        )

        status = controller._pipeline_status_text(
            _descriptor("velocity_analysis")
        )

        self.assertEqual("Available \u2014 estimators: Moments and Band ratio", status)

    def test_dag_exposes_transitive_upstream_and_downstream_relations(self) -> None:
        core = _descriptor(
            "core",
            produces=("core_output",),
            visibility="hidden",
        )
        velocity = _descriptor("velocity", requires=("core_output",))
        metrics = _descriptor("metrics", requires=("core_output",))
        report = _descriptor("report", requires=("velocity", "metrics"))
        dag = PipelineDAG((core, velocity, metrics, report))

        self.assertEqual(
            ("core", "velocity", "metrics"),
            dag.dependencies_of("report", transitive=True),
        )
        self.assertEqual(
            ("velocity", "metrics", "report"),
            dag.dependents_of("core", transitive=True),
        )

    def test_selected_targets_stay_explicit_while_upstream_is_derived(self) -> None:
        core = _descriptor("core", visibility="hidden")
        velocity = _descriptor("velocity", requires=("core",))
        metrics = _descriptor("metrics", requires=("core",))
        report = _descriptor("report", requires=("velocity", "metrics"))
        catalog = {
            item.name: item for item in (core, velocity, metrics, report)
        }
        app = self.make_app(
            pipeline_catalog=catalog,
            pipeline_dag=PipelineDAG(catalog.values()),
            pipeline_visibility={
                "velocity": False,
                "metrics": False,
                "report": False,
            },
            pipeline_visibility_vars={},
            pipeline_option_widgets={},
            pipeline_rows=[velocity, metrics, report],
            pipeline_options={},
        )
        controller = PipelineLibraryController(app)

        controller.update_summary = Mock()

        controller.set_visibility("report", True)

        self.assertEqual(
            {"velocity": False, "metrics": False, "report": True},
            app.settings.pipeline_visibility,
        )
        self.assertEqual({"velocity", "metrics"}, app.pipeline_required_names)

        controller.set_visibility("report", False)

        self.assertEqual(
            {"velocity": False, "metrics": False, "report": False},
            app.settings.pipeline_visibility,
        )
        self.assertEqual(set(), app.pipeline_required_names)

        app.settings_controller.save_preferences.reset_mock()
        controller.set_all_visibility(True)
        self.assertEqual({"velocity": True, "metrics": True, "report": True},
                         app.settings.pipeline_visibility)
        self.assertEqual({}, app.settings.pipeline_options)
        app.settings_controller.save_preferences.assert_called_once_with()

    def test_option_dependencies_replan_required_visible_pipelines(self) -> None:
        topology = _descriptor(
            "topology_core",
            produces=("prepared_topology",),
            visibility="hidden",
        )
        base_velocity = _descriptor(
            "velocity",
            produces=("velocity", "cardiac_cycles"),
            visibility="hidden",
        )
        gradient = _descriptor(
            "spatial_gradient_moment0",
            requires=("cardiac_cycles", "prepared_topology"),
            produces=("spatial_gradient_edges",),
        )
        velocity_analysis = _descriptor(
            "velocity_analysis",
            requires=("velocity", "prepared_topology"),
            produces=("velocity_analysis",),
            visibility="hidden",
        )
        bvr = _descriptor(
            "blood_volume_rate",
            requires=("velocity_analysis",),
            options=(
                PipelineOption(
                    "gradient_edges",
                    "Gradient",
                    dag_requires=("spatial_gradient_edges",),
                ),
                PipelineOption(
                    "masked_edges",
                    "Masked",
                ),
            ),
        )
        catalog = {
            item.name: item
            for item in (
                base_velocity,
                topology,
                gradient,
                velocity_analysis,
                bvr,
            )
        }
        app = self.make_app(
            pipeline_catalog=catalog,
            pipeline_dag=PipelineDAG(catalog.values()),
            pipeline_rows=[gradient, bvr],
            pipeline_visibility={
                "spatial_gradient_moment0": False,
                "blood_volume_rate": True,
            },
            pipeline_options={
                "blood_volume_rate": {
                    "gradient_edges": True,
                    "masked_edges": True,
                }
            },
            pipeline_visibility_vars={},
            pipeline_row_widgets={},
            pipeline_option_vars={"blood_volume_rate": {}},
            pipeline_option_widgets={},
        )
        controller = PipelineLibraryController(app)

        controller.update_summary = Mock()
        controller._refresh_required_pipelines()

        self.assertEqual(
            {"spatial_gradient_moment0"},
            app.pipeline_required_names,
        )
        controller.set_option_visibility(
            "blood_volume_rate",
            "gradient_edges",
            False,
        )
        self.assertEqual(set(), app.pipeline_required_names)
        self.assertFalse(app.settings.pipeline_visibility["spatial_gradient_moment0"])

    def test_child_option_selection_follows_declared_requirements(self) -> None:
        options = (
            PipelineOption(
                "profiles",
                "Profiles",
                requires=("per_beat",),
            ),
            PipelineOption("per_beat", "Per beat"),
            PipelineOption(
                "quadrants",
                "Quadrants",
                requires=("per_beat",),
            ),
        )

        self.assertEqual(
            ("per_beat", "quadrants"),
            option_selection_names(options, "quadrants", enabled=True),
        )
        self.assertEqual(
            ("profiles", "per_beat", "quadrants"),
            option_selection_names(options, "per_beat", enabled=False),
        )

        pipeline = _descriptor("velocity", options=options)
        app = self.make_app(
            pipeline_catalog={"velocity": pipeline},
            pipeline_options={
                "velocity": {
                    "profiles": False,
                    "per_beat": False,
                    "quadrants": False,
                }
            },
            pipeline_option_vars={"velocity": {}},
        )
        controller = PipelineLibraryController(app)

        controller.update_summary = Mock()

        controller.set_option_visibility("velocity", "quadrants", True)
        self.assertEqual(
            {"profiles": False, "per_beat": True, "quadrants": True},
            app.settings.pipeline_options["velocity"],
        )

        controller.set_option_visibility("velocity", "per_beat", False)
        self.assertEqual(
            {"profiles": False, "per_beat": False, "quadrants": False},
            app.settings.pipeline_options["velocity"],
        )

    def test_waveform_segment_substeps_follow_upstream_selection(self) -> None:
        velocity = _descriptor(
            "velocity_analysis",
            options=(
                PipelineOption("segments", "Segments"),
                PipelineOption(
                    "quadrants",
                    "Quadrants",
                    requires=("segments",),
                ),
            ),
        )
        shape = _descriptor(
            "waveform_shape_metrics",
            options=(
                PipelineOption("segments", "Segments"),
                PipelineOption("quadrants", "Quadrants"),
            ),
        )
        app = self.make_app(
            pipeline_catalog={
                velocity.name: velocity,
                shape.name: shape,
            },
            pipeline_options={
                "velocity_analysis": {
                    "segments": True,
                    "quadrants": True,
                },
                "waveform_shape_metrics": {
                    "segments": True,
                    "quadrants": True,
                },
            },
            pipeline_option_vars={"velocity_analysis": {}, "waveform_shape_metrics": {}},
        )
        controller = PipelineLibraryController(app)

        controller.update_summary = Mock()

        controller.set_option_visibility("velocity_analysis", "segments", False)

        self.assertFalse(
            app.settings.pipeline_options["velocity_analysis"]["segments"]
        )
        self.assertFalse(
            app.settings.pipeline_options["waveform_shape_metrics"]["segments"]
        )
        self.assertFalse(
            app.settings.pipeline_options["velocity_analysis"]["quadrants"]
        )
        self.assertFalse(
            app.settings.pipeline_options["waveform_shape_metrics"]["quadrants"]
        )

        controller.set_option_visibility("waveform_shape_metrics", "segments", True)

        self.assertTrue(
            app.settings.pipeline_options["velocity_analysis"]["segments"]
        )
        self.assertTrue(
            app.settings.pipeline_options["waveform_shape_metrics"]["segments"]
        )

    def test_pdf_report_does_not_mutate_shape_options(self) -> None:
        velocity = _descriptor(
            "velocity_analysis",
            options=(PipelineOption("segments", "Segments"),),
        )
        shape = _descriptor(
            "waveform_shape_metrics",
            options=(PipelineOption("segments", "Segments"),),
        )
        report = _descriptor(
            "pdf_report",
            requires=("velocity_analysis", "waveform_shape_metrics"),
        )
        app = self.make_app(
            pipeline_catalog={
                item.name: item for item in (velocity, shape, report)
            },
            pipeline_dag=PipelineDAG((velocity, shape, report)),
            pipeline_visibility={
                "velocity_analysis": False,
                "waveform_shape_metrics": False,
                "pdf_report": False,
            },
            pipeline_visibility_vars={},
            pipeline_option_widgets={},
            pipeline_options={
                "velocity_analysis": {"segments": False},
                "waveform_shape_metrics": {"segments": False},
            },
            pipeline_option_vars={
                "velocity_analysis": {},
                "waveform_shape_metrics": {},
            },
        )
        controller = PipelineLibraryController(app)


        controller.update_summary = Mock()

        controller.set_visibility("pdf_report", True)

        self.assertEqual(
            {"segments": False},
            app.settings.pipeline_options["waveform_shape_metrics"],
        )

    def test_stored_waveform_segment_selection_is_normalized_downstream(self) -> None:
        velocity = _descriptor(
            "velocity_analysis",
            options=(PipelineOption("segments", "Segments"),),
        )
        shape = _descriptor(
            "waveform_shape_metrics",
            options=(
                PipelineOption("segments", "Segments"),
                PipelineOption(
                    "quadrants",
                    "Quadrants",
                    requires=("segments",),
                ),
            ),
        )
        app = self.make_app(
            pipeline_catalog={
                velocity.name: velocity,
                shape.name: shape,
            },
            pipeline_options={
                "velocity_analysis": {"segments": False},
                "waveform_shape_metrics": {"segments": True, "quadrants": True},
            },
        )
        controller = PipelineLibraryController(app)


        controller.sync_options([velocity, shape])

        self.assertTrue(
            app.settings.pipeline_options["velocity_analysis"]["segments"]
        )
        self.assertTrue(
            app.settings.pipeline_options["waveform_shape_metrics"]["segments"]
        )
        self.assertTrue(
            app.settings.pipeline_options["waveform_shape_metrics"]["quadrants"]
        )


if __name__ == "__main__":
    unittest.main()
