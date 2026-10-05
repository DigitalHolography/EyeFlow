"""Tests for dependency-aware pipeline-library selection."""

from __future__ import annotations

import unittest
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
            _descriptor("waveform_velocity"),
        ]

        self.assertEqual(
            [
                "waveform_velocity",
                "waveform_shape_metrics",
                "pdf_report",
                "another_pipeline",
            ],
            [item.name for item in sorted(rows, key=pipeline_ui_sort_key)],
        )

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
        app = SimpleNamespace(
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
            pipeline_option_visibility={},
        )
        controller = PipelineLibraryController(app)
        controller.persist_visibility = Mock()
        controller.update_summary = Mock()

        controller.set_visibility("report", True)

        self.assertEqual(
            {"velocity": False, "metrics": False, "report": True},
            app.pipeline_visibility,
        )
        self.assertEqual({"velocity", "metrics"}, app.pipeline_required_names)

        controller.set_visibility("report", False)

        self.assertEqual(
            {"velocity": False, "metrics": False, "report": False},
            app.pipeline_visibility,
        )
        self.assertEqual(set(), app.pipeline_required_names)

    def test_option_dependencies_replan_required_visible_pipelines(self) -> None:
        topology = _descriptor(
            "topology_core",
            produces=("prepared_topology",),
            visibility="hidden",
        )
        retinal_velocity = _descriptor(
            "retinal_velocity",
            produces=("retinal_velocity", "cardiac_cycles"),
            visibility="hidden",
        )
        gradient = _descriptor(
            "spatial_gradient_moment0",
            requires=("cardiac_cycles", "prepared_topology"),
            produces=("spatial_gradient_edges",),
        )
        velocity = _descriptor(
            "waveform_velocity",
            requires=("retinal_velocity", "prepared_topology"),
            produces=("waveform_velocity",),
            visibility="hidden",
        )
        bvr = _descriptor(
            "blood_volume_rate",
            requires=("waveform_velocity",),
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
            for item in (retinal_velocity, topology, gradient, velocity, bvr)
        }
        app = SimpleNamespace(
            pipeline_catalog=catalog,
            pipeline_dag=PipelineDAG(catalog.values()),
            pipeline_rows=[gradient, bvr],
            pipeline_visibility={
                "spatial_gradient_moment0": False,
                "blood_volume_rate": True,
            },
            pipeline_option_visibility={
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
        controller.persist_options = Mock()
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
        self.assertFalse(app.pipeline_visibility["spatial_gradient_moment0"])

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
        app = SimpleNamespace(
            pipeline_catalog={"velocity": pipeline},
            pipeline_option_visibility={
                "velocity": {
                    "profiles": False,
                    "per_beat": False,
                    "quadrants": False,
                }
            },
            pipeline_option_vars={"velocity": {}},
        )
        controller = PipelineLibraryController(app)
        controller.persist_options = Mock()
        controller.update_summary = Mock()

        controller.set_option_visibility("velocity", "quadrants", True)
        self.assertEqual(
            {"profiles": False, "per_beat": True, "quadrants": True},
            app.pipeline_option_visibility["velocity"],
        )

        controller.set_option_visibility("velocity", "per_beat", False)
        self.assertEqual(
            {"profiles": False, "per_beat": False, "quadrants": False},
            app.pipeline_option_visibility["velocity"],
        )

    def test_waveform_segment_substeps_follow_upstream_selection(self) -> None:
        velocity = _descriptor(
            "waveform_velocity",
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
        app = SimpleNamespace(
            pipeline_catalog={
                velocity.name: velocity,
                shape.name: shape,
            },
            pipeline_option_visibility={
                "waveform_velocity": {
                    "segments": True,
                    "quadrants": True,
                },
                "waveform_shape_metrics": {
                    "segments": True,
                    "quadrants": True,
                },
            },
            pipeline_option_vars={"waveform_velocity": {}, "waveform_shape_metrics": {}},
        )
        controller = PipelineLibraryController(app)
        controller.persist_options = Mock()
        controller.update_summary = Mock()

        controller.set_option_visibility("waveform_velocity", "segments", False)

        self.assertFalse(
            app.pipeline_option_visibility["waveform_velocity"]["segments"]
        )
        self.assertFalse(
            app.pipeline_option_visibility["waveform_shape_metrics"]["segments"]
        )
        self.assertFalse(
            app.pipeline_option_visibility["waveform_velocity"]["quadrants"]
        )
        self.assertFalse(
            app.pipeline_option_visibility["waveform_shape_metrics"]["quadrants"]
        )

        controller.set_option_visibility("waveform_shape_metrics", "segments", True)

        self.assertTrue(
            app.pipeline_option_visibility["waveform_velocity"]["segments"]
        )
        self.assertTrue(
            app.pipeline_option_visibility["waveform_shape_metrics"]["segments"]
        )

    def test_pdf_report_does_not_mutate_shape_options(self) -> None:
        velocity = _descriptor(
            "waveform_velocity",
            options=(PipelineOption("segments", "Segments"),),
        )
        shape = _descriptor(
            "waveform_shape_metrics",
            options=(PipelineOption("segments", "Segments"),),
        )
        report = _descriptor(
            "pdf_report",
            requires=("waveform_velocity", "waveform_shape_metrics"),
        )
        app = SimpleNamespace(
            pipeline_catalog={
                item.name: item for item in (velocity, shape, report)
            },
            pipeline_dag=PipelineDAG((velocity, shape, report)),
            pipeline_visibility={
                "waveform_velocity": False,
                "waveform_shape_metrics": False,
                "pdf_report": False,
            },
            pipeline_visibility_vars={},
            pipeline_option_widgets={},
            pipeline_option_visibility={
                "waveform_velocity": {"segments": False},
                "waveform_shape_metrics": {"segments": False},
            },
            pipeline_option_vars={
                "waveform_velocity": {},
                "waveform_shape_metrics": {},
            },
        )
        controller = PipelineLibraryController(app)
        controller.persist_visibility = Mock()
        controller.persist_options = Mock()
        controller.update_summary = Mock()

        controller.set_visibility("pdf_report", True)

        self.assertEqual(
            {"segments": False},
            app.pipeline_option_visibility["waveform_shape_metrics"],
        )

    def test_stored_waveform_segment_selection_is_normalized_downstream(self) -> None:
        velocity = _descriptor(
            "waveform_velocity",
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
        app = SimpleNamespace(
            pipeline_catalog={
                velocity.name: velocity,
                shape.name: shape,
            },
            settings_store=SimpleNamespace(
                load_pipeline_options=lambda: {
                    "waveform_velocity": {"segments": False},
                    "waveform_shape_metrics": {
                        "segments": True,
                        "quadrants": True,
                    },
                }
            ),
        )
        controller = PipelineLibraryController(app)
        controller.persist_options = Mock()

        controller.sync_options([velocity, shape])

        self.assertTrue(
            app.pipeline_option_visibility["waveform_velocity"]["segments"]
        )
        self.assertTrue(
            app.pipeline_option_visibility["waveform_shape_metrics"]["segments"]
        )
        self.assertTrue(
            app.pipeline_option_visibility["waveform_shape_metrics"]["quadrants"]
        )


if __name__ == "__main__":
    unittest.main()
