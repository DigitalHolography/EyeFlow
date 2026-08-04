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

    def test_selecting_downstream_checks_upstream_and_unchecking_reverses(self) -> None:
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
        )
        controller = PipelineLibraryController(app)
        controller.persist_visibility = Mock()
        controller.update_summary = Mock()

        controller.set_visibility("report", True)

        self.assertEqual(
            {"velocity": True, "metrics": True, "report": True},
            app.pipeline_visibility,
        )

        controller.set_visibility("metrics", False)

        self.assertEqual(
            {"velocity": True, "metrics": False, "report": False},
            app.pipeline_visibility,
        )

    def test_child_option_selection_follows_declared_requirements(self) -> None:
        options = (
            PipelineOption(
                "profiles",
                "Profiles",
                requires=("per_beat",),
            ),
            PipelineOption("per_beat", "Per beat"),
            PipelineOption(
                "hemifield",
                "Hemifield",
                requires=("per_beat",),
            ),
        )

        self.assertEqual(
            ("per_beat", "hemifield"),
            option_selection_names(options, "hemifield", enabled=True),
        )
        self.assertEqual(
            ("profiles", "per_beat", "hemifield"),
            option_selection_names(options, "per_beat", enabled=False),
        )

        pipeline = _descriptor("velocity", options=options)
        app = SimpleNamespace(
            pipeline_catalog={"velocity": pipeline},
            pipeline_option_visibility={
                "velocity": {
                    "profiles": False,
                    "per_beat": False,
                    "hemifield": False,
                }
            },
            pipeline_option_vars={"velocity": {}},
        )
        controller = PipelineLibraryController(app)
        controller.persist_options = Mock()
        controller.update_summary = Mock()

        controller.set_option_visibility("velocity", "hemifield", True)
        self.assertEqual(
            {"profiles": False, "per_beat": True, "hemifield": True},
            app.pipeline_option_visibility["velocity"],
        )

        controller.set_option_visibility("velocity", "per_beat", False)
        self.assertEqual(
            {"profiles": False, "per_beat": False, "hemifield": False},
            app.pipeline_option_visibility["velocity"],
        )

    def test_waveform_segment_substeps_follow_upstream_selection(self) -> None:
        velocity = _descriptor(
            "waveform_velocity",
            options=(
                PipelineOption("segments", "Segments"),
                PipelineOption("per_beat", "Per beat"),
            ),
        )
        shape = _descriptor(
            "waveform_shape_metrics",
            options=(
                PipelineOption("per_beat", "Per beat"),
                PipelineOption(
                    "segments",
                    "Segments",
                    requires=("per_beat",),
                ),
            ),
        )
        app = SimpleNamespace(
            pipeline_catalog={
                velocity.name: velocity,
                shape.name: shape,
            },
            pipeline_option_visibility={
                "waveform_velocity": {"segments": True, "per_beat": True},
                "waveform_shape_metrics": {"segments": True, "per_beat": True},
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

        controller.set_option_visibility("waveform_shape_metrics", "segments", True)

        self.assertTrue(
            app.pipeline_option_visibility["waveform_velocity"]["segments"]
        )
        self.assertTrue(
            app.pipeline_option_visibility["waveform_shape_metrics"]["segments"]
        )

        controller.set_option_visibility("waveform_velocity", "per_beat", False)

        self.assertFalse(
            app.pipeline_option_visibility["waveform_shape_metrics"]["per_beat"]
        )
        self.assertFalse(
            app.pipeline_option_visibility["waveform_shape_metrics"]["segments"]
        )

        controller.set_option_visibility("waveform_shape_metrics", "per_beat", True)

        self.assertTrue(
            app.pipeline_option_visibility["waveform_velocity"]["per_beat"]
        )

    def test_pdf_report_keeps_upstream_per_beat_products_enabled(self) -> None:
        velocity = _descriptor(
            "waveform_velocity",
            options=(PipelineOption("per_beat", "Per beat"),),
        )
        shape = _descriptor(
            "waveform_shape_metrics",
            options=(PipelineOption("per_beat", "Per beat"),),
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
                "waveform_velocity": {"per_beat": False},
                "waveform_shape_metrics": {"per_beat": False},
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

        self.assertTrue(
            app.pipeline_option_visibility["waveform_velocity"]["per_beat"]
        )
        self.assertTrue(
            app.pipeline_option_visibility["waveform_shape_metrics"]["per_beat"]
        )

        controller.set_option_visibility("waveform_velocity", "per_beat", False)

        self.assertTrue(
            app.pipeline_option_visibility["waveform_velocity"]["per_beat"]
        )
        self.assertTrue(
            app.pipeline_option_visibility["waveform_shape_metrics"]["per_beat"]
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
                    "hemifield",
                    "Hemifield",
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
                        "hemifield": True,
                    },
                }
            ),
        )
        controller = PipelineLibraryController(app)
        controller.persist_options = Mock()

        controller.sync_options([velocity, shape])

        self.assertFalse(
            app.pipeline_option_visibility["waveform_shape_metrics"]["segments"]
        )
        self.assertFalse(
            app.pipeline_option_visibility["waveform_shape_metrics"]["hemifield"]
        )


if __name__ == "__main__":
    unittest.main()
