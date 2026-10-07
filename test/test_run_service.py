"""Tests for shared GUI/CLI run orchestration and direct outputs."""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
import zipfile
from contextlib import redirect_stderr, redirect_stdout
from io import StringIO
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock, patch

import h5py

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from app_settings import (  # noqa: E402
    AppSettingsStore,
    normalize_pipeline_options,
)
from input_output.archives import extracted_zip_tree  # noqa: E402
from input_output.output_manager import OutputType  # noqa: E402
from pipeline_engine import (  # noqa: E402
    PipelineDescriptor,
    PipelineOption,
    ProcessPipeline,
)
from pipeline_engine.run_service import (  # noqa: E402
    execute_run,
    expand_run_inputs,
    resolve_run_spec,
)
from ui.controllers.run import RunController  # noqa: E402


def _write_input(root: Path, stem: str = "scan") -> Path:
    root.mkdir(parents=True, exist_ok=True)
    holo = root / f"{stem}.holo"
    holo.write_text("holo", encoding="utf-8")
    hd = root / stem / f"{stem}_HD" / "h5" / f"{stem}_HD_output.h5"
    dv = root / stem / f"{stem}_DV" / "h5" / f"{stem}_DV.h5"
    hd.parent.mkdir(parents=True)
    dv.parent.mkdir(parents=True)
    with h5py.File(hd, "w"):
        pass
    with h5py.File(dv, "w"):
        pass
    return holo


class _NoopPipeline(ProcessPipeline):
    name = "sample"
    description = "sample"
    available = True
    requires = []
    missing_deps = []

    def run(self, ctx):
        del ctx
        return None


def test_run_controls_include_velocity_estimator_toggle() -> None:
    minimal_run = Mock()
    advanced_run = Mock()
    estimator_buttons = [Mock(), Mock()]
    controller = RunController.__new__(RunController)
    controller.app = SimpleNamespace(
        minimal_run_button=minimal_run,
        advanced_run_button=advanced_run,
        velocity_estimation_widgets=estimator_buttons,
    )

    controller._set_run_controls_enabled(False)

    minimal_run.configure.assert_called_once_with(state="disabled")
    advanced_run.configure.assert_called_once_with(state="disabled")
    for button in estimator_buttons:
        button.state.assert_called_once_with(["disabled"])


def _descriptor(*, visibility: str = "visible") -> PipelineDescriptor:
    return PipelineDescriptor(
        name="sample",
        description="sample",
        available=True,
        visibility=visibility,
        pipeline_factory=_NoopPipeline,
    )


def _named_descriptor(name: str) -> PipelineDescriptor:
    return PipelineDescriptor(
        name=name,
        description=name,
        available=True,
        visibility="visible",
        pipeline_factory=_NoopPipeline,
    )


class RunServiceTests(unittest.TestCase):
    def test_success_replaces_existing_output_with_direct_run(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            holo = _write_input(root)
            spec = resolve_run_spec(
                input_paths=[holo],
                target_names=["sample"],
                pipelines=[_descriptor()],
            )
            final_dir = spec.requests[0].output_manager.layout.ef_dir
            final_dir.mkdir(parents=True)
            (final_dir / "old.txt").write_text("old", encoding="utf-8")

            def fake_run(*, output_manager, **_kwargs):
                output_manager.prepare()
                (output_manager.layout.ef_dir / "new.txt").write_text(
                    "new", encoding="utf-8"
                )
                return output_manager.path_for(OutputType.H5)

            with patch(
                "pipeline_engine.run_service.run_pipelines_to_output",
                side_effect=fake_run,
            ):
                result = execute_run(spec)

            self.assertTrue(result.succeeded)
            self.assertFalse((final_dir / "old.txt").exists())
            self.assertEqual("new", (final_dir / "new.txt").read_text(encoding="utf-8"))
            self.assertFalse(list(final_dir.parent.glob(".*eyeflow-staging-*")))

    def test_failed_direct_run_leaves_partial_output(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            holo = _write_input(root)
            spec = resolve_run_spec(
                input_paths=[holo],
                target_names=["sample"],
                pipelines=[_descriptor()],
            )
            final_dir = spec.requests[0].output_manager.layout.ef_dir
            final_dir.mkdir(parents=True)
            marker = final_dir / "old.txt"
            marker.write_text("old", encoding="utf-8")

            def fake_failure(*, output_manager, **_kwargs):
                output_manager.prepare()
                (output_manager.layout.ef_dir / "partial.txt").write_text(
                    "partial", encoding="utf-8"
                )
                raise RuntimeError("analysis failed")

            with patch(
                "pipeline_engine.run_service.run_pipelines_to_output",
                side_effect=fake_failure,
            ):
                result = execute_run(spec)

            self.assertEqual(1, len(result.failures))
            self.assertFalse(marker.exists())
            self.assertEqual(
                "partial", (final_dir / "partial.txt").read_text(encoding="utf-8")
            )
            self.assertFalse(list(final_dir.parent.glob(".*eyeflow-staging-*")))

    def test_direct_run_refuses_to_replace_non_directory_output(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            holo = _write_input(root)
            spec = resolve_run_spec(
                input_paths=[holo],
                target_names=["sample"],
                pipelines=[_descriptor()],
            )
            final_dir = spec.requests[0].output_manager.layout.ef_dir
            final_dir.parent.mkdir(parents=True, exist_ok=True)
            final_dir.write_text("user file", encoding="utf-8")

            def fake_run(*, output_manager, **_kwargs):
                output_manager.prepare()

            with patch(
                "pipeline_engine.run_service.run_pipelines_to_output",
                side_effect=fake_run,
            ):
                result = execute_run(spec)

            self.assertEqual(1, len(result.failures))
            self.assertEqual("user file", final_dir.read_text(encoding="utf-8"))

    def test_hidden_pipeline_cannot_be_selected_directly(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            holo = _write_input(Path(temp_dir))
            with self.assertRaisesRegex(ValueError, "hidden"):
                resolve_run_spec(
                    input_paths=[holo],
                    target_names=["sample"],
                    pipelines=[_descriptor(visibility="hidden")],
                )

    def test_hidden_pipeline_can_satisfy_visible_target_dependency(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            holo = _write_input(Path(temp_dir))
            hidden = PipelineDescriptor(
                name="preparation",
                description="hidden preparation",
                available=True,
                visibility="hidden",
                dag_produces=("prepared",),
                pipeline_factory=_NoopPipeline,
            )
            visible = _descriptor()
            visible.dag_requires = ("prepared",)

            spec = resolve_run_spec(
                input_paths=[holo],
                target_names=["sample"],
                pipelines=[hidden, visible],
            )

            self.assertEqual(("preparation", "sample"), spec.plan.names)

    def test_velocity_estimation_method_is_validated_and_stored(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            holo = _write_input(Path(temp_dir))
            default_spec = resolve_run_spec(
                input_paths=[holo],
                target_names=["sample"],
                pipelines=[_descriptor()],
            )
            band_spec = resolve_run_spec(
                input_paths=[holo],
                target_names=["sample"],
                pipelines=[_descriptor()],
                velocity_estimation_method="frequency_bands",
                band_ratio_frequency_scale_hz=2.5,
            )

            self.assertEqual(
                "doppler_moments",
                default_spec.velocity_estimation_method,
            )
            self.assertEqual(
                "frequency_bands",
                band_spec.velocity_estimation_method,
            )
            self.assertEqual(2.5, band_spec.band_ratio_frequency_scale_hz)

            with self.assertRaisesRegex(ValueError, "velocity_estimation_method"):
                resolve_run_spec(
                    input_paths=[holo],
                    target_names=["sample"],
                    pipelines=[_descriptor()],
                    velocity_estimation_method="unknown",
                )
            with self.assertRaisesRegex(
                ValueError,
                "band_ratio_frequency_scale_hz",
            ):
                resolve_run_spec(
                    input_paths=[holo],
                    target_names=["sample"],
                    pipelines=[_descriptor()],
                    band_ratio_frequency_scale_hz=0.0,
                )

    def test_frequency_band_method_allows_physical_velocity_pipelines(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            holo = _write_input(Path(temp_dir))
            for pipeline_name in (
                "absolute_waveform_metrics",
                "blood_volume_rate",
            ):
                with self.subTest(pipeline=pipeline_name):
                    spec = resolve_run_spec(
                        input_paths=[holo],
                        target_names=[pipeline_name],
                        pipelines=[_named_descriptor(pipeline_name)],
                        velocity_estimation_method="frequency_bands",
                    )

                    self.assertEqual((pipeline_name,), spec.plan.targets)
                    self.assertEqual(
                        "frequency_bands",
                        spec.velocity_estimation_method,
                    )

    def test_gui_reads_velocity_method_when_building_run_spec(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            holo = _write_input(Path(temp_dir))
            progress_controller = SimpleNamespace(reset_run_log=Mock())
            settings_store = SimpleNamespace(
                load_velocity_estimation_method=Mock(
                    return_value="frequency_bands"
                ),
                load_band_ratio_frequency_scale_hz=Mock(return_value=3.0),
            )
            app = SimpleNamespace(
                input_controller=SimpleNamespace(
                    selected_holo_paths=Mock(return_value=[holo])
                ),
                pipeline_library_controller=SimpleNamespace(
                    selected_target_pipeline_names=Mock(return_value=["sample"]),
                    selected_pipeline_options=Mock(return_value={}),
                ),
                pipeline_catalog={"sample": _descriptor()},
                settings_store=settings_store,
                progress_controller=progress_controller,
                ui_services=SimpleNamespace(dialogs=Mock()),
            )
            controller = RunController.__new__(RunController)
            controller.app = app

            spec = controller._build_run_spec()

            self.assertIsNotNone(spec)
            assert spec is not None
            self.assertEqual(
                "frequency_bands",
                spec.velocity_estimation_method,
            )
            self.assertEqual(3.0, spec.band_ratio_frequency_scale_hz)
            settings_store.load_velocity_estimation_method.assert_called_once_with()
            settings_store.load_band_ratio_frequency_scale_hz.assert_called_once_with()
            progress_controller.reset_run_log.assert_called_once()

    def test_pipeline_options_default_validate_and_preserve_empty_selection(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            holo = _write_input(Path(temp_dir))
            descriptor = _descriptor()
            descriptor.options = (
                PipelineOption("default", "Default"),
                PipelineOption("opt_in", "Opt in", default_enabled=False),
            )

            default_spec = resolve_run_spec(
                input_paths=[holo],
                target_names=["sample"],
                pipelines=[descriptor],
            )
            empty_spec = resolve_run_spec(
                input_paths=[holo],
                target_names=["sample"],
                pipelines=[descriptor],
                pipeline_options={"sample": ()},
            )

            self.assertEqual(("default",), default_spec.pipeline_options["sample"])
            self.assertEqual((), empty_spec.pipeline_options["sample"])
            result = execute_run(empty_spec)
            self.assertTrue(result.succeeded)
            with h5py.File(result.outputs[0], "r") as output_h5:
                self.assertEqual(
                    {"sample": []},
                    json.loads(output_h5.attrs["pipeline_options"]),
                )
            with self.assertRaisesRegex(ValueError, "Unknown option"):
                resolve_run_spec(
                    input_paths=[holo],
                    target_names=["sample"],
                    pipelines=[descriptor],
                    pipeline_options={"sample": ("missing",)},
                )

    def test_pipeline_option_settings_are_normalized_and_persisted(self) -> None:
        options = {
            "velocity_analysis": (
                PipelineOption("segments", "Segments"),
                PipelineOption("quadrants", "Quadrants"),
            )
        }
        normalized, changed = normalize_pipeline_options(
            options,
            {
                "velocity_analysis": {"segments": False, "removed": True},
                "removed_pipeline": {"old": True},
            },
        )

        self.assertTrue(changed)
        self.assertEqual(
            {"velocity_analysis": {"segments": False, "quadrants": True}},
            normalized,
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            store = AppSettingsStore(
                path=Path(temp_dir) / "settings.json",
                default_template_path=None,
            )
            store.save_pipeline_options(normalized)
            self.assertEqual(normalized, store.load_pipeline_options())

    def test_duplicate_output_destinations_are_rejected_before_execution(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            holo = _write_input(Path(temp_dir))
            with self.assertRaisesRegex(ValueError, "same EyeFlow output"):
                resolve_run_spec(
                    input_paths=[holo, holo],
                    target_names=["sample"],
                    pipelines=[_descriptor()],
                )

    def test_runtime_output_does_not_contain_angioeye_trim_attribute(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            holo = _write_input(Path(temp_dir))
            spec = resolve_run_spec(
                input_paths=[holo],
                target_names=["sample"],
                pipelines=[_descriptor()],
            )

            result = execute_run(spec)

            self.assertTrue(result.succeeded)
            with h5py.File(result.outputs[0], "r") as output_h5:
                self.assertNotIn("trim_h5source", output_h5.attrs)
                self.assertEqual(["sample"], list(output_h5.attrs["pipeline_targets"]))
                self.assertEqual({}, json.loads(output_h5.attrs["pipeline_options"]))
                self.assertEqual(
                    "doppler_moments",
                    output_h5.attrs["velocity_estimation_method"],
                )
                self.assertEqual(
                    "physical_velocity",
                    output_h5.attrs["velocity_quantity"],
                )
                self.assertEqual("mm/s", output_h5.attrs["velocity_unit"])

    def test_frequency_band_method_is_persisted_in_output_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            holo = _write_input(Path(temp_dir))
            spec = resolve_run_spec(
                input_paths=[holo],
                target_names=["sample"],
                pipelines=[_descriptor()],
                velocity_estimation_method="frequency_bands",
                band_ratio_frequency_scale_hz=2.5,
            )

            result = execute_run(spec)

            self.assertTrue(result.succeeded)
            with h5py.File(result.outputs[0], "r") as output_h5:
                self.assertEqual(
                    "frequency_bands",
                    output_h5.attrs["velocity_estimation_method"],
                )
                self.assertEqual(
                    "physical_velocity",
                    output_h5.attrs["velocity_quantity"],
                )
                self.assertEqual("mm/s", output_h5.attrs["velocity_unit"])
                self.assertEqual(
                    2.5,
                    output_h5.attrs["band_ratio_frequency_scale_hz"],
                )
                self.assertEqual(
                    "linear_origin",
                    output_h5.attrs["band_ratio_calibration_model"],
                )
                self.assertEqual(
                    "eyeflow_setting",
                    output_h5.attrs["band_ratio_calibration_source"],
                )
                self.assertEqual(
                    "1",
                    output_h5.attrs["band_ratio_calibration_version"],
                )
                self.assertAlmostEqual(
                    8.52e-7,
                    output_h5.attrs["laser_wavelength_m"],
                )
                self.assertAlmostEqual(
                    0.76,
                    output_h5.attrs["numerical_aperture"],
                )
                self.assertEqual(
                    "/band_0_3000_9000",
                    output_h5.attrs["band_lf_source_path"],
                )
                self.assertEqual(
                    "/band_1_9000_18000",
                    output_h5.attrs["band_hf_source_path"],
                )

    def test_cli_reports_zip_creation_failure_with_nonzero_status(self) -> None:
        import cli

        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            holo = _write_input(root / "input")
            pipelines_file = root / "pipelines.txt"
            pipelines_file.write_text("sample\n", encoding="utf-8")
            output_root = root / "output"
            registry = {"sample": _descriptor()}

            with (
                patch("cli._build_pipeline_registry", return_value=registry),
                patch("cli._zip_output_dir", side_effect=OSError("zip failed")),
                redirect_stdout(StringIO()),
                redirect_stderr(StringIO()),
            ):
                status = cli.run_cli(
                    holo,
                    pipelines_file,
                    output_root,
                    zip_outputs=True,
                )

            self.assertEqual(1, status)

    def test_cli_uses_enabled_settings_and_source_adjacent_output_defaults(self) -> None:
        import cli

        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            holo = _write_input(root / "input")
            store = AppSettingsStore(
                path=root / "settings.json",
                default_template_path=None,
            )
            store.save(
                {
                    "pipeline_visibility": {"sample": True},
                    "velocity_estimation_method": "frequency_bands",
                    "band_ratio_frequency_scale_hz": 4.0,
                }
            )

            with (
                patch(
                    "cli._build_pipeline_registry",
                    return_value={"sample": _descriptor()},
                ),
                patch("cli.AppSettingsStore", return_value=store),
                redirect_stdout(StringIO()),
                redirect_stderr(StringIO()),
            ):
                status = cli.run_cli(holo)

            expected = root / "input" / "scan" / "scan_EF" / "h5" / "scan_EF.h5"
            self.assertEqual(0, status)
            self.assertTrue(expected.is_file())
            with h5py.File(expected, "r") as output_h5:
                self.assertEqual(
                    "frequency_bands",
                    output_h5.attrs["velocity_estimation_method"],
                )
                self.assertEqual(
                    4.0,
                    output_h5.attrs["band_ratio_frequency_scale_hz"],
                )

    def test_cli_requires_argument_when_no_pipeline_is_enabled(self) -> None:
        import cli

        with tempfile.TemporaryDirectory() as temp_dir:
            store = AppSettingsStore(
                path=Path(temp_dir) / "settings.json",
                default_template_path=None,
            )
            store.save({"pipeline_visibility": {"sample": False}})

            with self.assertRaisesRegex(ValueError, "provide --pipelines"):
                cli._load_configured_pipeline_targets(
                    {"sample": _descriptor()},
                    settings_store=store,
                )

    def test_zip_input_default_output_survives_temporary_extraction(self) -> None:
        import cli

        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            payload = root / "payload"
            _write_input(payload)
            archive_path = root / "input.zip"
            with zipfile.ZipFile(archive_path, "w") as archive:
                for path in sorted(payload.rglob("*")):
                    if path.is_file():
                        archive.write(path, path.relative_to(payload))
            pipelines_file = root / "pipelines.txt"
            pipelines_file.write_text("sample\n", encoding="utf-8")

            with (
                patch(
                    "cli._build_pipeline_registry",
                    return_value={"sample": _descriptor()},
                ),
                redirect_stdout(StringIO()),
                redirect_stderr(StringIO()),
            ):
                status = cli.run_cli(archive_path, pipelines_file=pipelines_file)

            expected = root / "scan" / "scan_EF" / "h5" / "scan_EF.h5"
            self.assertEqual(0, status)
            self.assertTrue(expected.is_file())

    def test_zip_creation_failure_preserves_previous_archive(self) -> None:
        import cli

        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            source = root / "source"
            source.mkdir()
            (source / "result.txt").write_text("result", encoding="utf-8")
            archive = root / "outputs.zip"
            archive.write_text("previous", encoding="utf-8")

            with patch("cli.create_zip_from_tree", side_effect=OSError("zip failed")):
                with self.assertRaisesRegex(OSError, "zip failed"):
                    cli._zip_output_dir(source, archive)

            self.assertEqual("previous", archive.read_text(encoding="utf-8"))
            self.assertFalse(list(root.glob(".*eyeflow-staging-*")))

    def test_recursive_folder_expansion_is_sorted(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            second = root / "b" / "second.holo"
            first = root / "a" / "first.holo"
            second.parent.mkdir()
            first.parent.mkdir()
            second.write_text("holo", encoding="utf-8")
            first.write_text("holo", encoding="utf-8")

            with expand_run_inputs(root) as expanded:
                self.assertEqual(
                    (Path("a/first.holo"), Path("b/second.holo")),
                    tuple(path.relative_to(expanded.batch_root) for path in expanded.paths),
                )
                self.assertEqual(root.resolve(), expanded.batch_root)

    def test_zip_extraction_rejects_parent_traversal(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            archive_path = Path(temp_dir) / "unsafe.zip"
            with zipfile.ZipFile(archive_path, "w") as archive:
                archive.writestr("../escaped.txt", "unsafe")

            with self.assertRaisesRegex(ValueError, "Unsafe ZIP"):
                with extracted_zip_tree(archive_path):
                    pass


if __name__ == "__main__":
    unittest.main()
