"""Tests for shared GUI/CLI run orchestration and transactional outputs."""

from __future__ import annotations

import sys
import tempfile
import unittest
import zipfile
from contextlib import redirect_stderr, redirect_stdout
from io import StringIO
from pathlib import Path
from unittest.mock import patch

import h5py

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from app_settings import AppSettingsStore  # noqa: E402
from input_output.archives import extracted_zip_tree  # noqa: E402
from input_output.output_manager import OutputType  # noqa: E402
from pipeline_engine import PipelineDescriptor, ProcessPipeline  # noqa: E402
from pipeline_engine.run_service import (  # noqa: E402
    execute_run,
    expand_run_inputs,
    resolve_run_spec,
)


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


def _descriptor(*, visibility: str = "visible") -> PipelineDescriptor:
    return PipelineDescriptor(
        name="sample",
        description="sample",
        available=True,
        visibility=visibility,
        pipeline_factory=_NoopPipeline,
    )


class RunServiceTests(unittest.TestCase):
    def test_success_replaces_existing_output_only_after_staged_run(self) -> None:
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
                output_manager.write_json({"new": True})
                return output_manager.path_for(OutputType.H5)

            with patch(
                "pipeline_engine.run_service.run_pipelines_to_output",
                side_effect=fake_run,
            ):
                result = execute_run(spec)

            self.assertTrue(result.succeeded)
            self.assertFalse((final_dir / "old.txt").exists())
            json_path = spec.requests[0].output_manager.path_for(OutputType.JSON)
            self.assertTrue(json_path.is_file())
            self.assertFalse(list(final_dir.parent.glob(".*eyeflow-staging-*")))

    def test_failed_staged_run_preserves_existing_output(self) -> None:
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
                output_manager.write_json({"partial": True})
                raise RuntimeError("analysis failed")

            with patch(
                "pipeline_engine.run_service.run_pipelines_to_output",
                side_effect=fake_failure,
            ):
                result = execute_run(spec)

            self.assertEqual(1, len(result.failures))
            self.assertEqual("old", marker.read_text(encoding="utf-8"))
            self.assertFalse(list(final_dir.parent.glob(".*eyeflow-staging-*")))

    def test_commit_refuses_to_replace_non_directory_output(self) -> None:
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
                output_manager.write_json({"new": True})

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
            store.save({"pipeline_visibility": {"sample": True}})

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
                self.assertEqual((first, second), expanded.paths)
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
