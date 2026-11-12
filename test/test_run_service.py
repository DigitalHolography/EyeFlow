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
from unittest.mock import patch

import h5py
import numpy as np

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from app_settings import AppSettingsStore  # noqa: E402
from input_output.archives import extracted_zip_tree  # noqa: E402
from input_output.output_manager import OutputType  # noqa: E402
from input_output.writers.h5 import open_h5 as real_open_h5  # noqa: E402
from pipeline_engine import (  # noqa: E402
    ExecutionVariant,
    PipelineDescriptor,
    PipelineOption,
    ProcessPipeline,
)
from pipeline_engine.run_service import (  # noqa: E402
    execute_run,
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


class _VariantProducer(ProcessPipeline):
    name = "variant_source"
    description = "variant source"
    available = True
    requires = []
    missing_deps = []

    def run(self, ctx):
        ctx.publish_execution_variants(
            [
                ExecutionVariant(
                    name="primary",
                    state={"variant_value": 1},
                    output_namespace="ProcessingMoments",
                    artifact_namespace="primary",
                    provenance={"method": "primary"},
                ),
                ExecutionVariant(
                    name="alternate",
                    state={"variant_value": 2},
                    output_namespace="ProcessingBandRatio",
                    artifact_namespace="alternate",
                    provenance={"method": "alternate"},
                ),
            ]
        )


class _VariantConsumer(ProcessPipeline):
    name = "variant_consumer"
    description = "variant consumer"
    available = True
    requires = []
    missing_deps = []

    def run(self, ctx):
        ctx.output.h5.write("ProcessingMoments/VariantValue/value", ctx.state["variant_value"])
        ctx.output.h5.write(
            "ProcessingMoments/ArtifactNamespace/value",
            ctx.output.manager.artifact_folder,
        )


def _descriptor(*, visibility: str = "visible") -> PipelineDescriptor:
    return PipelineDescriptor(
        name="sample",
        description="sample",
        available=True,
        visibility=visibility,
        pipeline_factory=_NoopPipeline,
    )


class RunServiceTests(unittest.TestCase):
    def test_runtime_opens_each_input_once_during_initialization(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            holo = _write_input(root)
            hd_path = root / "scan" / "scan_HD" / "h5" / "scan_HD_output.h5"
            dv_path = root / "scan" / "scan_DV" / "h5" / "scan_DV.h5"
            with h5py.File(hd_path, "a") as hd:
                hd.create_dataset("registration", data=[1, 2, 3])
            with h5py.File(dv_path, "a") as dv:
                dv.create_dataset("app_versions", data=json.dumps({"DV_version": "test"}))
            spec = resolve_run_spec(
                input_paths=[holo], target_names=["sample"], pipelines=[_descriptor()],
            )
            with (
                patch("pipeline_engine.runtime.open_h5", wraps=real_open_h5) as open_input,
                patch(
                    "input_output.writers.h5_initialization.open_h5",
                    side_effect=AssertionError("Initialization reopened an input"),
                ),
            ):
                result = execute_run(spec)
            self.assertTrue(result.succeeded, result.failures)
            self.assertEqual(2, open_input.call_count)
            with h5py.File(result.outputs[0], "r") as output:
                self.assertFalse(output.attrs["copy_hd_passthrough"])
                self.assertNotIn("Meta/registration", output)

    def test_declared_variant_producer_fans_out_its_dag_dependents(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            holo = _write_input(root)
            producer = PipelineDescriptor(
                name="variant_source",
                description="variant source",
                available=True,
                visibility="hidden",
                dag_produces=("variant_data",),
                produces_execution_variants=True,
                pipeline_factory=_VariantProducer,
            )
            consumer = PipelineDescriptor(
                name="variant_consumer",
                description="variant consumer",
                available=True,
                dag_requires=("variant_data",),
                pipeline_factory=_VariantConsumer,
            )
            spec = resolve_run_spec(
                input_paths=[holo],
                target_names=["variant_consumer"],
                pipelines=[producer, consumer],
            )

            result = execute_run(spec)

            self.assertTrue(result.succeeded, result.failures)
            with h5py.File(result.outputs[0], "r") as output:
                self.assertEqual(1, output["ProcessingMoments/VariantValue/value"][()])
                self.assertEqual(2, output["ProcessingBandRatio/VariantValue/value"][()])
                self.assertEqual(
                    "primary",
                    output["ProcessingMoments/ArtifactNamespace/value"][()].decode(),
                )
                self.assertEqual(
                    "alternate",
                    output["ProcessingBandRatio/ArtifactNamespace/value"][()].decode(),
                )
                self.assertEqual(
                    ["primary", "alternate"],
                    list(output.attrs["execution_variants_completed"]),
                )

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

    def test_band_ratio_calibration_is_validated_and_stored(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            holo = _write_input(Path(temp_dir))
            spec = resolve_run_spec(input_paths=[holo], target_names=["sample"],
                                    pipelines=[_descriptor()], band_ratio_frequency_scale_hz=2.5)
            self.assertEqual(2.5, spec.band_ratio_frequency_scale_hz)
            with self.assertRaisesRegex(ValueError, "band_ratio_frequency_scale_hz"):
                resolve_run_spec(input_paths=[holo], target_names=["sample"],
                                 pipelines=[_descriptor()], band_ratio_frequency_scale_hz=0.0)

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

    def test_duplicate_output_destinations_are_rejected_before_execution(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            holo = _write_input(Path(temp_dir))
            with self.assertRaisesRegex(ValueError, "same EyeFlow output"):
                resolve_run_spec(
                    input_paths=[holo, holo],
                    target_names=["sample"],
                    pipelines=[_descriptor()],
                )

    def test_runtime_records_run_and_velocity_metadata(self) -> None:
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
                self.assertEqual(["sample"], list(output_h5.attrs["pipeline_targets"]))
                self.assertEqual({}, json.loads(output_h5.attrs["pipeline_options"]))
                self.assertEqual(
                    ["doppler_moments", "frequency_bands"],
                    list(output_h5.attrs["velocity_estimation_methods"]),
                )
                self.assertEqual(
                    "physical_velocity",
                    output_h5.attrs["velocity_quantity"],
                )
                self.assertEqual("mm/s", output_h5.attrs["velocity_unit"])

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
                    "band_ratio_frequency_scale_hz": 4.0,
                    "velocity_estimation_method": "bandratio",
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
                    ["frequency_bands"],
                    list(output_h5.attrs["velocity_estimation_methods"]),
                )
                self.assertEqual(
                    4.0,
                    output_h5.attrs["band_ratio_frequency_scale_hz"],
                )

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
