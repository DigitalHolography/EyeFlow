"""Tests for the standalone waveform-fixture extraction utility."""

from __future__ import annotations

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "tools" / "extract_waveform_fixture.py"
SPEC = importlib.util.spec_from_file_location("extract_waveform_fixture", SCRIPT_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"Could not load {SCRIPT_PATH}")
extractor = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = extractor
SPEC.loader.exec_module(extractor)


class WaveformFixtureExtractionTests(unittest.TestCase):
    def test_extracts_allowlisted_arrays_and_strips_identifying_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir)
            source_path = root / "patient_identifier_EF.h5"
            output_dir = root / "fixtures"
            artery = np.linspace(10.0, 20.0, 201, dtype=np.float32)
            vein = np.linspace(5.0, 8.0, 201, dtype=np.float32)

            with h5py.File(source_path, "w") as source:
                source.attrs["patient_name"] = "Must Not Be Copied"
                analysis = source.create_group("analysis")
                analysis.create_dataset("retinal_artery_velocity_signal", data=artery)
                analysis.create_dataset("retinal_vein_velocity_signal", data=vein)
                analysis.create_dataset(
                    "beat_indices",
                    data=np.asarray([0, 100, 200], dtype=np.int32),
                )
                analysis.create_dataset(
                    "time_per_beat",
                    data=np.asarray([1.0, 1.0], dtype=np.float32),
                )
                analysis.create_dataset(
                    "velocitysignal_filtered",
                    data=artery - np.float32(0.5),
                )
                analysis.create_dataset(
                    "retinal_velocity_array",
                    data=np.ones((16, 16, 16), dtype=np.float32),
                )
                source.create_dataset(
                    "patient_name",
                    data=np.bytes_("Must Not Be Copied"),
                )

            result = extractor.extract_fixture(source_path, output_dir)

            self.assertTrue(result.created)
            self.assertNotIn("patient_identifier", result.output_path.name)
            with h5py.File(result.output_path, "r") as fixture:
                np.testing.assert_array_equal(
                    fixture["waveforms/artery/raw"][()],
                    artery,
                )
                np.testing.assert_array_equal(
                    fixture["waveforms/vein/raw"][()],
                    vein,
                )
                self.assertAlmostEqual(float(fixture["beats/dt_seconds"][()]), 0.01)
                self.assertNotIn("retinal_velocity_array", fixture)
                self.assertNotIn("patient_name", fixture)
                self.assertNotIn("patient_name", fixture.attrs)
                self.assertTrue(bool(fixture.attrs["metadata_stripped"]))
                self.assertTrue(bool(fixture.attrs["contains_patient_derived_data"]))

    def test_supports_slim_output_schema_aliases(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir)
            source_path = root / "slim.h5"
            output_dir = root / "fixtures"

            with h5py.File(source_path, "w") as source:
                source.create_dataset(
                    "artery/velocity/signal/value",
                    data=np.arange(12, dtype=np.float32),
                )
                source.create_dataset(
                    "vein/velocity/signal/value",
                    data=np.arange(12, dtype=np.float32) + 20,
                )
                source.create_dataset(
                    "perbeat/beat_indices/value",
                    data=np.asarray([0, 5, 10], dtype=np.int32),
                )
                source.create_dataset(
                    "perbeat/time_per_beat/value",
                    data=np.asarray([0.5, 0.5], dtype=np.float32),
                )
                source.create_dataset(
                    "artery/velocity/perbeat/signal/value",
                    data=np.ones((2, 8), dtype=np.float32),
                )

            result = extractor.extract_fixture(source_path, output_dir)

            with h5py.File(result.output_path, "r") as fixture:
                self.assertIn("per_beat/artery/raw", fixture)
                self.assertEqual(
                    fixture["beats/boundary_indices"].attrs["index_base"],
                    0,
                )
                self.assertAlmostEqual(float(fixture["beats/dt_seconds"][()]), 0.1)

    def test_reuses_identical_fixture_without_overwriting(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir)
            source_path = root / "source.h5"
            output_dir = root / "fixtures"
            _write_minimal_source(source_path)

            first = extractor.extract_fixture(source_path, output_dir)
            second = extractor.extract_fixture(source_path, output_dir)

            self.assertTrue(first.created)
            self.assertFalse(second.created)
            self.assertEqual(first.output_path, second.output_path)

    def test_rejects_missing_required_waveform(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            source_path = Path(tmp_dir) / "invalid.h5"
            with h5py.File(source_path, "w") as source:
                source.create_dataset(
                    "analysis/retinal_artery_velocity_signal",
                    data=np.arange(10, dtype=np.float32),
                )

            with self.assertRaisesRegex(
                extractor.FixtureExtractionError,
                "Missing required EyeFlow waveform datasets",
            ):
                extractor.prepare_fixture(source_path, max_output_mb=64)

    def test_recursive_discovery_excludes_output_directory(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir)
            nested = root / "inputs" / "nested"
            nested.mkdir(parents=True)
            source_path = nested / "source.h5"
            _write_minimal_source(source_path)
            output_dir = root / "inputs" / "fixtures"
            output_dir.mkdir()
            _write_minimal_source(output_dir / "waveform_fixture_old.h5")

            discovered = extractor.discover_h5_files(
                [root / "inputs"],
                recursive=True,
                output_dir=output_dir,
            )

            self.assertEqual(discovered, [source_path.resolve()])

    def test_discovery_allows_source_file_in_output_directory(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            output_dir = Path(tmp_dir)
            source_path = output_dir / "acquisition.h5"
            generated_path = output_dir / "waveform_fixture_existing.h5"
            _write_minimal_source(source_path)
            _write_minimal_source(generated_path)

            discovered = extractor.discover_h5_files(
                [output_dir],
                recursive=False,
                output_dir=output_dir,
            )

            self.assertEqual(discovered, [source_path.resolve()])


def _write_minimal_source(path: Path) -> None:
    with h5py.File(path, "w") as source:
        source.create_dataset(
            "analysis/retinal_artery_velocity_signal",
            data=np.arange(12, dtype=np.float32),
        )
        source.create_dataset(
            "analysis/retinal_vein_velocity_signal",
            data=np.arange(12, dtype=np.float32) + 20,
        )
        source.create_dataset(
            "analysis/beat_indices",
            data=np.asarray([0, 5, 10], dtype=np.int32),
        )
        source.create_dataset(
            "analysis/time_per_beat",
            data=np.asarray([0.5, 0.5], dtype=np.float32),
        )


if __name__ == "__main__":
    unittest.main()
