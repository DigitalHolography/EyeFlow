"""Compatibility tests for DopplerView exports consumed by EyeFlow."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from input_output import load_h5_sidecar_config  # noqa: E402
from input_output.schema import DopplerViewSource, HolodopplerSource  # noqa: E402
from pipeline_engine.context import RawH5SourceReader  # noqa: E402
from pipelines.waveform_velocity_core.runner import (  # noqa: E402
    _segment_ring_settings,
)
from pipelines.waveform_velocity_core.sources import WaveformVelocitySources  # noqa: E402


class DopplerViewCompatibilityTests(unittest.TestCase):
    def test_dopplerview_json_sidecar_is_loaded(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            h5_path = Path(tmp_dir) / "sample_DV" / "h5" / "sample_DV.h5"
            config_path = h5_path.parent.parent / "json" / "DV_params.json"
            h5_path.parent.mkdir(parents=True)
            config_path.parent.mkdir()
            config_path.write_text(
                '{"Velocity Estimation": {"Local Background Dist": 5}}',
                encoding="utf-8",
            )
            with h5py.File(h5_path, "w") as h5:
                config = load_h5_sidecar_config(h5, source="dv")

            self.assertEqual(5, config["VelocityEstimation"]["LocalBackgroundDist"])

    def test_missing_analysis_is_allowed_and_spatial_masks_align_to_hd(self) -> None:
        artery_raw = np.array(
            [[1, 0], [0, 0], [1, 1], [0, 1]],
            dtype=bool,
        )
        vein_raw = ~artery_raw

        with self._source_pair() as (hd_source, dv_source):
            self._write_hd(hd_source)
            with h5py.File(dv_source, "w") as dv:
                self._write_segmentation(dv, artery_raw, vein_raw)

            source_data = self._load_sources(hd_source, dv_source)

        self.assertIsNone(source_data.dopplerview_analysis)
        np.testing.assert_array_equal(source_data.retinal_artery_mask, artery_raw.T)
        np.testing.assert_array_equal(source_data.retinal_vein_mask, vein_raw.T)
        expected_optic_disc_mask = np.zeros_like(artery_raw)
        expected_optic_disc_mask[1:3, 0] = True
        np.testing.assert_array_equal(
            source_data.optic_disc_mask,
            expected_optic_disc_mask.T,
        )
        np.testing.assert_array_equal(source_data.optic_disc_center, [2.0, 1.0])
        self.assertTrue(
            source_data.provenance["dv_spatial_axes_swapped_to_match_hd"]
        )
        self.assertFalse(source_data.provenance["dopplerview_analysis_available"])
        self.assertTrue(source_data.provenance["has_optic_disc_mask"])

    def test_existing_dopplerview_analysis_is_ignored(self) -> None:
        artery_raw = np.array(
            [[1, 0], [0, 0], [1, 1], [0, 1]],
            dtype=bool,
        )
        vein_raw = ~artery_raw
        velocity_raw = np.arange(24, dtype=np.float32).reshape(3, 4, 2)

        with self._source_pair() as (hd_source, dv_source):
            self._write_hd(hd_source)
            with h5py.File(dv_source, "w") as dv:
                self._write_segmentation(dv, artery_raw, vein_raw)
                self._write_analysis(dv, velocity_raw)

            source_data = self._load_sources(hd_source, dv_source)

        self.assertIsNone(source_data.dopplerview_analysis)
        self.assertFalse(source_data.provenance["dopplerview_analysis_available"])

    def test_precomputed_holodoppler_flat_field_moments_are_reused(
        self,
    ) -> None:
        with self._source_pair() as (hd_source, dv_source):
            self._write_hd(hd_source)
            with h5py.File(hd_source, "a") as hd:
                hd.create_dataset(
                    "moment0ff",
                    data=np.full((3, 2, 4), 7.0, dtype=np.float32),
                )
            with h5py.File(dv_source, "w") as dv:
                self._write_segmentation(
                    dv,
                    np.zeros((2, 4), dtype=bool),
                    np.ones((2, 4), dtype=bool),
                )
            source_data = self._load_sources(hd_source, dv_source)
            with h5py.File(hd_source, "r") as hd:
                source = HolodopplerSource(
                    RawH5SourceReader(h5file=hd, label="HD"),
                )
                selected_moment0 = np.asarray(source.moment0_flat_field_dataset())
                selected_moment2 = source.moment2_flat_field_dataset()

        np.testing.assert_array_equal(
            selected_moment0,
            np.full((3, 2, 4), 7.0, dtype=np.float32),
        )
        self.assertIsNone(selected_moment2)
        self.assertEqual(
            "holodoppler_precomputed_flat_field",
            source_data.moment0_flat_field_source,
        )
        self.assertEqual(
            "holodoppler_raw_moment2",
            source_data.moment2_flat_field_source,
        )

    def test_eyeflow_analysis_settings_are_not_read_from_source_configs(self) -> None:
        artery = np.zeros((2, 4), dtype=bool)
        vein = ~artery
        with self._source_pair() as (hd_source, dv_source):
            self._write_hd(hd_source)
            with h5py.File(dv_source, "w") as dv:
                self._write_segmentation(dv, artery, vein)

            source_data = self._load_sources(
                hd_source,
                dv_source,
                hd_config={
                    "SizeOfField": {"SmallRadiusRatio": 0.45},
                    "generateCrossSectionSignals": {
                        "NumberOfCircles": 3,
                        "SegmentsLength": 0.40,
                        "ScaleFactorWidth": 99.0,
                        "HydrodynamicDiameters": False,
                        "velocityProfileThreshold": 0.99,
                        "RotateFromMask": True,
                        "RefPapillaSize": 99.0,
                        "DefaultPixelSize": 99.0,
                    },
                    "Preprocess": {"InterpolationFactor": 7.0},
                },
                dv_config={
                    "PeripapillaryRingAnalysis": {
                        "RingsNumber": 4,
                        "RingsWidth": 0.30,
                    },
                    "PeripapillaryVascularZone": {
                        "InnerRadius": 0.35,
                        "OuterRadius": 0.45,
                    },
                    "FlatFieldCorrection": {
                        "GWRatio": 0.11,
                        "Border": 0.22,
                    },
                    "VelocityEstimation": {"LocalBackgroundDist": 7},
                },
            )

        ring_settings = _segment_ring_settings()
        cross_section = source_data.cross_section_settings
        self.assertEqual(10, ring_settings.ring_count)
        self.assertEqual(0.10, ring_settings.inner_radius_frac)
        self.assertEqual(0.025, ring_settings.segment_length_frac)
        self.assertEqual(3.0, cross_section.scale_factor_width)
        self.assertTrue(cross_section.hydrodynamic_diameters)
        self.assertEqual(0.5, cross_section.velocity_profile_threshold)
        self.assertFalse(cross_section.rotate_from_mask)
        self.assertAlmostEqual(1.91 / 3.5, cross_section.pixel_size_mm)
        self.assertEqual(7, source_data.local_background_dist)
        self.assertEqual(0.11, source_data.flat_field_gaussian_ratio)
        self.assertEqual(0.22, source_data.flat_field_border)

    @staticmethod
    def _write_hd(path: Path) -> None:
        with h5py.File(path, "w") as hd:
            hd.create_dataset("moment0", data=np.ones((3, 2, 4), dtype=np.float32))
            hd.create_dataset("moment2", data=np.ones((3, 2, 4), dtype=np.float32))
            hd.create_dataset("sampling_freq", data=np.float32(100.0))
            hd.create_dataset("batch_stride", data=np.float32(10.0))

    @staticmethod
    def _write_segmentation(
        h5: h5py.File,
        artery_mask: np.ndarray,
        vein_mask: np.ndarray,
    ) -> None:
        retina = h5.create_group("segmentation/Retina")
        retina.create_dataset("artery_mask", data=artery_mask)
        retina.create_dataset("vein_mask", data=vein_mask)
        retina.create_dataset(
            "labeled_vessels",
            data=np.arange(8, dtype=np.int32).reshape(4, 2),
        )
        optic_disc = h5.create_group("segmentation/OpticDisc")
        optic_disc_mask = np.zeros_like(artery_mask)
        optic_disc_mask[1:3, 0] = True
        optic_disc.create_dataset("mask", data=optic_disc_mask)
        optic_disc.create_dataset("center", data=np.asarray([1.0, 2.0], dtype=np.float32))
        optic_disc.create_dataset("width", data=np.float32(3.0))
        optic_disc.create_dataset("height", data=np.float32(4.0))

    @staticmethod
    def _write_analysis(h5: h5py.File, velocity_raw: np.ndarray) -> None:
        analysis = h5.create_group("analysis")
        analysis.create_dataset("retinal_velocity_array", data=velocity_raw)
        analysis.create_dataset(
            "retinal_artery_velocity_signal",
            data=np.asarray([1000.0, 2000.0, 3000.0], dtype=np.float32),
        )
        analysis.create_dataset(
            "retinal_vein_velocity_signal",
            data=np.asarray([4000.0, 5000.0, 6000.0], dtype=np.float32),
        )
        analysis.create_dataset("velocity_map_avg", data=np.mean(velocity_raw, axis=0))
        analysis.create_dataset("fRMS_avg", data=np.mean(velocity_raw, axis=0))
        analysis.create_dataset("fRMS_bkg_avg", data=np.mean(velocity_raw, axis=0))
        analysis.create_dataset(
            "velocitysignal_per_beat",
            data=np.ones((2, 3), dtype=np.float32),
        )
        analysis.create_dataset(
            "velocitysignal_filtered",
            data=np.asarray([1000.0, 2000.0, 3000.0], dtype=np.float32),
        )
        analysis.create_dataset("beat_indices", data=np.asarray([0, 1, 2], dtype=np.int32))
        analysis.create_dataset("time_per_beat", data=np.asarray([0.1, 0.1], dtype=np.float32))

    @staticmethod
    def _load_sources(
        hd_path: Path,
        dv_path: Path,
        *,
        hd_config: dict[str, object] | None = None,
        dv_config: dict[str, object] | None = None,
    ):
        hd_file = h5py.File(hd_path, "r")
        dv_file = h5py.File(dv_path, "r")
        try:
            sources = WaveformVelocitySources(
                hd=HolodopplerSource(
                    RawH5SourceReader(h5file=hd_file, label="HD"),
                    hd_config,
                ),
                dv=DopplerViewSource(
                    RawH5SourceReader(h5file=dv_file, label="DV"),
                    dv_config,
                ),
            )
            return sources.load()
        finally:
            hd_file.close()
            dv_file.close()

    @staticmethod
    def _source_pair():
        class _SourcePair:
            def __enter__(self):
                self.tmp = tempfile.TemporaryDirectory()
                root = Path(self.tmp.name)
                return root / "sample_HD.h5", root / "sample_DV.h5"

            def __exit__(self, exc_type, exc, tb):
                self.tmp.cleanup()

        return _SourcePair()


if __name__ == "__main__":
    unittest.main()
