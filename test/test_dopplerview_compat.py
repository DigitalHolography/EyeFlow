"""Compatibility tests for DopplerView exports consumed by EyeFlow."""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from input_output.inputs import load_h5_sidecar_config
from input_output.schema import DopplerViewSource, HolodopplerSource
from pipeline_engine.context import RawH5SourceReader
from pipelines.velocity_analysis.sources import (
    VelocityAnalysisSources,
    _load_moment_pair,
)


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

        self.assertFalse(hasattr(source_data, "velocity_analysis"))
        source = source_data.source
        np.testing.assert_array_equal(source.segmentation.vessels.artery, artery_raw.T)
        np.testing.assert_array_equal(source.segmentation.vessels.vein, vein_raw.T)
        expected_optic_disc_mask = np.zeros_like(artery_raw)
        expected_optic_disc_mask[1:3, 0] = True
        np.testing.assert_array_equal(
            source.segmentation.optic_disc.mask,
            expected_optic_disc_mask.T,
        )
        self.assertEqual(source.segmentation.optic_disc.center, (2.0, 1.0))
        self.assertTrue(source_data.provenance["dv_spatial_axes_swapped_to_match_hd"])
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

        self.assertFalse(hasattr(source_data, "velocity_analysis"))

    def test_missing_optic_disc_uses_centered_frame_relative_fallback(self) -> None:
        artery_raw = np.zeros((4, 2), dtype=bool)
        vein_raw = ~artery_raw
        with self._source_pair() as (hd_source, dv_source):
            self._write_hd(hd_source)
            with h5py.File(dv_source, "w") as dv:
                retina = dv.create_group("segmentation/Retina")
                retina.create_dataset("artery_mask", data=artery_raw)
                retina.create_dataset("vein_mask", data=vein_raw)

            source_data = self._load_sources(hd_source, dv_source)

        disc = source_data.source.segmentation.optic_disc
        vessels = source_data.source.segmentation.vessels
        radius_scale = np.hypot(0.5, 1.5)
        self.assertTrue(disc.is_fallback)
        self.assertEqual((2.0, 1.0), disc.center)
        self.assertAlmostEqual(0.20 * radius_scale, disc.width)
        self.assertAlmostEqual(0.20 * radius_scale, disc.height)
        self.assertEqual((2, 4), disc.mask.shape)
        self.assertFalse(np.any(vessels.vein))
        np.testing.assert_array_equal(
            vessels.velocity_background,
            vein_raw.T,
        )

    def test_nonfinite_optic_disc_mask_uses_the_same_fallback(self) -> None:
        artery_raw = np.zeros((4, 2), dtype=bool)
        vein_raw = ~artery_raw
        with self._source_pair() as (hd_source, dv_source):
            self._write_hd(hd_source)
            with h5py.File(dv_source, "w") as dv:
                retina = dv.create_group("segmentation/Retina")
                retina.create_dataset("artery_mask", data=artery_raw)
                retina.create_dataset("vein_mask", data=vein_raw)
                disc = dv.create_group("segmentation/OpticDisc")
                disc.create_dataset("mask", data=np.full((4, 2), np.nan))
                disc.create_dataset("center", data=np.asarray([1.0, 2.0]))
                disc.create_dataset("width", data=np.float32(10.0))
                disc.create_dataset("height", data=np.float32(12.0))

            source_data = self._load_sources(hd_source, dv_source)

        self.assertTrue(source_data.source.segmentation.optic_disc.is_fallback)
        self.assertEqual((2.0, 1.0), source_data.source.segmentation.optic_disc.center)
        self.assertFalse(np.any(source_data.source.segmentation.vessels.vein))

    def test_hd_pixel_pitch_has_no_legacy_fallback_and_must_be_isotropic(self) -> None:
        with self._source_pair() as (hd_source, _):
            self._write_hd(hd_source)
            with h5py.File(hd_source, "a") as hd:
                del hd["HD_parameters"]
            with h5py.File(hd_source, "r") as hd:
                source = HolodopplerSource(RawH5SourceReader(h5file=hd, label="HD"))
                with self.assertRaisesRegex(KeyError, "HD_parameters"):
                    source.pixel_pitch()

        with self._source_pair() as (hd_source, _):
            self._write_hd(hd_source)
            with h5py.File(hd_source, "a") as hd:
                del hd["HD_parameters"]
                hd.create_dataset(
                    "HD_parameters",
                    data=json.dumps({"pixel_pitch": [20e-6, 21e-6]}),
                )
            with h5py.File(hd_source, "r") as hd:
                source = HolodopplerSource(RawH5SourceReader(h5file=hd, label="HD"))
                pitch = source.pixel_pitch()
                self.assertEqual((20e-6, 21e-6), pitch.xy_m)
                with self.assertRaisesRegex(ValueError, "approximately equal"):
                    _ = pitch.isotropic_m

    def test_hd_timing_is_read_from_parameters_mapping(self) -> None:
        with self._source_pair() as (hd_source, _):
            self._write_hd(hd_source)
            with h5py.File(hd_source, "r") as hd:
                timing = HolodopplerSource(RawH5SourceReader(h5file=hd, label="HD")).timing()

        self.assertEqual((100.0, 10.0), (timing.sampling_freq, timing.batch_stride))

    def test_velocity_analysis_uses_one_coherent_raw_moment_mode(
        self,
    ) -> None:
        with self._source_pair() as (hd_source, dv_source):
            self._write_hd(hd_source)
            with h5py.File(hd_source, "a") as hd:
                hd.create_dataset(
                    "moment0ff",
                    data=np.full((3, 2, 4), 7.0, dtype=np.float32),
                )
                hd.create_dataset(
                    "moment2ff",
                    data=np.full((3, 2, 4), 9.0, dtype=np.float32),
                )
            with h5py.File(dv_source, "w") as dv:
                self._write_segmentation(
                    dv,
                    np.zeros((2, 4), dtype=bool),
                    np.ones((2, 4), dtype=bool),
                )
            with h5py.File(hd_source, "r") as hd:
                source = HolodopplerSource(
                    RawH5SourceReader(h5file=hd, label="HD"),
                )
                moment0, moment2 = _load_moment_pair(source)
                selected_moment0 = np.asarray(moment0)
                selected_moment2 = np.asarray(moment2)

        np.testing.assert_array_equal(
            selected_moment0,
            np.ones((3, 2, 4), dtype=np.float32),
        )
        np.testing.assert_array_equal(
            selected_moment2,
            np.ones((3, 2, 4), dtype=np.float32),
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
                    "VelocityEstimation": {"LocalBackgroundDist": 7},
                },
            )

        ring_settings = source_data.source.segmentation.optic_disc.annulus_geometry((200, 400))
        cross_section = source_data.profile_settings
        radius_scale = np.hypot(99.5, 199.5)
        expected_width = 400 / 25 / radius_scale
        self.assertAlmostEqual(2.0 / radius_scale, ring_settings.inner_radius_frac)
        self.assertAlmostEqual(expected_width, ring_settings.ring_width_frac)
        self.assertAlmostEqual(expected_width, ring_settings.segment_length_frac)
        self.assertAlmostEqual(0.02, cross_section.pixel_size_mm)
        self.assertEqual(
            (20e-6, 20e-6),
            source_data.source.holodoppler.pixel_pitch.xy_m,
        )
        self.assertEqual(512.0, cross_section.working_memory_mb)
        self.assertEqual(0.95, cross_section.submask_size_percentile_kept)
        self.assertEqual(7, source_data.source.doppler_view.local_background_dist)

    @staticmethod
    def _write_hd(path: Path) -> None:
        with h5py.File(path, "w") as hd:
            hd.create_dataset("moment0", data=np.ones((3, 2, 4), dtype=np.float32))
            hd.create_dataset("moment2", data=np.ones((3, 2, 4), dtype=np.float32))
            hd.create_dataset(
                "HD_parameters",
                data=json.dumps(
                    {
                        "pixel_pitch": [20e-6, 20e-6],
                        "sampling_freq": 100.0,
                        "batch_stride": 10.0,
                    }
                ),
            )

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
            data=np.arange(artery_mask.size, dtype=np.int32).reshape(artery_mask.shape),
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
            sources = VelocityAnalysisSources(
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
