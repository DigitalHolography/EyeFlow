"""Tests for displacement-map pipeline inputs and persisted outputs."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import h5py
import numpy as np

from input_output.holo_run_layout import HoloRunLayout
from input_output.output_manager import OutputManager, OutputType
from pipeline_engine import PipelineContext
from pipelines.displacement_map import registration
from pipelines.displacement_map.runner import (
    ARTERY_MASK_PATH,
    DISPLACEMENT_MAP_PATH,
    LABELED_VESSELS_PATH,
    MAGNITUDE_VIDEO_FILENAME,
    VEIN_MASK_PATH,
    VESSEL_MASK_PATH,
    resolve_moment_dataset,
    resolve_retina_mask,
    run_displacement_map,
)


class DisplacementRegistrationTests(unittest.TestCase):
    @unittest.skipIf(
        registration.cv2 is None or registration.sitk is None,
        "OpenCV and SimpleITK are optional pipeline dependencies.",
    )
    def test_identical_frames_produce_a_finite_zero_field(self) -> None:
        y, x = np.mgrid[:16, :16]
        image = np.exp(-((x - 8) ** 2 + (y - 8) ** 2) / 20).astype(np.float32)

        field, metric = registration.estimate_registration_field(
            fixed_full=image,
            moving_full=image,
            mask_soft_full=np.ones_like(image),
            initial_full=None,
            method="symmetric_forces_demons",
            scale=0.5,
            iterations=2,
            field_sigma=1.0,
            update_sigma=0.0,
            metric_radius=4,
            learning_rate=1.0,
            bspline_mesh_size=8,
        )

        self.assertEqual((16, 16, 2), field.shape)
        self.assertEqual(np.float32, field.dtype)
        self.assertTrue(np.isfinite(field).all())
        np.testing.assert_allclose(field, 0.0, atol=1.0e-6)
        self.assertTrue(np.isfinite(metric))


class DisplacementMapInputTests(unittest.TestCase):
    def test_resolves_root_moment0_alias(self) -> None:
        with h5py.File("moment_alias.h5", "w", driver="core", backing_store=False) as hd:
            expected = hd.create_dataset("M0", data=np.ones((3, 2, 4), np.float32))

            resolved = resolve_moment_dataset(hd)

            self.assertEqual(expected.name, resolved.name)

    def test_resolves_an_explicit_alternate_root_moment(self) -> None:
        with h5py.File("moment.h5", "w", driver="core", backing_store=False) as hd:
            expected = hd.create_dataset(
                "moment2",
                data=np.ones((3, 2, 4), np.float32),
            )

            resolved = resolve_moment_dataset(hd, "moment2")

            self.assertEqual(expected.name, resolved.name)

    def test_combined_mask_prefers_vessel_mask_and_aligns_transposed_axes(self) -> None:
        vessel = np.asarray(
            [[1, 0], [0, 1], [0, 0], [1, 0]],
            dtype=np.uint8,
        )
        with h5py.File("mask.h5", "w", driver="core", backing_store=False) as dv:
            dv.create_dataset(VESSEL_MASK_PATH, data=vessel)
            dv.create_dataset(LABELED_VESSELS_PATH, data=np.ones_like(vessel))

            mask, source = resolve_retina_mask(dv, (2, 4))

        self.assertEqual(VESSEL_MASK_PATH, source)
        np.testing.assert_array_equal(mask, vessel.T.astype(bool))

    def test_combined_mask_falls_back_to_labeled_vessels(self) -> None:
        labeled = np.asarray([[0, 3], [7, 0]], dtype=np.int32)
        with h5py.File("mask.h5", "w", driver="core", backing_store=False) as dv:
            dv.create_dataset(LABELED_VESSELS_PATH, data=labeled)

            mask, source = resolve_retina_mask(dv, labeled.shape)

        self.assertEqual(LABELED_VESSELS_PATH, source)
        np.testing.assert_array_equal(mask, labeled != 0)

    def test_combined_mask_falls_back_to_artery_vein_union(self) -> None:
        artery = np.asarray([[1, 0], [0, 0]], dtype=np.uint8)
        vein = np.asarray([[0, 0], [0, 1]], dtype=np.uint8)
        with h5py.File("mask.h5", "w", driver="core", backing_store=False) as dv:
            dv.create_dataset(ARTERY_MASK_PATH, data=artery)
            dv.create_dataset(VEIN_MASK_PATH, data=vein)

            mask, source = resolve_retina_mask(dv, artery.shape)

        self.assertEqual(f"{ARTERY_MASK_PATH}+{VEIN_MASK_PATH}", source)
        np.testing.assert_array_equal(mask, (artery | vein).astype(bool))

    def test_explicit_artery_vein_and_labeled_modes_use_requested_dataset(self) -> None:
        values = {
            ARTERY_MASK_PATH: np.asarray([[1, 0], [0, 0]], dtype=np.uint8),
            VEIN_MASK_PATH: np.asarray([[0, 1], [0, 0]], dtype=np.uint8),
            LABELED_VESSELS_PATH: np.asarray([[0, 0], [5, 0]], dtype=np.int32),
        }
        with h5py.File("mask.h5", "w", driver="core", backing_store=False) as dv:
            for path, value in values.items():
                dv.create_dataset(path, data=value)

            for mode, path in (
                ("artery", ARTERY_MASK_PATH),
                ("vein", VEIN_MASK_PATH),
                ("labeled", LABELED_VESSELS_PATH),
            ):
                with self.subTest(mode=mode):
                    mask, source = resolve_retina_mask(dv, (2, 2), mode)
                    self.assertEqual(path, source)
                    np.testing.assert_array_equal(mask, values[path] != 0)


class DisplacementMapRunnerTests(unittest.TestCase):
    def test_writes_dense_field_and_magnitude_video_to_managed_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir)
            hd_path = root / "scan_HD.h5"
            dv_path = root / "scan_DV.h5"
            output_h5_path = root / "work.h5"
            moment = np.ones((3, 2, 4), dtype=np.float32)
            vessel = np.asarray(
                [[1, 0, 1, 0], [0, 1, 0, 1]],
                dtype=np.uint8,
            )
            with h5py.File(hd_path, "w") as hd:
                hd.create_dataset("moment0", data=moment)
                hd.create_dataset("sampling_freq", data=np.float32(100.0))
                hd.create_dataset("batch_stride", data=np.float32(10.0))
            with h5py.File(dv_path, "w") as dv:
                dv.create_dataset(VESSEL_MASK_PATH, data=vessel)

            manager = OutputManager(
                HoloRunLayout.from_holo(
                    root / "scan.holo",
                    output_root=root / "outputs",
                )
            )
            expected_field = np.arange(48, dtype=np.float32).reshape(3, 2, 4, 2)

            def fake_motion_map(config, *, analysis_mask_array, magnitude_video_path):
                self.assertEqual("moment0", config.h5_dataset)
                self.assertEqual(10.0, config.h5_fps)
                np.testing.assert_array_equal(analysis_mask_array, vessel.astype(bool))
                magnitude_video_path.write_bytes(b"fake mp4")
                field_path = config.output_dir / "displacement_field.npy"
                np.save(field_path, expected_field)
                return {"displacement_field": field_path}

            with (
                h5py.File(hd_path, "r") as hd,
                h5py.File(dv_path, "r") as dv,
                h5py.File(output_h5_path, "w") as output_h5,
                patch(
                    "pipelines.displacement_map.runner.create_retinal_motion_map",
                    side_effect=fake_motion_map,
                ),
            ):
                ctx = PipelineContext(
                    work_h5=output_h5,
                    holodoppler_h5=hd,
                    doppler_vision_h5=dv,
                    output_manager=manager,
                    pipeline_name="displacement_map",
                )
                run_displacement_map(ctx)
                saved = output_h5[DISPLACEMENT_MAP_PATH]
                np.testing.assert_array_equal(saved[()], expected_field)
                self.assertEqual("pixels", saved.attrs["units"])
                self.assertEqual("/moment0", saved.attrs["source_moment_path"])
                self.assertEqual(VESSEL_MASK_PATH, saved.attrs["source_mask_path"])

            video_path = manager.path_for(OutputType.MP4, MAGNITUDE_VIDEO_FILENAME)
            self.assertEqual(manager.layout.ef_dir / "mp4", video_path.parent)
            self.assertEqual(b"fake mp4", video_path.read_bytes())


if __name__ == "__main__":
    unittest.main()
