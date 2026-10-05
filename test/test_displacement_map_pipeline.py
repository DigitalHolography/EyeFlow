"""Tests for displacement-map pipeline inputs and persisted outputs."""

from __future__ import annotations

import tempfile
import unittest
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import h5py
import numpy as np

from input_output.inputs import HoloRunLayout
from input_output.output_manager import OutputManager, OutputType
from pipeline_engine import PipelineContext
from pipelines.displacement_map import registration
from input_output.frame_sequences import FrameSequence
from pipelines.displacement_map.runner import (
    ARTERY_MASK_PATH,
    DISPLACEMENT_MAP_STATE,
    LABELED_VESSELS_PATH,
    MAGNITUDE_VIDEO_FILENAME,
    VEIN_MASK_PATH,
    VESSEL_MASK_PATH,
    DisplacementMapArtifacts,
    DisplacementMapInputs,
    DisplacementMaskInput,
    DisplacementMapPipelineConfig,
    attach_displacement_segment_profiles,
    resolve_moment_dataset,
    resolve_retina_mask,
    resolve_retina_masks,
    run_displacement_map,
)


@dataclass(frozen=True)
class _SegmentProfiles:
    topology: object
    displacements: dict[str, object]


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
        )

        self.assertEqual((16, 16, 2), field.shape)
        self.assertEqual(np.float32, field.dtype)
        self.assertTrue(np.isfinite(field).all())
        np.testing.assert_allclose(field, 0.0, atol=1.0e-6)
        self.assertTrue(np.isfinite(metric))


class DisplacementMapInputTests(unittest.TestCase):
    def test_frame_sequence_reuses_open_h5_dataset(self) -> None:
        with h5py.File("already_open.h5", "w", driver="core", backing_store=False) as hd:
            dataset = hd.create_dataset(
                "moment0",
                data=np.arange(24, dtype=np.float32).reshape(3, 2, 4),
            )
            with patch(
                "input_output.frame_sequences.h5py.File",
                side_effect=AssertionError("HDF5 input was reopened"),
            ):
                sequence = FrameSequence(
                    Path("not-on-disk.h5"), "moment0", 0, 10.0, 1.0, 99.5,
                    h5_source=dataset,
                )
                frames = list(sequence.iter_frames())

        self.assertEqual(3, sequence.frame_count)
        self.assertEqual(3, len(frames))
        self.assertEqual((2, 4, 3), frames[0].shape)

    def test_displacement_pipeline_attaches_segment_results(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            field_path = Path(temp_dir) / "field.npy"
            np.save(field_path, np.zeros((2, 3, 4, 2), dtype=np.float32))
            artifacts = DisplacementMapArtifacts(
                registration_method="method",
                field_paths_by_vessel={
                    "artery": field_path,
                    "vein": field_path,
                },
                temporary_directory=SimpleNamespace(cleanup=lambda: None),
            )
            ctx = SimpleNamespace(
                pipeline_scheduled=lambda name: name == "displacement_map",
                state=SimpleNamespace(get=lambda key: artifacts),
            )
            profiles = {
                name: _SegmentProfiles(
                    topology=SimpleNamespace(prepared_topology=f"{name} topology"),
                    displacements={},
                )
                for name in ("artery", "vein")
            }

            with patch(
                "pipelines.displacement_map.runner.analyze_displacement_segments",
                side_effect=lambda maps, topology, **kwargs: {
                    "method": (next(iter(maps)), topology, kwargs["retain_maps"])
                },
            ):
                attached = attach_displacement_segment_profiles(
                    ctx,
                    profiles,
                    retain_maps=True,
                    profile_settings=SimpleNamespace(working_memory_mb=64.0),
                )

        self.assertEqual(
            ("method", "artery topology", True),
            attached["artery"].displacements["method"],
        )
        self.assertEqual(
            ("method", "vein topology", True),
            attached["vein"].displacements["method"],
        )

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

    def test_both_mode_resolves_separate_artery_and_vein_masks(self) -> None:
        artery = np.asarray([[1, 0], [0, 0]], dtype=np.uint8)
        vein = np.asarray([[0, 0], [0, 1]], dtype=np.uint8)
        with h5py.File("mask.h5", "w", driver="core", backing_store=False) as dv:
            dv.create_dataset(ARTERY_MASK_PATH, data=artery)
            dv.create_dataset(VEIN_MASK_PATH, data=vein)

            masks = resolve_retina_masks(dv, artery.shape, "both")

        self.assertEqual("both", DisplacementMapPipelineConfig().mask_mode)
        self.assertEqual(["artery", "vein"], [item.name for item in masks])
        np.testing.assert_array_equal(masks[0].mask, artery != 0)
        np.testing.assert_array_equal(masks[1].mask, vein != 0)


class DisplacementMapRunnerTests(unittest.TestCase):
    def test_both_masks_keep_two_separate_videos_in_displacement_avi_folder(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir)
            hd_path = root / "scan_HD.h5"
            manager = OutputManager(
                HoloRunLayout.from_holo(root / "scan.holo", output_root=root / "outputs")
            )
            written: list[Path] = []

            def fake_motion_map(config, *, analysis_mask_array, magnitude_video_path, h5_source):
                self.assertEqual("/moment0", h5_source.name)
                self.assertEqual((2, 4), analysis_mask_array.shape)
                magnitude_video_path.write_bytes(b"fake avi")
                written.append(magnitude_video_path)
                field_path = config.output_dir / "displacement_field.npy"
                np.save(field_path, np.zeros((3, 2, 4, 2), dtype=np.float32))
                return {"displacement_field": field_path}

            with (
                h5py.File(hd_path, "w") as hd,
                h5py.File(root / "work.h5", "w") as work_h5,
            ):
                moment = hd.create_dataset("moment0", data=np.ones((3, 2, 4), np.float32))
                masks = tuple(
                    DisplacementMaskInput(
                        name=vessel,
                        vessels=(vessel,),
                        mask=np.ones((2, 4), dtype=bool),
                        source=f"{vessel}_mask",
                    )
                    for vessel in ("artery", "vein")
                )
                ctx = PipelineContext(
                    work_h5=work_h5,
                    holodoppler_h5=hd,
                    doppler_vision_h5=None,
                    output_manager=manager,
                )
                with (
                    patch(
                        "pipelines.displacement_map.runner.load_displacement_map_inputs",
                        return_value=DisplacementMapInputs(moment, masks, 10.0),
                    ),
                    patch(
                        "pipelines.displacement_map.runner.create_retinal_motion_map",
                        side_effect=fake_motion_map,
                    ),
                ):
                    run_displacement_map(ctx)
                artifacts = ctx.state.get(DISPLACEMENT_MAP_STATE)
                self.assertEqual({"artery", "vein"}, set(artifacts.field_paths_by_vessel))
                artifacts.cleanup()

            expected = [
                manager.path_for(
                    OutputType.AVI,
                    f"displacement_maps/{vessel}_{MAGNITUDE_VIDEO_FILENAME}",
                )
                for vessel in ("artery", "vein")
            ]
            self.assertEqual(expected, written)

    def test_missing_optic_disc_prepares_only_artery_field(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir)
            hd_path = root / "scan_HD.h5"
            dv_path = root / "scan_DV.h5"
            output_h5_path = root / "work.h5"
            moment = np.ones((3, 2, 4), dtype=np.float32)
            artery = np.asarray(
                [[1, 0, 1, 0], [0, 1, 0, 1]],
                dtype=np.uint8,
            )
            vein = np.asarray(
                [[0, 1, 0, 1], [1, 0, 1, 0]],
                dtype=np.uint8,
            )
            with h5py.File(hd_path, "w") as hd:
                hd.create_dataset("moment0", data=moment)
                hd.create_dataset("sampling_freq", data=np.float32(100.0))
                hd.create_dataset("batch_stride", data=np.float32(10.0))
            with h5py.File(dv_path, "w") as dv:
                dv.create_dataset(ARTERY_MASK_PATH, data=artery)
                dv.create_dataset(VEIN_MASK_PATH, data=vein)

            manager = OutputManager(
                HoloRunLayout.from_holo(
                    root / "scan.holo",
                    output_root=root / "outputs",
                )
            )
            def fake_motion_map(config, *, analysis_mask_array, magnitude_video_path, h5_source):
                self.assertEqual("moment0", config.h5_dataset)
                self.assertEqual(10.0, config.h5_fps)
                self.assertEqual("/moment0", h5_source.name)
                if np.array_equal(analysis_mask_array, artery.astype(bool)):
                    field_value = 1.0
                elif np.array_equal(analysis_mask_array, vein.astype(bool)):
                    field_value = 2.0
                else:
                    self.fail("Unexpected displacement analysis mask.")
                magnitude_video_path.write_bytes(b"fake avi")
                field_path = config.output_dir / "displacement_field.npy"
                np.save(
                    field_path,
                    np.full((3, 2, 4, 2), field_value, dtype=np.float32),
                )
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
                self.assertNotIn("Processing/Displacement/Map", output_h5)
                artifacts = ctx.state.get(DISPLACEMENT_MAP_STATE)
                self.assertIsInstance(artifacts, DisplacementMapArtifacts)
                try:
                    np.testing.assert_array_equal(
                        np.load(artifacts.field_paths_by_vessel["artery"]),
                        1.0,
                    )
                    self.assertNotIn("vein", artifacts.field_paths_by_vessel)
                finally:
                    artifacts.cleanup()

            video_path = manager.path_for(
                OutputType.AVI,
                f"displacement_maps/{MAGNITUDE_VIDEO_FILENAME}",
            )
            self.assertEqual(manager.layout.ef_dir / "avi" / "displacement_maps", video_path.parent)
            self.assertEqual(b"fake avi", video_path.read_bytes())


if __name__ == "__main__":
    unittest.main()
