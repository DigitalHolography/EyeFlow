from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import h5py
import numpy as np
from calculations.dopplerview_analysis.vessel_velocity_estimator import (
    _bounded_inpaint_result,
    _inpaint_frame_batch,
    _signed_rms_difference,
    run_chunked_velocity_estimator,
)
from input_output.schema import EyeFlowOutputPaths
from pipelines.waveform_velocity.continuous import pack_continuous_velocity_outputs
from pipelines.waveform_velocity_core.dopplerview.outputs import (
    pack_dopplerview_shared_outputs,
)
from pipelines.waveform_velocity_core.scratch import waveform_scratch_h5


class ScratchAndSchemaTests(unittest.TestCase):
    def test_inpaint_result_is_finite_and_bounded_by_each_frame_background(self) -> None:
        source = np.asarray(
            [
                [[1.0, 2.0], [3.0, 4.0]],
                [[10.0, 20.0], [30.0, 40.0]],
            ],
            dtype=np.float32,
        )
        mask = np.asarray([[True, False], [False, False]])
        unstable = np.asarray(
            [
                [[np.inf, -100.0], [3.0, 100.0]],
                [[np.nan, -100.0], [30.0, 100.0]],
            ],
            dtype=np.float32,
        )

        actual = _bounded_inpaint_result(unstable, source, mask)

        self.assertTrue(np.all(np.isfinite(actual)))
        self.assertGreaterEqual(float(actual[0].min()), 2.0)
        self.assertLessEqual(float(actual[0].max()), 4.0)
        self.assertGreaterEqual(float(actual[1].min()), 20.0)
        self.assertLessEqual(float(actual[1].max()), 40.0)

    def test_signed_rms_difference_avoids_float32_square_overflow(self) -> None:
        foreground = np.asarray([3.0e20], dtype=np.float32)
        background = np.asarray([1.0e20], dtype=np.float32)

        with np.errstate(over="raise"):
            actual = _signed_rms_difference(foreground, background)

        self.assertTrue(np.all(np.isfinite(actual)))
        np.testing.assert_allclose(
            actual,
            np.sqrt(8.0) * np.float32(1.0e20),
            rtol=1e-6,
        )

    def test_batched_inpainting_matches_independent_frames(self) -> None:
        from skimage.restoration import inpaint

        rng = np.random.default_rng(8)
        frames = rng.random((4, 12, 10), dtype=np.float32)
        mask = np.zeros((12, 10), dtype=bool)
        mask[4:8, 3:7] = True

        expected = np.stack(
            [inpaint.inpaint_biharmonic(frame, mask) for frame in frames],
            axis=0,
        )
        actual = _inpaint_frame_batch(frames, mask, inpaint)

        np.testing.assert_allclose(actual, expected, rtol=1e-6, atol=1e-6)

    def test_velocity_estimator_uses_summary_only_frequency_intermediates(
        self,
    ) -> None:
        rng = np.random.default_rng(9)
        moment0 = (1.0 + rng.random((3, 16, 16))).astype(np.float32)
        moment2 = (2.0 + rng.random((3, 16, 16))).astype(np.float32)
        artery = np.zeros((16, 16), dtype=bool)
        vein = np.zeros_like(artery)
        artery[6, 6] = True
        vein[9, 9] = True

        with h5py.File("scratch.h5", "w", driver="core", backing_store=False) as h5:
            result = run_chunked_velocity_estimator(
                moment0=moment0,
                moment2=moment2,
                artery_mask=artery,
                vein_mask=vein,
                local_background_dist=1,
                scratch_h5=h5,
                retain_velocity_video=False,
            )

            self.assertEqual([], list(h5["waveform"].keys()))
        self.assertIsNone(result["retinal_vessel_velocity"])
        self.assertIsNone(result["fRMS"])
        self.assertIsNone(result["fRMS_bkg"])
        self.assertIsNone(result["deltafRMS"])
        self.assertEqual((16, 16), result["deltafRMS_avg"].shape)
        self.assertEqual((3,), result["retinal_artery_fRMS_signal"].shape)

        with h5py.File("scratch.h5", "w", driver="core", backing_store=False) as h5:
            retained = run_chunked_velocity_estimator(
                moment0=moment0,
                moment2=moment2,
                artery_mask=artery,
                vein_mask=vein,
                local_background_dist=1,
                scratch_h5=h5,
                retain_velocity_video=True,
            )
            dataset = h5["waveform/velocity"]
            self.assertEqual(["velocity"], list(h5["waveform"].keys()))
            self.assertEqual(
                retained["retinal_vessel_velocity"].name,
                dataset.name,
            )
            self.assertIsNone(dataset.compression)

    def test_scratch_h5_is_removed_after_context(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            output_path = Path(temp_dir) / "output.h5"
            with h5py.File(output_path, "w") as output:
                ctx = SimpleNamespace(runtime=SimpleNamespace(work_h5=output))
                with waveform_scratch_h5(ctx) as scratch:
                    scratch_path = Path(scratch.filename)
                    scratch.create_dataset("large", data=np.ones((2, 3, 4)))
                    self.assertTrue(scratch_path.exists())
                self.assertFalse(scratch_path.exists())

    def test_active_schema_has_no_published_velocity_video_or_analysis_group(self) -> None:
        schema = EyeFlowOutputPaths.active()
        self.assertEqual("eyeflow_v2", schema.name)
        self.assertEqual(
            "Processing/Velocity/global/Artery/BandLimited/value",
            schema.analysis.retinal_artery_velocity_signal_band_limited,
        )
        self.assertEqual(
            "Processing/Velocity/global/Vein/BandLimited/value",
            schema.analysis.retinal_vein_velocity_signal_band_limited,
        )
        self.assertEqual(
            "Processing/Velocity/segments/Artery/Raw/value",
            schema.artery_segments.velocity_signal,
        )
        self.assertEqual(
            "Processing/Velocity/segments/Artery/BandLimited/value",
            schema.artery_segments.velocity_signal_band_limited,
        )
        self.assertEqual(
            "Processing/Velocity/segments/Vein/Raw/value",
            schema.vein_segments.velocity_signal,
        )
        self.assertEqual(
            "Processing/Velocity/segments/Vein/BandLimited/value",
            schema.vein_segments.velocity_signal_band_limited,
        )
        self.assertIsNone(schema.analysis.retinal_velocity_array)
        self.assertIsNone(schema.analysis.velocity_map_avg)
        self.assertTrue(
            schema.analysis.fRMS_avg.startswith("Processing/FrequencyMaps/")
        )
        self.assertFalse(hasattr(schema, "topology"))
        self.assertTrue(
            schema.segmentation.artery.branch_label_map.startswith("Segmentation/")
        )
        analysis = {
            "retinal_artery_velocity_signal": np.arange(8, dtype=np.float32),
            "retinal_vein_velocity_signal": np.arange(8, dtype=np.float32),
            "retinal_artery_velocity_signal_filtered": np.arange(8, dtype=np.float32),
            "retinal_vein_velocity_signal_filtered": np.arange(8, dtype=np.float32),
            "retinal_artery_velocity_signal_filtered_perbeat": np.ones((1, 8)),
            "retinal_vessel_velocity": np.ones((8, 4, 4)),
            "velocity_map_avg": np.ones((4, 4)),
            "fRMS_avg": np.ones((4, 4)),
            "fRMS_bkg_avg": np.ones((4, 4)),
            "beat_indices": np.asarray([1, 5], dtype=np.int32),
            "time_per_beat": np.asarray([0.4], dtype=np.float32),
        }
        shared = pack_dopplerview_shared_outputs(analysis)
        velocity = pack_continuous_velocity_outputs(analysis)
        metrics = {**shared, **velocity}

        self.assertFalse(any(path.startswith("analysis/") for path in metrics))
        self.assertFalse(
            any(value is analysis["retinal_vessel_velocity"] for value in metrics.values())
        )
        self.assertNotIn(schema.analysis.retinal_artery_velocity_signal, shared)
        self.assertNotIn(schema.analysis.retinal_vein_velocity_signal, shared)
        self.assertEqual(
            {
                schema.analysis.retinal_artery_velocity_signal,
                schema.analysis.retinal_vein_velocity_signal,
                schema.analysis.retinal_artery_velocity_signal_filtered,
                schema.analysis.retinal_vein_velocity_signal_filtered,
            },
            set(velocity),
        )

        slim_schema = EyeFlowOutputPaths.active("slim_temp")
        slim_shared = pack_dopplerview_shared_outputs(analysis, slim_schema)
        _, velocity_map_attrs = slim_shared[slim_schema.analysis.velocity_map_avg]
        self.assertEqual("mm/s", velocity_map_attrs["unit"])


if __name__ == "__main__":
    unittest.main()
