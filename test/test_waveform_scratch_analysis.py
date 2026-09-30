from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import h5py
import numpy as np

from pipelines.retinal_velocity.estimation import (
    _bounded_inpaint_result,
    _inpaint_frame_batch,
    _signed_rms_difference,
    estimate_retinal_velocity,
)
from input_output.schema import EyeFlowOutputPaths
from pipelines.waveform_velocity.continuous import pack_continuous_velocity_outputs
from pipelines.retinal_velocity.models import (
    RetinalVelocity,
    RetinalVelocityMaps,
    VesselVelocity,
    VesselVelocitySignals,
)
from pipelines.retinal_velocity.outputs import (
    pack_retinal_velocity_outputs,
)
from pipelines.retinal_velocity.scratch import retinal_velocity_scratch_h5


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

    def test_inpaint_fallback_only_bounds_invalid_frames(self) -> None:
        frames = np.asarray(
            [
                [[1.0, 2.0], [3.0, 4.0]],
                [[10.0, 20.0], [30.0, 40.0]],
            ],
            dtype=np.float32,
        )
        mask = np.asarray([[True, False], [False, False]])
        inpainted = np.asarray(
            [
                [[100.0, 2.0], [3.0, 4.0]],
                [[np.inf, 20.0], [30.0, 40.0]],
            ],
            dtype=np.float32,
        )
        fake_inpaint = SimpleNamespace(
            inpaint_biharmonic=lambda *args, **kwargs: np.moveaxis(
                inpainted,
                0,
                -1,
            )
        )

        actual = _inpaint_frame_batch(frames, mask, fake_inpaint)

        np.testing.assert_array_equal(actual[0], inpainted[0])
        self.assertTrue(np.all(np.isfinite(actual[1])))
        self.assertGreaterEqual(float(actual[1].min()), 20.0)
        self.assertLessEqual(float(actual[1].max()), 40.0)

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
        optic_disc_center = (8.0, 8.0)

        with h5py.File("scratch.h5", "w", driver="core", backing_store=False) as h5:
            result = estimate_retinal_velocity(
                moment0=moment0,
                moment2=moment2,
                artery_mask=artery,
                vein_mask=vein,
                optic_disc_center=optic_disc_center,
                local_background_dist=1,
                scratch_h5=h5,
                retain_velocity_video=False,
            )

            self.assertEqual([], list(h5["waveform"].keys()))
        self.assertIsNone(result.maps.velocity)
        self.assertEqual((16, 16), result.maps.delta_frms_average.shape)
        self.assertEqual((3,), result.artery.frms.shape)

        with h5py.File("scratch.h5", "w", driver="core", backing_store=False) as h5:
            retained = estimate_retinal_velocity(
                moment0=moment0,
                moment2=moment2,
                artery_mask=artery,
                vein_mask=vein,
                optic_disc_center=optic_disc_center,
                local_background_dist=1,
                scratch_h5=h5,
                retain_velocity_video=True,
            )
            dataset = h5["waveform/velocity"]
            expected_velocity = np.asarray(dataset)
            self.assertEqual(["velocity"], list(h5["waveform"].keys()))
            self.assertEqual(
                retained.maps.velocity.name,
                dataset.name,
            )
            self.assertIsNone(dataset.compression)

        velocity_output = np.empty_like(moment0)
        with h5py.File("scratch.h5", "w", driver="core", backing_store=False) as h5:
            buffered = estimate_retinal_velocity(
                moment0=moment0,
                moment2=moment2,
                artery_mask=artery,
                vein_mask=vein,
                optic_disc_center=optic_disc_center,
                local_background_dist=1,
                scratch_h5=h5,
                retain_velocity_video=True,
                velocity_video_output=velocity_output,
            )
            self.assertEqual([], list(h5["waveform"].keys()))

        self.assertIs(buffered.maps.velocity, velocity_output)
        np.testing.assert_array_equal(velocity_output, expected_velocity)
        for field, buffered_value, retained_value in (
            (
                "velocity_average",
                buffered.maps.velocity_average,
                retained.maps.velocity_average,
            ),
            ("frms_average", buffered.maps.frms_average, retained.maps.frms_average),
            (
                "frms_background_average",
                buffered.maps.frms_background_average,
                retained.maps.frms_background_average,
            ),
            (
                "delta_frms_average",
                buffered.maps.delta_frms_average,
                retained.maps.delta_frms_average,
            ),
            ("artery_velocity", buffered.artery.velocity, retained.artery.velocity),
            ("vein_velocity", buffered.vein.velocity, retained.vein.velocity),
        ):
            np.testing.assert_array_equal(
                buffered_value,
                retained_value,
                err_msg=field,
            )

    def test_disabled_vein_keeps_artery_background_and_returns_nan_signals(
        self,
    ) -> None:
        rng = np.random.default_rng(19)
        moment0 = (1.0 + rng.random((3, 16, 16))).astype(np.float32)
        moment2 = (2.0 + rng.random((3, 16, 16))).astype(np.float32)
        artery = np.zeros((16, 16), dtype=bool)
        source_vein = np.zeros_like(artery)
        artery[6, 6] = True
        source_vein[9, 9] = True
        background = artery | source_vein

        with h5py.File("scratch.h5", "w", driver="core", backing_store=False) as h5:
            expected = estimate_retinal_velocity(
                moment0=moment0,
                moment2=moment2,
                artery_mask=artery,
                vein_mask=source_vein,
                optic_disc_center=(8.0, 8.0),
                local_background_dist=1,
                scratch_h5=h5,
                retain_velocity_video=False,
            )
        with h5py.File("scratch.h5", "w", driver="core", backing_store=False) as h5:
            actual = estimate_retinal_velocity(
                moment0=moment0,
                moment2=moment2,
                artery_mask=artery,
                vein_mask=np.zeros_like(source_vein),
                background_mask=background,
                optic_disc_center=(8.0, 8.0),
                local_background_dist=1,
                scratch_h5=h5,
                retain_velocity_video=False,
            )

        np.testing.assert_array_equal(
            actual.artery.velocity,
            expected.artery.velocity,
        )
        np.testing.assert_array_equal(
            actual.maps.velocity_average,
            expected.maps.velocity_average,
        )
        for values in (
            actual.vein.velocity,
            actual.vein.frms,
            actual.vein.frms_background,
            actual.vein.delta_frms,
        ):
            self.assertTrue(np.all(np.isnan(values)))

    def test_velocity_estimator_is_independent_of_frame_chunk_size(self) -> None:
        rng = np.random.default_rng(10)
        shape = (17, 20, 18)
        moment0 = (0.5 + 2.0 * rng.random(shape)).astype(np.float32)
        frequency = (
            np.linspace(2.0e5, 8.0e5, shape[0], dtype=np.float32)[:, None, None]
            * (0.8 + 0.4 * rng.random(shape, dtype=np.float32))
        )
        moment2 = (
            np.mean(moment0, axis=(-1, -2), keepdims=True, dtype=np.float32)
            * np.square(frequency)
        ).astype(np.float32)
        artery = np.zeros(shape[1:], dtype=bool)
        vein = np.zeros_like(artery)
        artery[6:9, 5:8] = True
        vein[12:15, 11:14] = True
        optic_disc_center = (9.0, 10.0)

        results = []
        for chunk_size in (1, 2, 7, 32):
            velocity_output = np.empty(shape, dtype=np.float32)
            with (
                patch(
                    "pipelines.retinal_velocity.estimation."
                    "SCRATCH_FRAME_CHUNK_SIZE",
                    chunk_size,
                ),
                h5py.File(
                    "scratch.h5",
                    "w",
                    driver="core",
                    backing_store=False,
                ) as h5,
            ):
                result = estimate_retinal_velocity(
                    moment0=moment0,
                    moment2=moment2,
                    artery_mask=artery,
                    vein_mask=vein,
                    optic_disc_center=optic_disc_center,
                    local_background_dist=2,
                    scratch_h5=h5,
                    velocity_video_output=velocity_output,
                )
                results.append(_retinal_velocity_data_arrays(result))

        expected = results[0]
        for actual in results[1:]:
            for key in expected:
                np.testing.assert_array_equal(actual[key], expected[key])

    def test_scratch_h5_is_memory_backed(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            output_path = Path(temp_dir) / "output.h5"
            with h5py.File(output_path, "w") as output:
                ctx = SimpleNamespace(runtime=SimpleNamespace(work_h5=output))
                with retinal_velocity_scratch_h5(ctx) as scratch:
                    scratch_path = Path(scratch.filename)
                    scratch.create_dataset("large", data=np.ones((2, 3, 4)))
                    self.assertEqual("core", scratch.driver)
                    self.assertEqual("memory", scratch.attrs["storage"])
                    self.assertFalse(scratch_path.exists())
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
        self.assertEqual(
            "Processing/Maps/VelocityAverage/value",
            schema.analysis.velocity_map_avg,
        )
        self.assertEqual(
            "Processing/Maps/FRMSAverage/value",
            schema.analysis.fRMS_avg,
        )
        self.assertEqual(
            "Processing/Maps/FRMSBackgroundAverage/value",
            schema.analysis.fRMS_bkg_avg,
        )
        self.assertEqual(
            "Processing/Maps/DeltaFRMSAverage/value",
            schema.analysis.delta_fRMS_avg,
        )
        self.assertFalse(hasattr(schema, "topology"))
        self.assertTrue(
            schema.segmentation.artery.branch_label_map.startswith("Segmentation/")
        )
        typed = _typed_velocity()
        shared = pack_retinal_velocity_outputs(typed)
        velocity = pack_continuous_velocity_outputs(typed)
        metrics = {**shared, **velocity}

        self.assertFalse(any(path.startswith("analysis/") for path in metrics))
        self.assertFalse(
            any(value is typed.maps.velocity for value in metrics.values())
        )
        self.assertNotIn(schema.analysis.retinal_artery_velocity_signal, shared)
        self.assertNotIn(schema.analysis.retinal_vein_velocity_signal, shared)
        self.assertEqual(
            {
                schema.analysis.retinal_artery_velocity_signal,
                schema.analysis.retinal_vein_velocity_signal,
                schema.analysis.retinal_artery_velocity_signal_band_limited,
                schema.analysis.retinal_vein_velocity_signal_band_limited,
            },
            set(velocity),
        )

def _typed_velocity() -> RetinalVelocity:
    signal = np.arange(8, dtype=np.float32)
    zeros = np.zeros_like(signal)
    cardiac_cycle = SimpleNamespace(
        systole=SimpleNamespace(
            systole_indexes=np.asarray([1, 5], dtype=np.int32),
            min_peak_distance=1,
            min_peak_height=np.float32(0.0),
        ),
        spectral=SimpleNamespace(
            fundamental_hz=1.0,
            heart_rate_bpm=60.0,
            heart_rate_ste_bpm=0.0,
            period_seconds=1.0,
        ),
    )
    return RetinalVelocity(
        maps=RetinalVelocityMaps(
            velocity=np.ones((8, 4, 4), dtype=np.float32),
            moment0_average=np.ones((4, 4), dtype=np.float32),
            velocity_average=np.ones((4, 4), dtype=np.float32),
            frms_average=np.ones((4, 4), dtype=np.float32),
            frms_background_average=np.ones((4, 4), dtype=np.float32),
            delta_frms_average=np.ones((4, 4), dtype=np.float32),
            section_mask=np.ones((4, 4), dtype=bool),
        ),
        artery=VesselVelocity(
            signals=VesselVelocitySignals(signal, zeros, zeros, zeros),
            velocity_filtered=signal,
        ),
        vein=VesselVelocity(
            signals=VesselVelocitySignals(signal, zeros, zeros, zeros),
            velocity_filtered=signal,
        ),
        vessel_frms_background=zeros,
        cardiac_cycle=cardiac_cycle,
        cardiac_cycle_source="artery",
        dt_seconds=0.1,
    )


def _retinal_velocity_data_arrays(result) -> dict[str, np.ndarray]:
    values = {
        "velocity": result.maps.velocity,
        "moment0_average": result.maps.moment0_average,
        "velocity_average": result.maps.velocity_average,
        "frms_average": result.maps.frms_average,
        "frms_background_average": result.maps.frms_background_average,
        "delta_frms_average": result.maps.delta_frms_average,
        "artery_velocity": result.artery.velocity,
        "vein_velocity": result.vein.velocity,
        "artery_frms": result.artery.frms,
        "vein_frms": result.vein.frms,
        "artery_frms_background": result.artery.frms_background,
        "vein_frms_background": result.vein.frms_background,
        "vessel_frms_background": result.vessel_frms_background,
        "artery_delta_frms": result.artery.delta_frms,
        "vein_delta_frms": result.vein.delta_frms,
    }
    return {name: np.asarray(value).copy() for name, value in values.items()}


if __name__ == "__main__":
    unittest.main()
