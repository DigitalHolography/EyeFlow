from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import h5py
import numpy as np
from scipy.ndimage import gaussian_filter

from calculations.dopplerview_analysis.flat_field import (
    corrected_flat_field_chunk,
    fit_flat_field_parameters,
)
from calculations.dopplerview_analysis.vessel_velocity_estimator import (
    VesselVelocityEstimatorStep,
)
from input_output.schema import EyeFlowOutputPaths
from pipelines.waveform_velocity.continuous import pack_continuous_velocity_outputs
from pipelines.waveform_velocity_core.dopplerview.models import DopplerViewStepContext
from pipelines.waveform_velocity_core.dopplerview.outputs import (
    pack_dopplerview_shared_outputs,
)
from pipelines.waveform_velocity_core.scratch import waveform_scratch_h5


class FlatFieldTests(unittest.TestCase):
    def test_chunked_flat_field_matches_dopplerview_formula(self) -> None:
        rng = np.random.default_rng(4)
        source = (10.0 + rng.random((5, 12, 10)) * 30.0).astype(np.float32)
        width = 1.4
        border = 0.15

        params = fit_flat_field_parameters(
            source,
            gaussian_width=width,
            border_amount=border,
            frame_chunk_size=2,
        )
        actual = np.concatenate(
            [
                corrected_flat_field_chunk(source, slice(0, 2), params),
                corrected_flat_field_chunk(source, slice(2, 5), params),
            ]
        )

        source64 = source.astype(np.float64)
        minimum = source64.min()
        maximum = source64.max()
        normalized = (source64 - minimum) / (maximum - minimum)
        ratio = normalized / (
            gaussian_filter(
                normalized,
                sigma=(0, width, width),
                mode="reflect",
                truncate=2.0,
            )
            + 1e-8
        )
        y0, y1 = int(np.ceil(12 * border)), int(np.floor(12 * (1 - border)))
        x0, x1 = int(np.ceil(10 * border)), int(np.floor(10 * (1 - border)))
        scale = (
            normalized[:, y0:y1, x0:x1].sum()
            / ratio[:, y0:y1, x0:x1].sum()
        )
        expected = minimum + (maximum - minimum) * scale * ratio

        np.testing.assert_allclose(actual, expected, rtol=2e-6, atol=2e-5)


class VelocitySignalTests(unittest.TestCase):
    def test_estimator_publishes_unfiltered_raw_signals(self) -> None:
        frames = 48
        moment0 = np.ones((frames, 8, 8), dtype=np.float32)
        moment2 = np.ones_like(moment0)
        artery = np.zeros((8, 8), dtype=bool)
        vein = np.zeros((8, 8), dtype=bool)
        artery[2:4, 2:4] = True
        vein[4:6, 4:6] = True
        modulation = (1.5 + np.sin(np.arange(frames) * 1.7)).astype(np.float32)
        moment2[:, artery] = modulation[:, None]
        cache = {
            "moment0": moment0,
            "moment2": moment2,
            "retinal_artery_mask": artery,
            "retinal_vein_mask": vein,
        }
        ctx = DopplerViewStepContext(
            cache=cache,
            holodoppler_config={"sampling_freq": 100.0, "batch_stride": 1.0},
            dopplerview_config={
                "VelocityEstimation": {"LocalBackgroundDist": 1},
                "PulseAnalysis": {"FilterSignals": True, "LowpassFreqHz": 15.0},
            },
        )

        with patch(
            "calculations.dopplerview_analysis.vessel_velocity_estimator."
            "_run_in_parallel",
            side_effect=lambda func, iterable, **kwargs: np.zeros_like(iterable),
        ):
            VesselVelocityEstimatorStep().run(ctx)

        expected = np.nanmean(
            cache["retinal_vessel_velocity"][:, artery],
            axis=1,
        )
        np.testing.assert_allclose(
            cache["retinal_artery_velocity_signal"],
            expected,
        )


class ScratchAndSchemaTests(unittest.TestCase):
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


if __name__ == "__main__":
    unittest.main()
