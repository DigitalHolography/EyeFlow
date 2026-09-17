"""Tests for the hidden moment0 spatial-gradient pipeline."""

from __future__ import annotations

from io import BytesIO
import struct
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np
from PIL import Image

import pipelines
from input_output.output_manager import OutputManager
from pipeline_engine import PIPELINE_REGISTRY, PipelineContext, PipelineDAG
from pipelines.spatial_gradient_moment0.runner import (
    AVI_FILENAME,
    CONTRAST_GAMMA,
    CONTRAST_HIGH_PERCENTILE,
    run_spatial_gradient_moment0,
    spatial_gradient,
)


class SpatialGradientMoment0Tests(unittest.TestCase):
    def test_spatial_gradient_uses_sobel_magnitude(self) -> None:
        ramp = np.tile(np.arange(5, dtype=np.float32), (4, 1))

        gradient = spatial_gradient(ramp)

        np.testing.assert_allclose(gradient[:, 1:-1], 8.0)
        np.testing.assert_allclose(gradient[:, (0, -1)], 4.0)

    def test_temporal_filter_reads_current_shared_window(self) -> None:
        from unittest.mock import patch
        from pipelines.spatial_gradient_moment0 import runner
        frames = [np.full((2, 2), value, np.float32) for value in (0, 10, 0)]
        with patch.object(runner, "TEMPORAL_MEDIAN_WINDOW", 1):
            np.testing.assert_array_equal(list(runner.temporal_median_filter(frames))[1], 10)
        with patch.object(runner, "TEMPORAL_MEDIAN_WINDOW", 3):
            np.testing.assert_array_equal(list(runner.temporal_median_filter(frames))[1], 0)

    def test_all_pre_and_post_filter_combinations_match_reference(self):
        from unittest.mock import patch
        from scipy import ndimage
        from pipelines.spatial_gradient_moment0 import runner
        stack = np.random.default_rng(12).normal(size=(21, 4, 5)).astype(np.float32)
        def reference_filter(values, kind):
            if kind == "median":
                return ndimage.median_filter(values, size=(17, 1, 1), mode="nearest")
            if kind == "gaussian":
                return ndimage.gaussian_filter1d(values, sigma=2, axis=0, mode="nearest")
            return values
        for pre in ("none", "median", "gaussian"):
            for post in ("none", "median", "gaussian"):
                with self.subTest(pre=pre, post=post), \
                     patch.object(runner, "PRE_SOBEL_TEMPORAL_FILTER", pre), \
                     patch.object(runner, "POST_SOBEL_TEMPORAL_FILTER", post):
                    expected = reference_filter(stack, pre)
                    expected = np.stack([runner.spatial_gradient(frame) for frame in expected])
                    expected = reference_filter(expected, post)
                    actual = np.stack(list(runner.filtered_segment_gradients(stack)))
                    np.testing.assert_allclose(actual, expected, rtol=1e-6, atol=1e-6)
                    metadata = runner._temporal_filter_metadata()
                    self.assertEqual(pre, metadata["pre_sobel_temporal_filter"])
                    self.assertEqual(post, metadata["post_sobel_temporal_filter"])

    def test_gaussian_preserves_padding_and_does_not_blur_spatially(self):
        from pipelines.spatial_gradient_moment0 import runner
        stack = np.zeros((3, 3, 4), np.float32)
        stack[:, :, 0] = np.nan
        stack[:, 1, 2] = 10
        result = np.stack(list(runner.apply_temporal_filter(stack, "gaussian")))
        np.testing.assert_allclose(result, stack, equal_nan=True)
        self.assertEqual([], list(runner.apply_temporal_filter([], "gaussian")))

    def test_dispatcher_rejects_unknown_types_and_accepts_new_filters(self):
        from unittest.mock import patch
        from pipelines.spatial_gradient_moment0 import runner
        with self.assertRaisesRegex(ValueError, "Unknown temporal filter"):
            runner.apply_temporal_filter([], "typo")
        with patch.dict(runner.TEMPORAL_FILTERS, {"double": lambda frames: (f * 2 for f in frames)}):
            result = list(runner.apply_temporal_filter([np.ones((2, 2))], "double"))
        np.testing.assert_array_equal(result[0], 2)

    def test_pipeline_is_hidden_and_required_by_waveform_velocity(self) -> None:
        pipelines.load_pipeline_catalog()

        descriptor = PIPELINE_REGISTRY["spatial_gradient_moment0"]
        self.assertEqual("hidden", descriptor.visibility)
        self.assertIn(
            "spatial_gradient_moment0",
            PIPELINE_REGISTRY["waveform_velocity"].dag_requires,
        )
        plan = PipelineDAG(PIPELINE_REGISTRY.values()).resolve_targets(
            ["waveform_velocity"]
        )
        self.assertLess(plan.names.index("waveform_velocity_core"),
                        plan.names.index("spatial_gradient_moment0"))
        self.assertLess(
            plan.names.index("spatial_gradient_moment0"),
            plan.names.index("waveform_velocity"),
        )

    def test_exports_only_one_processed_segment_video(self) -> None:
        from types import SimpleNamespace
        from unittest.mock import patch
        from pipelines.spatial_gradient_moment0 import runner
        from pipelines.waveform_velocity_core.runner import WAVEFORM_CONTEXT_STATE
        frames = np.zeros((19, 5, 5), np.float32)
        frames[:, :, 2:] = 10
        frames[9, 2, 2] = 1000
        segments = SimpleNamespace(
            velocity_profiles=np.zeros((1, 1, 19, 181), np.float32),
            topology=SimpleNamespace(valid_segments=np.ones((1, 1), bool)),
            segment_center_xy=np.array([[[2, 2]]]),
            profile_window_bounds_xyxy=np.array([[[0, 5, 0, 5]]]),
            profile_window_side_pixels=5, profile_rotation_degrees=np.array([[0.0]]),
            segment_masks=np.ones((1, 1, 181, 181), bool),
            labels=np.ones((5, 5)), branch_ids=np.array([1]),
        )
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            with h5py.File(root / "hd.h5", "w") as hd, h5py.File(root / "work.h5", "w") as work:
                hd.create_dataset("moment0ff", data=frames)
                hd.create_dataset("sampling_freq", data=100.0)
                hd.create_dataset("batch_stride", data=2.0)
                output = OutputManager.from_holo(root / "scan.holo", output_root=root)
                output.prepare()
                ctx = PipelineContext(work_h5=work, holodoppler_h5=hd,
                    doppler_vision_h5=None, pipeline_name="spatial_gradient_moment0",
                    output_manager=output)
                ctx.state.set(WAVEFORM_CONTEXT_STATE, SimpleNamespace(
                    artery_segment_result=segments, vein_segment_result=segments))
                with patch.object(runner, "_export_debug_segment", wraps=runner._export_debug_segment) as export:
                    results = run_spatial_gradient_moment0(ctx)
                    self.assertEqual(1, export.call_count)
                    debug_stack = export.call_args.args[1]
                    np.testing.assert_allclose(debug_stack[0], debug_stack[9], equal_nan=True)
                    self.assertEqual(19, len(_jpeg_frames(output.path_for(
                        runner.OutputType.AVI, AVI_FILENAME).read_bytes())))
                self.assertIs(results, ctx.state.get(runner.STATE_KEY))
                self.assertEqual([], list(root.rglob("*.png")))



def _jpeg_frames(contents: bytes) -> list[np.ndarray]:
    position = contents.index(b"movi") + 4
    index_position = contents.index(b"idx1", position)
    frames: list[np.ndarray] = []
    while position < index_position:
        chunk_size = struct.unpack_from("<I", contents, position + 4)[0]
        payload_start = position + 8
        if contents[position : position + 4] == b"00dc":
            with Image.open(BytesIO(contents[payload_start : payload_start + chunk_size])) as image:
                frames.append(np.asarray(image.convert("RGB")))
        position = payload_start + chunk_size + chunk_size % 2
    return frames


if __name__ == "__main__":
    unittest.main()
