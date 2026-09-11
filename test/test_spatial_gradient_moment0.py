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
    PNG_FILENAME,
    run_spatial_gradient_moment0,
    spatial_gradient,
)


class SpatialGradientMoment0Tests(unittest.TestCase):
    def test_spatial_gradient_uses_sobel_magnitude(self) -> None:
        ramp = np.tile(np.arange(5, dtype=np.float32), (4, 1))

        gradient = spatial_gradient(ramp)

        np.testing.assert_allclose(gradient[:, 1:-1], 8.0)
        np.testing.assert_allclose(gradient[:, (0, -1)], 4.0)

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
        self.assertLess(
            plan.names.index("spatial_gradient_moment0"),
            plan.names.index("waveform_velocity"),
        )

    def test_exports_named_avi_and_temporal_mean_png(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            hd_path = root / "input.h5"
            work_path = root / "work.h5"
            frames = np.stack(
                [
                    np.tile(np.arange(6, dtype=np.float32), (4, 1)),
                    np.tile(np.arange(6, dtype=np.float32) * 2.0, (4, 1)),
                ]
            )
            with h5py.File(hd_path, "w") as hd, h5py.File(work_path, "w") as work:
                # A different raw moment0 proves that the flat-field input is used.
                hd.create_dataset("moment0", data=np.zeros_like(frames))
                hd.create_dataset("moment0ff", data=frames)
                hd.create_dataset("sampling_freq", data=100.0)
                hd.create_dataset("batch_stride", data=2.0)
                output = OutputManager.from_holo(root / "scan.holo", output_root=root)
                output.prepare()
                ctx = PipelineContext(
                    work_h5=work,
                    holodoppler_h5=hd,
                    doppler_vision_h5=None,
                    pipeline_name="spatial_gradient_moment0",
                    output_manager=output,
                )

                artifacts = run_spatial_gradient_moment0(ctx)

            self.assertEqual(AVI_FILENAME, artifacts.avi_path.name)
            self.assertEqual(PNG_FILENAME, artifacts.mean_png_path.name)
            self.assertEqual("avi", artifacts.avi_path.parent.name)
            self.assertEqual("png", artifacts.mean_png_path.parent.name)
            self.assertEqual(2, len(_jpeg_frames(artifacts.avi_path.read_bytes())))
            with Image.open(artifacts.mean_png_path) as image:
                mean_png = np.asarray(image)
            gradients = [spatial_gradient(frame) for frame in frames]
            mean_gradient = np.mean(gradients, axis=0)
            self.assertAlmostEqual(
                float(np.percentile(mean_gradient, CONTRAST_HIGH_PERCENTILE)),
                artifacts.display_maximum,
            )
            displayed = [
                np.rint(
                    np.power(
                        np.clip(gradient / artifacts.display_maximum, 0.0, 1.0),
                        CONTRAST_GAMMA,
                    )
                    * 255.0
                ).astype(np.uint8)
                for gradient in gradients
            ]
            expected = np.rint(np.mean(displayed, axis=0)).astype(np.uint8)
            np.testing.assert_array_equal(expected, mean_png)
            gradient_video = np.load(artifacts.gradient_path, mmap_mode="r")
            try:
                np.testing.assert_allclose(gradient_video, gradients)
            finally:
                gradient_video._mmap.close()
                artifacts.cleanup()


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
