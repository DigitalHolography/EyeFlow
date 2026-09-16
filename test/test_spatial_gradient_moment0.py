"""Tests for the hidden moment0 spatial-gradient pipeline."""

from __future__ import annotations

from io import BytesIO
import struct
import tempfile
import unittest
from unittest.mock import patch
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
    def test_imports_do_not_require_other_pipelines_or_cross_section_generator(self) -> None:
        import subprocess
        import sys
        source = Path(__file__).resolve().parents[1] / "src"
        script = f"""
import sys
sys.path.insert(0, {str(source)!r})
for module in (
    'pipelines.waveform_velocity', 'pipelines.waveform_velocity_core',
    'pipelines.displacement_map',
    'calculations.blood_flow_velocity.cross_section.generate_cross_section_signals',
):
    sys.modules[module] = None
import pipelines.spatial_gradient_moment0.profiles
import pipelines.heartbeat_core.runner
"""
        result = subprocess.run(
            [sys.executable, "-c", script], capture_output=True, text=True,
        )
        self.assertEqual(0, result.returncode, result.stderr)

    def test_exports_profiles_without_waveform_or_displacement(self) -> None:
        from pipelines.heartbeat_core.runner import HEARTBEAT_RESULT_STATE, HeartbeatResult
        from pipelines.spatial_gradient_moment0.runner import spatial_gradient_profile_products
        from calculations.topology import (
            prepare_topologies, run_topology_cache, segment_ring_settings, topology_source_id,
        )
        from utils.logger import Logger
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with h5py.File(root / "hd.h5", "w") as hd, h5py.File(
                root / "dv.h5", "w",
            ) as dv, h5py.File(root / "work.h5", "w") as work:
                frame = np.tile(np.arange(61, dtype=np.float32), (61, 1))
                hd.create_dataset("moment0ff", data=np.stack([frame * (i + 1) for i in range(5)]))
                hd.create_dataset("sampling_freq", data=100.)
                hd.create_dataset("batch_stride", data=2.)
                artery = np.zeros((61, 61), bool)
                artery[27:34, 5:56] = True
                vein = np.zeros_like(artery)
                disc = np.zeros_like(artery)
                disc[27:34, 27:34] = True
                dv.create_dataset("segmentation/Retina/artery_mask", data=artery)
                dv.create_dataset("segmentation/Retina/vein_mask", data=vein)
                dv.create_dataset("segmentation/OpticDisc/mask", data=disc)
                dv.create_dataset("segmentation/OpticDisc/center", data=[30., 30.])
                dv.create_dataset("segmentation/OpticDisc/width", data=7.)
                dv.create_dataset("segmentation/OpticDisc/height", data=7.)
                output = OutputManager.from_holo(root / "scan.holo", output_root=root)
                output.prepare()
                ctx = PipelineContext(
                    work_h5=work, holodoppler_h5=hd, doppler_vision_h5=dv,
                    pipeline_name="spatial_gradient_moment0", output_manager=output,
                )
                ctx.state.set(HEARTBEAT_RESULT_STATE, HeartbeatResult(
                    cycle_boundary_indexes=np.array([0, 2, 4], np.int32), index_base=0,
                ))
                prepared = prepare_topologies(
                    {"artery": artery, "vein": vein}, disc,
                    segment_ring_settings(7., 7., image_shape=(61, 61)),
                    source_id=topology_source_id(hd.filename, dv.filename),
                    cache=run_topology_cache(ctx.state.raw), optic_disc_center=[30., 30.],
                )
                with patch("calculations.topology.workflow.prepare_topology", side_effect=AssertionError(
                    "cached topology must be reused",
                )), patch.object(Logger, "log") as log:
                    artifacts = run_spatial_gradient_moment0(ctx)
                products = spatial_gradient_profile_products(ctx)
                self.assertIs(prepared["artery"], products.artery_segments.prepared_topology)
                self.assertTrue(products.outputs)
                self.assertIn("Processing/SpatialGradientProfiles/Artery/Transverse/Masked/SpatialGradientProfile/value", work)
                self.assertIn("Processing/SpatialGradientMetrics/Artery/Transverse/Masked/tbkr/lumen/size", work)
                self.assertFalse(artifacts.gradient_path.exists())
                self.assertIsNone(ctx.state.get("waveform_velocity_context"))
                messages = [call.args[0] for call in log.call_args_list]
                self.assertTrue(any("Topology cache hit: artery" in message for message in messages))
                self.assertTrue(any("fused_transform=" in message for message in messages))

    def test_spatial_gradient_uses_sobel_magnitude(self) -> None:
        ramp = np.tile(np.arange(5, dtype=np.float32), (4, 1))

        gradient = spatial_gradient(ramp)

        np.testing.assert_allclose(gradient[:, 1:-1], 8.0)
        np.testing.assert_allclose(gradient[:, (0, -1)], 4.0)

    def test_pipeline_is_independent_and_consumed_by_waveform_velocity(self) -> None:
        pipelines.load_pipeline_catalog()

        descriptor = PIPELINE_REGISTRY["spatial_gradient_moment0"]
        self.assertEqual("visible", descriptor.visibility)
        self.assertEqual(("heartbeat",), descriptor.dag_requires)
        self.assertEqual("both", descriptor.input_slot)
        own_plan = PipelineDAG(PIPELINE_REGISTRY.values()).resolve_targets(["spatial_gradient_moment0"])
        self.assertEqual(("heartbeat_core", "spatial_gradient_moment0"), own_plan.names)
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

                artifacts = run_spatial_gradient_moment0(ctx, profiles=False)

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
