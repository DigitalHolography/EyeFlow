"""Tests for branch lumen-size PNGs and branch-based quartile selection."""

from __future__ import annotations

import tempfile
import unittest
import warnings
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import h5py
import numpy as np
from matplotlib.figure import Figure
from PIL import Image

from input_output.output_manager import OutputManager
from pipeline_engine import DatasetValue, PipelineContext
from pipeline_engine.context import apply_pipeline_result
from pipelines.spatial_gradient_moment0.lumen_size import export_lumen_size_pngs
from pipelines.waveform_velocity.spatial_gradient_profiles import (
    pack_spatial_gradient_profile_outputs,
)


class LumenSizePngTests(unittest.TestCase):
    def test_exports_branch_curves_and_mean_of_branches_ranked_by_temporal_median(self):
        # The first branch has a large transient but a small temporal median.
        # The last branch has no data; one selected branch has a missing phase.
        signals = np.tile(np.arange(1.0, 10.0, dtype=np.float32), (4, 1))
        signals[-1, 0] = 1000.0
        signals[:, 6] = [6.0, 7.0, 7.0, 100.0]
        signals[:, 7] = [8.0, np.nan, 8.0, -100.0]
        signals[:, 8] = np.nan
        values = np.empty((4, 2, 9, 3), dtype=np.float32)
        # Unequal valid sample counts make sequential medians differ from
        # the requested joint median over beat and radius.
        values[:, 0] = signals[:, :, None] + [-25.0, -5.0, 5.0]
        values[:, 1] = signals[:, :, None] + [15.0, np.nan, np.nan]
        branch_ids = np.arange(11, 20)
        profiles = np.ones((3, 9, 5, 11), dtype=np.float32)
        segments = SimpleNamespace(
            velocity_profiles=profiles,
            transverse_velocity_profiles_masked=profiles,
        )
        with patch(
            "pipelines.waveform_velocity.spatial_gradient_profiles._spatial_gradient_peak_metrics",
            return_value={
                "Masked/tbkr/lumen/size": DatasetValue(values),
                "Masked/tbkr/lumen/size_qc": DatasetValue(np.zeros_like(values)),
            },
        ):
            outputs = pack_spatial_gradient_profile_outputs(segments, segments, [0, 2, 4])
        snapshots = []
        savefig = Figure.savefig

        def capture(fig, *args, **kwargs):
            ax = fig.axes[0]
            snapshots.append(
                (
                    [(line.get_xdata(), line.get_ydata(), line.get_label()) for line in ax.lines],
                    ax.get_xlabel(),
                    ax.get_ylabel(),
                    ax.get_xlim(),
                )
            )
            return savefig(fig, *args, **kwargs)

        with tempfile.TemporaryDirectory() as temp_dir:
            output = OutputManager.from_holo(
                Path(temp_dir) / "sample.holo", output_root=Path(temp_dir)
            )
            h5_path = Path(temp_dir) / "output.h5"
            with h5py.File(h5_path, "w") as work:
                ctx = PipelineContext(
                    work_h5=work,
                    holodoppler_h5=None,
                    doppler_vision_h5=None,
                )
                apply_pipeline_result(ctx, outputs)
            with h5py.File(h5_path, "r") as work:
                for vessel in ("Artery", "Vein"):
                    source_path = (
                        f"/Processing/SpatialGradientMetrics/{vessel}/Transverse/"
                        "Masked/tbkr/lumen/size"
                    )
                    path = (
                        f"/Processing/SpatialGradientMetrics/{vessel}/Transverse/"
                        "Masked/tk/lumen_size"
                    )
                    dataset = work[path]
                    np.testing.assert_allclose(dataset[()], signals, equal_nan=True)
                    self.assertEqual((4, 9), dataset.shape)
                    self.assertEqual(np.dtype(np.float32), dataset.dtype)
                    self.assertEqual(["time", "branch"], list(dataset.attrs["dimDesc"]))
                    self.assertEqual("pixels", dataset.attrs["unit"])
                    self.assertEqual(source_path, dataset.attrs["source_metric"])
                    self.assertEqual(0, dataset.attrs["qc_applied"])
            with patch.object(Figure, "savefig", new=capture):
                paths = []
                for vessel in ("Artery", "Vein"):
                    path = (
                        f"Processing/SpatialGradientMetrics/{vessel}/Transverse/"
                        "Masked/tk/lumen_size"
                    )
                    paths.extend(
                        export_lumen_size_pngs(
                            output,
                            outputs[path].data,
                            branch_ids,
                            vessel_name=vessel,
                            period_seconds=0.8,
                        )
                    )

            self.assertEqual(
                {
                    "lumen_size_by_branch_artery.png",
                    "lumen_size_by_branch_top_quartile_artery.png",
                    "lumen_size_by_branch_vein.png",
                    "lumen_size_by_branch_top_quartile_vein.png",
                },
                {path.name for path in paths},
            )
            for path in paths:
                self.assertEqual("lumen_size", path.parent.name)
                self.assertEqual("png", path.parent.parent.name)
                with Image.open(path) as image:
                    self.assertEqual("PNG", image.format)
                    self.assertGreater(image.width, 500)
                    self.assertGreater(image.height, 300)
                    image.verify()

        for by_branch, top_quartile in (snapshots[:2], snapshots[2:]):
            curves, xlabel, ylabel, xlim = by_branch
            self.assertEqual(8, len(curves))
            self.assertEqual("Time (s)", xlabel)
            self.assertEqual("Lumen size (pixels)", ylabel)
            self.assertEqual((0.0, 0.8), xlim)
            for index, (time_seconds, signal, label) in enumerate(curves):
                np.testing.assert_allclose(time_seconds, [0.0, 0.2, 0.4, 0.6])
                np.testing.assert_allclose(signal, signals[:, index], equal_nan=True)
                self.assertEqual(f"Branch {branch_ids[index]}", label)
            mean_curves = top_quartile[0]
            self.assertEqual(1, len(mean_curves))
            np.testing.assert_allclose(mean_curves[0][0], [0.0, 0.2, 0.4, 0.6])
            self.assertEqual("Time (s)", top_quartile[1])
            self.assertEqual((0.0, 0.8), top_quartile[3])
            np.testing.assert_allclose(mean_curves[0][1], [7.0, 7.0, 7.5, 0.0])
            self.assertIn("n=2", mean_curves[0][2])

    def test_includes_ties_at_quartile_boundary(self):
        values = np.asarray([1.0, 2.0, 3.0, 4.0, 4.0], dtype=np.float32)
        values = values.reshape(1, 5)
        with patch(
            "pipelines.spatial_gradient_moment0.lumen_size._save_lumen_size_plot"
        ) as save_plot:
            export_lumen_size_pngs(
                None, values, np.arange(5), vessel_name="Artery", period_seconds=0.8
            )
        np.testing.assert_array_equal(save_plot.call_args_list[1].args[1], [[4.0]])
        self.assertIn("n=2", save_plot.call_args_list[1].args[2][0])

    def test_all_nan_and_empty_inputs_export_without_reduction_warnings(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            output = OutputManager.from_holo(
                Path(temp_dir) / "sample.holo", output_root=Path(temp_dir)
            )
            for shape in ((4, 3), (4, 0), (0, 3)):
                with self.subTest(shape=shape), warnings.catch_warnings():
                    warnings.simplefilter("error", category=RuntimeWarning)
                    paths = export_lumen_size_pngs(
                        output,
                        np.full(shape, np.nan, dtype=np.float32),
                        np.arange(shape[1]),
                        vessel_name="Artery",
                        period_seconds=0.8,
                    )
                    self.assertTrue(all(path.is_file() for path in paths))


if __name__ == "__main__":
    unittest.main()
