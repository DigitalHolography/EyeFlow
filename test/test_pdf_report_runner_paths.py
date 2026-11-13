"""Tests for PDF report runner output paths."""

from __future__ import annotations

import sys
import tempfile
import unittest
import re
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import h5py
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

matplotlib.use("Agg", force=True)

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from input_output.runs.layout import HoloRunLayout  # noqa: E402
from input_output.output_manager import OutputManager, OutputType  # noqa: E402
from input_output.reports.pdf_report import write_pdf_report  # noqa: E402
from input_output.schema import EyeFlowOutputPaths  # noqa: E402
from input_output.writers.h5_output import PipelineH5Output  # noqa: E402
from pipeline_engine.dag import PipelineDAG  # noqa: E402
from pipelines import load_pipeline_catalog  # noqa: E402
from pipelines.pdf_report.runner import run_pdf_report  # noqa: E402
from pipelines.pdf_report.artifacts import export_moments_report_figures  # noqa: E402
from pipelines.pdf_report.data import (  # noqa: E402
    ReportData, VolumeRateReportData, extract_report_data,
)
from pipelines.pdf_report.pages import build_velocity_bvr_pages  # noqa: E402
from pipelines.pdf_report.spectra import (  # noqa: E402
    SpectraGapsReportData,
)


class PdfReportRunnerPathTests(unittest.TestCase):
    def test_report_requires_mask_geometry_but_not_optional_gradient_work(self) -> None:
        available, _ = load_pipeline_catalog()
        dag = PipelineDAG(available)
        report_only = dag.resolve_targets(["pdf_report"]).names
        self.assertIn("blood_volume_rate", report_only)
        self.assertNotIn("low_freq_spectra", report_only)
        self.assertNotIn("displacement_map", report_only)
        self.assertNotIn("displacement_profiles", report_only)
        self.assertNotIn("spatial_gradient_moment0", report_only)
        for target in ("displacement_profiles", "spatial_gradient_moment0"):
            selected = dag.resolve_targets([target]).names
            self.assertIn("pdf_width_report", selected)
            self.assertLess(selected.index(target), selected.index("pdf_width_report"))
            self.assertNotIn("pdf_report", selected)
        self.assertIn("pdf_width_report", report_only)
        self.assertLess(report_only.index("pdf_report"), report_only.index("pdf_width_report"))
        self.assertLess(report_only.index("blood_volume_rate"), report_only.index("pdf_report"))
        with_optional = dag.resolve_targets(
            ["pdf_report", "blood_volume_rate", "low_freq_spectra",
             "displacement_profiles", "spatial_gradient_moment0"]
        ).names
        self.assertLess(with_optional.index("blood_volume_rate"), with_optional.index("pdf_report"))
        self.assertLess(with_optional.index("low_freq_spectra"), with_optional.index("pdf_report"))
        self.assertLess(with_optional.index("displacement_profiles"),
                        with_optional.index("pdf_width_report"))
        self.assertLess(with_optional.index("spatial_gradient_moment0"),
                        with_optional.index("pdf_width_report"))

    def test_uses_layout_stem_and_pdf_report_dir(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            input_root = root / "input"
            output_root = root / "output"
            holo_path = input_root / "scan.holo"
            hd_h5_path = input_root / "scan" / "scan_HD" / "h5" / "scan_HD_output.h5"
            manager = OutputManager(
                HoloRunLayout.from_holo(holo_path, output_root=output_root)
            )
            output_h5_path = manager.path_for(OutputType.H5)
            (manager.layout.ef_dir / "png").mkdir(parents=True)
            ctx = _fake_context(manager, output_h5_path, hd_h5_path)
            warnings = []
            ctx.log_warning = warnings.append

            with patch("pipelines.pdf_report.runner.extract_report_data", return_value=ReportData()) as extract, \
                 patch("pipelines.pdf_report.runner.extract_spectra_gaps_report_data", return_value=None), \
                 patch("pipelines.pdf_report.runner.export_moments_report_figures") as export, \
                 patch("pipelines.pdf_report.runner.build_velocity_bvr_pages", return_value=[]) as build, \
                 patch("pipelines.pdf_report.runner.write_pdf_report") as write:
                result = run_pdf_report(ctx)

            extract.assert_called_once_with(ctx.output.h5)
            export.assert_called_once_with(ctx.output, extract.return_value)
            self.assertEqual(manager.layout.ef_dir / "png", build.call_args.args[2])
            pdf_path, pages = write.call_args.args
            self.assertEqual(
                manager.layout.ef_dir / "pdf" / "scan_velocityBvrReport.pdf",
                pdf_path,
            )
            self.assertEqual([], pages)
            self.assertNotIn("_HD_output_EF", str(pdf_path))
            self.assertEqual(2, len(warnings))
            self.assertIsNone(result)

    def test_extracts_active_waveform_output_schema(self) -> None:
        schema = EyeFlowOutputPaths.active()
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "scan_EF.h5"
            with h5py.File(path, "w") as output_h5:
                output_h5.create_dataset(
                    schema.artery_per_beat.velocity_signal,
                    data=np.full((2, 3), 99.0),
                )
                output_h5.create_dataset(
                    schema.artery_per_beat_filtered.velocity_signal_band_limited,
                    data=np.full((2, 3), 77.0),
                )
                output_h5.create_dataset(
                    schema.artery_per_beat_filtered.velocity_signal,
                    data=np.array([[-4.0, 4.0, -4.0], [4.0, -4.0, 4.0]]),
                )
                output_h5.create_dataset(
                    schema.cardiac_cycle.systolic_cycle_duration_seconds,
                    data=np.asarray([0.5, 0.5]),
                )
                output_h5.create_dataset(
                    schema.cardiac_cycle.spectral_heart_rate_bpm,
                    data=np.asarray(96.795),
                )
                output_h5.create_dataset(
                    f"{schema.waveform_shape_metrics_root}/artery/global/raw/RI",
                    data=np.asarray([0.2, 0.4]),
                )
                output_h5.create_dataset(
                    f"{schema.waveform_shape_metrics_root}/artery/global/bandlimited/RI",
                    data=np.asarray([0.9, 0.9]),
                )

            parameters = _read_report(path, include_bvr=False).parameters

        self.assertEqual(4.0, parameters["Average_Arterial_Velocity"]["value"])
        self.assertEqual(4.0, parameters["Min_Arterial_Velocity"]["value"])
        self.assertAlmostEqual(96.795, parameters["heart_beat"]["value"])
        self.assertAlmostEqual(0.3, parameters["ARI"]["value"])

    def test_bandratio_writes_only_the_shared_spectra_report(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            manager = OutputManager(
                HoloRunLayout.from_holo(root / "scan.holo")
            ).for_artifact_namespace("bandratio")
            ctx = _fake_context(
                manager, root / "scan_EF.h5", root / "scan_HD_output.h5",
            )
            ctx.execution_variant = SimpleNamespace(name="frequency_bands")
            spectra = SpectraGapsReportData(sei=0.4)
            with patch("pipelines.pdf_report.runner.extract_report_data") as extract, \
                 patch("pipelines.pdf_report.runner.extract_spectra_gaps_report_data", return_value=spectra), \
                 patch("pipelines.pdf_report.runner.build_spectra_gaps_pages", return_value=[]) as build, \
                 patch("pipelines.pdf_report.runner.write_pdf_report") as write:
                run_pdf_report(ctx)

            extract.assert_not_called()
            build.assert_called_once_with(spectra, "scan")
            path, pages = write.call_args.args
            self.assertEqual(
                manager.layout.ef_dir / "pdf" / "scan_spectraGapsReport.pdf",
                path,
            )
            self.assertEqual([], pages)
            self.assertNotIn("bandratio", path.parts)

    def test_velocity_page_embeds_maps_and_systole_artifacts(self) -> None:
        paths = {
            (kind, vessel): Path(f"{kind}_{vessel}.png")
            for kind in ("vessel_map", "systole")
            for vessel in ("artery", "vein")
        }
        image = np.zeros((8, 12, 3), dtype=np.uint8)
        def load_registered(candidates):
            candidates = tuple(candidates)
            return image if candidates and candidates[0] in paths.values() else None

        with patch(
            "pipelines.pdf_report.pages.load_first_image", side_effect=load_registered
        ) as load:
            fig = build_velocity_bvr_pages(ReportData(), "scan", None, paths)[0]
        try:
            self.assertEqual(12, load.call_count)
            loaded = [
                call.args[0][0] for call in load.call_args_list
                if call.args[0]
            ]
            self.assertCountEqual(paths.values(), loaded)
            self.assertEqual(4, sum(bool(ax.images) for ax in fig.axes))
            titles = [ax.get_title() for ax in fig.axes]
            self.assertIn("Artery filtered per-beat velocity", titles)
            self.assertIn("Vein filtered per-beat velocity", titles)
        finally:
            plt.close(fig)

    def test_filtered_per_beat_plot_is_finite_nan_safe_and_black(self) -> None:
        values = np.array([[1.0, 2.0, np.nan, 4.0], [2.0, 4.0, 6.0, 8.0]])
        report = ReportData(
            parameters={
                "ARI": {"value": 0.4, "unit": ""},
                "VRI": {"value": 0.6, "unit": ""},
            },
            velocity_per_beat={"artery": values, "vein": values},
            velocity_units={"artery": "mm/s", "vein": "mm/s"},
        )
        fig = build_velocity_bvr_pages(report, "scan", None, {})[0]
        try:
            plot = next(
                ax for ax in fig.axes
                if ax.get_title() == "Artery filtered per-beat velocity"
            )
            mean = next(line for line in plot.lines if line.get_label() == "Mean beat")
            np.testing.assert_allclose(mean.get_ydata(), [1.5, 3.0, 6.0, 6.0])
            self.assertEqual("black", mean.get_color())
            envelope = next(
                collection for collection in plot.collections
                if collection.get_label() == "Beat min/max envelope"
            )
            self.assertAlmostEqual(0.18, envelope.get_alpha())
            self.assertIsNone(plot.get_legend())
            self.assertFalse(any("Mean per-beat RI" in text.get_text() for text in plot.texts))
            vein_plot = next(
                ax for ax in fig.axes
                if ax.get_title() == "Vein filtered per-beat velocity"
            )
            self.assertFalse(
                any("Mean per-beat RI" in text.get_text() for text in vein_plot.texts)
            )
            map_plot = next(
                ax for ax in fig.axes if ax.get_title() == "Artery RMS-frequency map"
            )
            self.assertEqual(map_plot.title.get_fontsize(), plot.title.get_fontsize())
            self.assertEqual(map_plot.title.get_fontweight(), plot.title.get_fontweight())
            self.assertEqual(map_plot.title.get_fontfamily(), plot.title.get_fontfamily())
            systole_plot = next(
                ax for ax in fig.axes if ax.get_title() == "Artery systole indices"
            )
            self.assertGreater(systole_plot.get_position().y0, plot.get_position().y0)
            self.assertAlmostEqual(
                systole_plot.get_position().width,
                plot.get_position().width,
                places=3,
            )
            self.assertAlmostEqual(
                systole_plot.get_position().height,
                plot.get_position().height,
                places=3,
            )
            page_width, page_height = fig.get_size_inches()
            plot_ratio = (
                plot.get_position().width * page_width
                / (plot.get_position().height * page_height)
            )
            self.assertAlmostEqual(1.618, plot_ratio, places=2)
            vein_systole_plot = next(
                ax for ax in fig.axes if ax.get_title() == "Vein systole indices"
            )
            self.assertGreater(vein_systole_plot.get_position().y0, vein_plot.get_position().y0)
            summary = next(
                ax for ax in fig.axes
                if any(text.get_text() == "Summary" for text in ax.texts)
            )
            for filtered_plot in (plot, vein_plot):
                self.assertGreater(
                    filtered_plot.get_position().y0,
                    summary.get_position().y1 + 0.04,
                )
            self.assertEqual(map_plot.title.get_fontsize(), systole_plot.title.get_fontsize())
            self.assertEqual(map_plot.title.get_fontfamily(), systole_plot.title.get_fontfamily())
            page_title = next(text for text in fig.texts if text.get_text() == "Velocity")
            self.assertGreater(page_title.get_fontsize(), plot.title.get_fontsize())
        finally:
            plt.close(fig)

    def test_velocity_summary_pairs_artery_and_vein_metrics(self) -> None:
        parameters = {
            "Min_Arterial_Velocity": {"value": 1.0, "unit": "mm/s"},
            "Average_Arterial_Velocity": {"value": 2.0, "unit": "mm/s"},
            "Max_Arterial_Velocity": {"value": 3.0, "unit": "mm/s"},
            "ARI": {"value": 0.4, "unit": ""},
            "API": {"value": 0.5, "unit": ""},
            "Min_Venous_Velocity": {"value": 6.0, "unit": "mm/s"},
            "Average_Venous_Velocity": {"value": 7.0, "unit": "mm/s"},
            "Max_Venous_Velocity": {"value": 8.0, "unit": "mm/s"},
            "VRI": {"value": 0.9, "unit": ""},
            "VPI": {"value": 1.0, "unit": ""},
            "heart_beat": {"value": 72.5, "unit": "bpm"},
            "not_requested": {"value": 99.0, "unit": ""},
        }
        fig = build_velocity_bvr_pages(ReportData(parameters=parameters), "scan", None, {})[0]
        try:
            summary = next(
                ax for ax in fig.axes
                if any(text.get_text() == "Summary" for text in ax.texts)
            )
            texts = {text.get_text(): text for text in summary.texts}
            self.assertEqual("black", texts["Artery"].get_color())
            self.assertEqual("black", texts["Vein"].get_color())
            velocity_unit = r"$\mathrm{mm.s^{-1}}$"
            self.assertLess(texts[f"1.00 {velocity_unit}"].get_position()[0], 0.5)
            self.assertGreater(texts[f"6.00 {velocity_unit}"].get_position()[0], 0.5)
            self.assertIn("Heart rate  |  72.5 bpm", texts)
            self.assertIsNone(texts["Heart rate  |  72.5 bpm"].get_bbox_patch())
            self.assertNotIn("99.00", texts)
            self.assertEqual(0, len(summary.patches))
        finally:
            plt.close(fig)

    def test_second_page_uses_vessel_columns_and_mask_geometry_artifacts(self) -> None:
        artery = VolumeRateReportData(
            beat_curves=np.array([[1.0, 1.2], [2.0, 2.2]]),
            phase_curve=np.array([1.1, 2.1]), overall_value=1.6,
            unit="mm^3/s",
        )
        vein = VolumeRateReportData(
            beat_curves=np.array([[3.0, 3.2], [4.0, 4.2]]),
            phase_curve=np.array([3.1, 4.1]), overall_value=3.6,
            unit="mm^3/s",
        )
        report = ReportData(
            mask_volume_rates={"artery": artery, "vein": vein},
        )
        image = np.zeros((10, 16, 3), dtype=np.uint8)
        with tempfile.TemporaryDirectory() as temp_dir, patch(
            "pipelines.pdf_report.pages.load_first_image", return_value=image,
        ) as load:
            pages = build_velocity_bvr_pages(
                report, "scan", Path(temp_dir) / "png" / "doppler_moments", {},
            )
        try:
            self.assertEqual(2, len(pages))
            labels = {text.get_text() for text in pages[1].texts}
            self.assertIn("Vessel Geometry and Blood-Volume Rate", labels)
            self.assertIn("Artery", labels)
            self.assertIn("Vein", labels)
            titles = [ax.get_title() for ax in pages[1].axes]
            self.assertEqual(2, titles.count("Velocity histogram"))
            self.assertEqual(2, titles.count("Branch identity"))
            self.assertEqual(2, titles.count("Lumen diameter distribution"))
            self.assertEqual(2, titles.count("Mask-derived geometry"))
            loaded = [call.args[0][0].as_posix() for call in load.call_args_list]
            for vessel in ("artery", "vein"):
                self.assertTrue(any(f"histogramVelocity{vessel}.png" in path for path in loaded))
                self.assertTrue(any(
                    f"{vessel}_13_substack_boxes_on_labels_with_rings.png" in path
                    for path in loaded
                ))
                self.assertTrue(any(
                    f"{vessel}_lumen_diameter_distribution.png" in path
                    for path in loaded
                ))
                self.assertTrue(any(
                    f"{vessel}_blood_volume_rate.png" in path for path in loaded
                ))
            value_text = " ".join(
                text.get_text() for ax in pages[1].axes for text in ax.texts
            )
            self.assertIn(r"$\mathrm{mm^3.s^{-1}}$", value_text)
            self.assertNotIn("Spatial-Gradient Edge", value_text)
        finally:
            for page in pages:
                plt.close(page)

    def test_report_extracts_only_mask_derived_volume_rate(self) -> None:
        schema = EyeFlowOutputPaths.active()
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "scan_EF.h5"
            with h5py.File(path, "w") as output_h5:
                output_h5.create_dataset(
                    schema.blood_volume_rate.artery.dynamic_edges,
                    data=np.full((2, 2, 1, 2), 99.0),
                )
                output_h5.create_dataset(
                    schema.blood_volume_rate.artery.static_edges,
                    data=np.array([
                        [
                            [[1.0, 10.0], [2.0, 20.0]],
                            [[3.0, np.nan], [4.0, 40.0]],
                        ],
                        [
                            [[5.0, 50.0], [6.0, 60.0]],
                            [[7.0, 70.0], [8.0, 80.0]],
                        ],
                    ]),
                )
                output_h5.create_dataset(
                    schema.blood_volume_rate.artery.total_masked_edges,
                    data=np.array([[10.0, 12.0], [14.0, 16.0]]),
                )
            rates = _read_report(path).mask_volume_rates

        self.assertEqual(["artery"], list(rates))
        np.testing.assert_allclose(
            rates["artery"].phase_curve, [11.0, 15.0],
        )

    def test_filtered_velocity_exports_matching_png_and_eps(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            manager = OutputManager(
                HoloRunLayout.from_holo(Path(temp_dir) / "scan.holo")
            ).for_artifact_namespace("doppler_moments")
            report = ReportData(
                velocity_per_beat={
                    "artery": np.array([[1.0, 2.0], [2.0, 3.0]]),
                    "vein": np.array([[3.0, 2.0], [2.0, 1.0]]),
                },
            )
            paths = export_moments_report_figures(manager, report)

            self.assertEqual(4, len(paths))
            self.assertEqual({".png", ".eps"}, {path.suffix for path in paths})
            self.assertTrue(all(path.is_file() for path in paths))

    def test_renders_two_page_pdf_with_optional_products(self) -> None:
        schema = EyeFlowOutputPaths.active()
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            path = root / "scan_EF.h5"
            with h5py.File(path, "w") as output_h5:
                output_h5.create_dataset(
                    schema.artery_per_beat_filtered.velocity_signal,
                    data=np.array([[1., 2., 3., 2.], [2., 3., 4., 3.]]),
                )
                output_h5.create_dataset(
                    f"{schema.waveform_shape_metrics_root}/artery/global/raw/RI",
                    data=np.array([0.75, 0.75]),
                )
                output_h5.create_dataset(
                    f"{schema.spectrograms_masked_metrics_root}/SE",
                    data=np.array([0.1, 0.2, 0.3, 0.4]),
                )
                output_h5.create_dataset(
                    schema.blood_volume_rate.artery.total_masked_edges,
                    data=np.array([[1., 2.], [2., 3.], [3., 4.], [2., 3.]]),
                )
            report = _read_report(path)
            pdf_path = write_pdf_report(
                root / "pdf" / "scan_velocityBvrReport.pdf",
                build_velocity_bvr_pages(report, "scan", None, {}),
            )
            content = pdf_path.read_bytes()
            self.assertTrue(content.startswith(b"%PDF-"))
            self.assertEqual(2, len(re.findall(rb"/Type /Page(?!s)", content)))


def _read_report(path: Path, *, processing_root="ProcessingMoments", include_bvr=True):
    with h5py.File(path, "r") as h5file:
        output = PipelineH5Output(h5file, processing_root=processing_root)
        return extract_report_data(output, include_bvr=include_bvr)


def _fake_context(manager: OutputManager, output_h5_path: Path, hd_h5_path: Path):
    logs: list[str] = []
    return SimpleNamespace(
        require_inputs=lambda *names: None,
        log=logs.append,
        log_warning=logs.append,
        state=SimpleNamespace(get=lambda _key, default=None: default),
        output=SimpleNamespace(
            manager=manager,
            h5=SimpleNamespace(filename=str(output_h5_path)),
        ),
        inputs=SimpleNamespace(
            hd=SimpleNamespace(filename=str(hd_h5_path)),
            dv=SimpleNamespace(filename="dv.h5"),
        ),
    )


if __name__ == "__main__":
    unittest.main()
