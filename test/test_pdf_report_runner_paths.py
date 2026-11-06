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
from pipelines.pdf_report.data import (  # noqa: E402
    ReportData, SpectralReportData, VolumeRateReportData, extract_report_data,
)
from pipelines.pdf_report.pages import build_report_pages  # noqa: E402
from pipelines.pdf_report.width_data import PhaseSummary, WidthReportData  # noqa: E402
from pipelines.pdf_report.width_pages import build_width_page  # noqa: E402


class PdfReportRunnerPathTests(unittest.TestCase):
    def test_optional_pipelines_precede_report_without_becoming_dependencies(self) -> None:
        available, _ = load_pipeline_catalog()
        dag = PipelineDAG(available)
        report_only = dag.resolve_targets(["pdf_report"]).names
        self.assertNotIn("blood_volume_rate", report_only)
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

    def test_report_inherits_velocity_unit_with_dataset_override(self) -> None:
        from pipelines.pdf_report.data import _velocity_dataset_unit

        with h5py.File("report_units.h5", "w", driver="core", backing_store=False) as f:
            f.attrs["velocity_unit"] = "wrong root fallback"
            output = PipelineH5Output(f, processing_root="ProcessingBandRatio",
                                      provenance={"velocity_unit": "mm/s"})
            output.write("ProcessingMoments/signal", [1])
            dataset = output.get("ProcessingMoments/signal")
            self.assertEqual("mm/s", _velocity_dataset_unit(f, dataset))
            dataset.attrs["unit"] = "local override"
            self.assertEqual("local override", _velocity_dataset_unit(f, dataset))

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
                 patch("pipelines.pdf_report.runner.build_report_pages", return_value=[]) as build, \
                 patch("pipelines.pdf_report.runner.write_pdf_report") as write:
                result = run_pdf_report(ctx)

            extract.assert_called_once_with(ctx.output.h5)
            self.assertEqual(manager.layout.ef_dir / "png", build.call_args.args[2])
            pdf_path, pages = write.call_args.args
            self.assertEqual(manager.layout.ef_dir / "pdf" / "scan_report.pdf", pdf_path)
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

            parameters = _read_report(path, include_optional=False).parameters

        self.assertEqual(4.0, parameters["Average_Arterial_Velocity"]["value"])
        self.assertEqual(4.0, parameters["Min_Arterial_Velocity"]["value"])
        self.assertAlmostEqual(96.795, parameters["heart_beat"]["value"])
        self.assertAlmostEqual(0.3, parameters["ARI"]["value"])

    def test_missing_filtered_velocity_does_not_fall_back_to_unfiltered(self) -> None:
        schema = EyeFlowOutputPaths.active()
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "scan_EF.h5"
            with h5py.File(path, "w") as output_h5:
                output_h5.create_dataset(
                    schema.artery_per_beat.velocity_signal,
                    data=np.full((2, 3), 99.0),
                )
                output_h5.create_dataset(
                    schema.cardiac_cycle.systolic_cycle_duration_seconds,
                    data=np.array([0.5, 0.5]),
                )
                output_h5.create_dataset(
                    f"{schema.waveform_shape_metrics_root}/artery/global/raw/RI",
                    data=np.array([0.5, 0.5]),
                )
            parameters = _read_report(path, include_optional=False).parameters
        self.assertNotIn("Average_Arterial_Velocity", parameters)
        self.assertNotIn("heart_beat", parameters)
        self.assertNotIn("ARI", parameters)

    def test_optional_page_uses_shared_spectra_and_variant_volume_rate(self) -> None:
        schema = EyeFlowOutputPaths.active()
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            path = root / "scan_EF.h5"
            with h5py.File(path, "w") as output_h5:
                for processing_root, velocity, flow in (
                    ("ProcessingMoments", 2.0, 4.0),
                    ("ProcessingBandRatio", 5.0, 7.0),
                ):
                    output_h5.create_dataset(
                        schema.artery_per_beat_filtered.velocity_signal.replace(
                            "ProcessingMoments", processing_root
                        ),
                        data=np.full((2, 4), velocity),
                    )
                    output_h5.create_dataset(
                        schema.blood_volume_rate.artery.total_masked_edges.replace(
                            "ProcessingMoments", processing_root
                        ),
                        data=np.full((4, 2), flow),
                    )
                output_h5.create_dataset(
                    f"{schema.spectrograms_masked_metrics_root}/SE",
                    data=np.array([0.1, 0.3, 0.5, 0.7]),
                )
                output_h5.create_dataset(
                    f"{schema.spectrograms_masked_metrics_root}/SEI", data=0.4,
                )
                output_h5.create_dataset(
                    f"{schema.spectrograms_masked_qc.data_completeness}/AcquisitionValidFraction",
                    data=0.8,
                )

            moments = _read_report(path, processing_root="ProcessingMoments")
            bands = _read_report(path, processing_root="ProcessingBandRatio")
            self.assertEqual(2.0, moments.parameters["Average_Arterial_Velocity"]["value"])
            self.assertEqual(5.0, bands.parameters["Average_Arterial_Velocity"]["value"])
            self.assertEqual(4.0, moments.volume_rates["artery"].overall_value)
            self.assertEqual(7.0, bands.volume_rates["artery"].overall_value)
            self.assertEqual(0.4, moments.spectral.sei)
            self.assertEqual(0.4, bands.spectral.sei)

            with patch("input_output.reports.pdf_report.PdfPages") as pages:
                pages.return_value.__enter__.return_value = pages.return_value
                write_pdf_report(root / "pdf" / "scan_report.pdf", build_report_pages(moments, "scan", None, {}))
                self.assertEqual(2, pages.return_value.savefig.call_count)

    def test_report_without_optional_products_has_one_page(self) -> None:
        schema = EyeFlowOutputPaths.active()
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            path = root / "scan_EF.h5"
            with h5py.File(path, "w") as output_h5:
                output_h5.create_dataset(
                    schema.artery_per_beat_filtered.velocity_signal,
                    data=np.full((2, 4), 2.0),
                )
            report = _read_report(path)
            report_pages = build_report_pages(report, "scan", None, {})
            page_labels = " ".join(
                [axis.get_title() for axis in report_pages[0].axes]
                + [text.get_text() for axis in report_pages[0].axes for text in axis.texts]
            )
            self.assertIn("Artery filtered per-beat velocity", page_labels)
            self.assertIn("Vein filtered per-beat velocity", page_labels)
            self.assertIn("Vein systole indices unavailable", page_labels)
            with patch("input_output.reports.pdf_report.PdfPages") as pages:
                pages.return_value.__enter__.return_value = pages.return_value
                write_pdf_report(root / "pdf" / "scan_report.pdf", report_pages)
                self.assertEqual(1, pages.return_value.savefig.call_count)
                self.assertNotIn("vein", report.velocity_per_beat)

    def test_velocity_page_embeds_maps_and_systole_artifacts(self) -> None:
        paths = {
            (kind, vessel): Path(f"{kind}_{vessel}.png")
            for kind in ("vessel_map", "systole")
            for vessel in ("artery", "vein")
        }
        image = np.zeros((8, 12, 3), dtype=np.uint8)
        with patch(
            "pipelines.pdf_report.pages.load_first_image", return_value=image
        ) as load:
            fig = build_report_pages(ReportData(), "scan", None, paths)[0]
        try:
            self.assertEqual(4, load.call_count)
            loaded = [call.args[0][0] for call in load.call_args_list]
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
            parameters={"ARI": {"value": 0.4, "unit": ""}},
            velocity_per_beat={"artery": values},
            velocity_units={"artery": "mm/s"},
        )
        fig = build_report_pages(report, "scan", None, {})[0]
        try:
            plot = next(
                ax for ax in fig.axes
                if ax.get_title() == "Artery filtered per-beat velocity"
            )
            mean = next(line for line in plot.lines if line.get_label() == "Mean beat")
            np.testing.assert_allclose(mean.get_ydata(), [1.5, 3.0, 6.0, 6.0])
            self.assertEqual("black", mean.get_color())
            self.assertIn("Mean per-beat RI = 0.40", [text.get_text() for text in plot.texts])
            map_plot = next(
                ax for ax in fig.axes if ax.get_title() == "Artery RMS-frequency map"
            )
            self.assertEqual(map_plot.title.get_fontsize(), plot.title.get_fontsize())
            self.assertEqual(map_plot.title.get_fontweight(), plot.title.get_fontweight())
            self.assertEqual(map_plot.title.get_fontfamily(), plot.title.get_fontfamily())
            systole_plot = next(
                ax for ax in fig.axes if ax.get_title() == "Artery systole indices"
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
        fig = build_report_pages(ReportData(parameters=parameters), "scan", None, {})[0]
        try:
            summary = next(
                ax for ax in fig.axes
                if any(text.get_text() == "Hemodynamic summary" for text in ax.texts)
            )
            texts = {text.get_text(): text for text in summary.texts}
            self.assertEqual("tab:red", texts["Artery"].get_color())
            self.assertEqual("tab:blue", texts["Vein"].get_color())
            self.assertLess(texts["1.00 mm/s"].get_position()[0], 0.5)
            self.assertGreater(texts["6.00 mm/s"].get_position()[0], 0.5)
            self.assertIn("Heart rate  |  72.5 bpm", texts)
            self.assertNotIn("99.00", texts)
            self.assertEqual(2, len(summary.patches))
        finally:
            plt.close(fig)

    def test_optional_page_names_its_results(self) -> None:
        report = ReportData(
            spectral=SpectralReportData(),
            volume_rates={
                "artery": VolumeRateReportData(
                    phase_curve=np.array([1.0]), overall_value=1.0,
                    valid_fraction=1.0, beat_count=1, unit="mm^3/s", method="mask",
                ),
                "vein": VolumeRateReportData(
                    phase_curve=np.array([2.0]), overall_value=2.0,
                    valid_fraction=1.0, beat_count=1, unit="mm^3/s", method="mask",
                ),
            },
        )
        pages = build_report_pages(report, "scan", None, {})
        try:
            labels = [text.get_text() for text in pages[1].texts]
            self.assertIn("Low-Frequency Spectra and Blood-Volume Rate", labels)
            self.assertNotIn("Additional results", " ".join(labels))
            flow_ax = next(
                ax for ax in pages[1].axes if ax.get_ylabel() == "Flow rate (mm$^3$/s)"
            )
            self.assertEqual(["tab:red", "tab:blue"], [line.get_color() for line in flow_ax.lines[:2]])
        finally:
            for page in pages:
                plt.close(page)

    def test_width_report_uses_red_arteries_and_blue_veins(self) -> None:
        def summary(label, value):
            return PhaseSummary(
                label=label, curve=np.array([value]), median=value,
                valid_fraction=1.0, branch_count=1, beat_count=1, unit="pixels",
            )

        fig = build_width_page(
            WidthReportData(
                displacement=(summary("Rigid Artery", 1.0), summary("Rigid Vein", 2.0)),
                gradient=(),
            ),
            "scan",
        )
        try:
            plot = next(ax for ax in fig.axes if ax.get_ylabel() == "Radial movement (pixels)")
            self.assertEqual(["tab:red", "tab:blue"], [line.get_color() for line in plot.lines])
        finally:
            plt.close(fig)

    def test_gradient_only_volume_rate_is_not_called_total_flow(self) -> None:
        schema = EyeFlowOutputPaths.active()
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "scan_EF.h5"
            with h5py.File(path, "w") as output_h5:
                output_h5.create_dataset(
                    schema.blood_volume_rate.artery.dynamic_edges,
                    data=np.full((4, 2, 3, 2), 6.0),
                )
            report = _read_report(path)
        artery = report.volume_rates["artery"]
        self.assertEqual("dynamic gradient-edge local", artery.method)
        self.assertEqual(6.0, artery.overall_value)
        self.assertEqual(1.0, artery.valid_fraction)

    def test_volume_headline_reduces_the_displayed_phase_curve(self) -> None:
        schema = EyeFlowOutputPaths.active()
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "scan_EF.h5"
            with h5py.File(path, "w") as output_h5:
                output_h5.create_dataset(
                    schema.blood_volume_rate.artery.total_masked_edges,
                    data=np.array([[1.0, np.nan], [3.0, 4.0]]),
                )
            rate = _read_report(path).volume_rates["artery"]
        np.testing.assert_allclose(rate.phase_curve, [1.0, 3.5])
        self.assertEqual(2.25, rate.overall_value)

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
                root / "pdf" / "scan_report.pdf", build_report_pages(report, "scan", None, {}),
            )
            content = pdf_path.read_bytes()
            self.assertTrue(content.startswith(b"%PDF-"))
            self.assertEqual(2, len(re.findall(rb"/Type /Page(?!s)", content)))


def _read_report(path: Path, *, processing_root="ProcessingMoments", include_optional=True):
    with h5py.File(path, "r") as h5file:
        output = PipelineH5Output(h5file, processing_root=processing_root)
        return extract_report_data(output, include_optional=include_optional)


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
