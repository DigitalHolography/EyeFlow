"""Focused contracts for moments report data, layout, and destinations."""

import re
from types import MappingProxyType, SimpleNamespace
from unittest.mock import patch

import h5py
import numpy as np

from input_output.output_manager import OutputManager, OutputType
from input_output.reports.pdf_report import write_pdf_report
from input_output.runs.layout import HoloRunLayout
from input_output.schema import EyeFlowOutputPaths
from input_output.writers.h5_output import PipelineH5Output
from pipeline_engine import PIPELINE_REGISTRY, PipelineDAG
from pipelines import load_pipeline_catalog
from pipelines.pdf_report.artifacts import export_moments_report_figures
from pipelines.pdf_report.data import ReportData, extract_report_data
from pipelines.pdf_report.pages import build_velocity_bvr_pages
from pipelines.pdf_report.runner import _html_report, _pdf_report, pdf_report_graph
from plotting.headless import close_figure


def test_report_plan_requires_mask_geometry_without_gradient_pipeline() -> None:
    load_pipeline_catalog()
    names = PipelineDAG(PIPELINE_REGISTRY.values()).resolve_targets(
        ["pdf_report"]
    ).names

    assert "blood_volume_rate" in names
    assert "spatial_gradient_moment0" not in names
    assert names.index("blood_volume_rate") < names.index("pdf_report")


def test_report_graph_keeps_pdf_html_and_figure_exports_independent() -> None:
    graph = pdf_report_graph()

    assert [stage.id for stage in graph.stages] == ["data", "figures", "pdf", "html"]
    assert all(stage.artifact for stage in graph.stages[1:])
    assert all(stage.requires == ("report_data",) for stage in graph.stages[1:])


def test_report_extracts_filtered_and_unfiltered_moments_data(tmp_path) -> None:
    schema = EyeFlowOutputPaths.active()
    path = tmp_path / "scan_EF.h5"
    with h5py.File(path, "w") as h5file:
        filtered = h5file.create_dataset(
            schema.artery_per_beat_filtered.velocity_signal,
            data=np.asarray([[-4.0, 4.0], [2.0, -2.0]]),
        )
        filtered.attrs["unit"] = "mm/s"
        h5file.create_dataset(
            schema.artery_per_beat.velocity_signal,
            data=np.asarray([[8.0, 9.0], [10.0, 11.0]]),
        )
        h5file.create_dataset(
            schema.cardiac_cycle.spectral_heart_rate_bpm,
            data=np.asarray(72.0),
        )

    with h5py.File(path, "r") as h5file:
        report = extract_report_data(
            PipelineH5Output(h5file, processing_root="ProcessingMoments"),
            include_bvr=False,
        )

    np.testing.assert_array_equal(
        report.velocity_per_beat["artery"], [[4.0, 4.0], [2.0, 2.0]]
    )
    np.testing.assert_array_equal(
        report.velocity_per_beat_unfiltered["artery"], [[8.0, 9.0], [10.0, 11.0]]
    )
    assert report.parameters["Average_Arterial_Velocity"]["value"] == 3.0
    assert report.parameters["heart_beat"]["value"] == 72.0
    assert report.velocity_units["artery"] == "mm/s"


def test_report_callbacks_use_shared_pdf_and_html_folders(tmp_path) -> None:
    manager = OutputManager(
        HoloRunLayout.from_holo(tmp_path / "scan.holo")
    ).for_artifact_namespace("doppler_moments")
    report = ReportData()
    messages = []
    ctx = SimpleNamespace(
        dag=SimpleNamespace(require=lambda name: report),
        output=SimpleNamespace(manager=manager),
        state=SimpleNamespace(get=lambda _key, default=None: default),
        log=messages.append,
    )

    with patch(
        "pipelines.pdf_report.runner.build_velocity_bvr_pages", return_value=[]
    ), patch("pipelines.pdf_report.runner.write_pdf_report") as write_pdf:
        _pdf_report(ctx)
    with patch(
        "pipelines.pdf_report.runner.build_velocity_bvr_html", return_value="<html></html>"
    ), patch("pipelines.pdf_report.runner.write_html_report") as write_html:
        _html_report(ctx)

    assert write_pdf.call_args.args[0] == (
        manager.layout.ef_dir / "pdf" / "scan_velocityBvrReport.pdf"
    )
    assert write_html.call_args.args[0] == (
        manager.layout.ef_dir / "html" / "scan_velocityBvrReport.html"
    )
    assert write_pdf.call_args.args[1] == []
    assert write_html.call_args.args[1] == "<html></html>"


def test_report_graphs_keep_requested_visual_contracts() -> None:
    filtered = np.asarray([[1.0, 2.0, 3.0], [2.0, 3.0, 4.0]])
    unfiltered = np.asarray(
        [[0.5, 1.5, 2.5], [1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
    )
    report = ReportData(
        velocity_per_beat={"artery": filtered, "vein": filtered},
        velocity_per_beat_unfiltered={"artery": unfiltered, "vein": unfiltered},
        velocity_units={"artery": "mm/s", "vein": "mm/s"},
    )
    pages = build_velocity_bvr_pages(report, "scan", None, {})
    try:
        velocity = next(
            ax for ax in pages[0].axes if ax.get_title() == "Artery Per-Beat Velocity"
        )
        systole = next(
            ax for ax in pages[0].axes if ax.get_title() == "Artery Systole Indices"
        )
        lumen = next(
            ax for ax in pages[1].axes
            if ax.get_title() == "Artery Lumen-Diameter Distribution"
        )
        volume_rate = next(
            ax for ax in pages[1].axes
            if ax.get_title() == "Artery Mask-Derived Blood-Volume Rate"
        )

        assert {line.get_label() for line in velocity.lines} == {
            "Mean beat", "Unfiltered median"
        }
        assert all(line.get_color() == "black" for line in velocity.lines)
        assert velocity.get_legend() is None
        np.testing.assert_allclose(
            velocity.get_position().bounds[2:], systole.get_position().bounds[2:]
        )
        np.testing.assert_allclose(
            lumen.get_position().bounds[2:], volume_rate.get_position().bounds[2:]
        )
    finally:
        for page in pages:
            close_figure(page)


def test_page_one_summary_accepts_managed_immutable_parameters() -> None:
    def parameter(value, unit=""):
        return MappingProxyType({"value": value, "unit": unit})

    parameters = MappingProxyType({
        "Min_Arterial_Velocity": parameter(1.0, "mm/s"),
        "Average_Arterial_Velocity": parameter(2.0, "mm/s"),
        "Max_Arterial_Velocity": parameter(3.0, "mm/s"),
        "ARI": parameter(0.4),
        "API": parameter(0.5),
        "Min_Venous_Velocity": parameter(4.0, "mm/s"),
        "Average_Venous_Velocity": parameter(5.0, "mm/s"),
        "Max_Venous_Velocity": parameter(6.0, "mm/s"),
        "VRI": parameter(0.6),
        "VPI": parameter(0.7),
        "heart_beat": parameter(72.0, "bpm"),
    })
    pages = build_velocity_bvr_pages(
        ReportData(parameters=parameters), "scan", None, {}
    )
    try:
        summary = next(
            axis for axis in pages[0].axes
            if any(text.get_text() == "Summary" for text in axis.texts)
        )
        text = {item.get_text() for item in summary.texts}
        assert "Unavailable" not in text
        assert "Heart rate  |  72.0 bpm" in text
        assert "2.00 $\\mathrm{mm.s^{-1}}$" in text
        assert "5.00 $\\mathrm{mm.s^{-1}}$" in text
    finally:
        for page in pages:
            close_figure(page)


def test_report_native_graph_exports_matching_png_and_eps(tmp_path) -> None:
    manager = OutputManager(
        HoloRunLayout.from_holo(tmp_path / "scan.holo")
    ).for_artifact_namespace("doppler_moments")
    report = ReportData(
        velocity_per_beat={
            "artery": np.asarray([[1.0, 2.0], [2.0, 3.0]]),
            "vein": np.asarray([[3.0, 2.0], [2.0, 1.0]]),
        }
    )

    paths = export_moments_report_figures(manager, report)

    assert len(paths) == 4
    assert {path.suffix for path in paths} == {".png", ".eps"}
    assert all(path.is_file() for path in paths)


def test_renders_two_page_pdf_with_optional_products(tmp_path) -> None:
    schema = EyeFlowOutputPaths.active()
    h5_path = tmp_path / "scan_EF.h5"
    with h5py.File(h5_path, "w") as h5file:
        h5file.create_dataset(
            schema.artery_per_beat_filtered.velocity_signal,
            data=np.asarray([[1.0, 2.0, 3.0, 2.0], [2.0, 3.0, 4.0, 3.0]]),
        )
        h5file.create_dataset(
            f"{schema.waveform_shape_metrics_root}/artery/global/raw/RI",
            data=np.asarray([0.75, 0.75]),
        )
        h5file.create_dataset(
            schema.blood_volume_rate.artery.total_masked_edges,
            data=np.asarray([[1.0, 2.0], [2.0, 3.0], [3.0, 4.0], [2.0, 3.0]]),
        )

    with h5py.File(h5_path, "r") as h5file:
        report = extract_report_data(
            PipelineH5Output(h5file, processing_root="ProcessingMoments")
        )
    pdf_path = write_pdf_report(
        tmp_path / "pdf" / "scan_velocityBvrReport.pdf",
        build_velocity_bvr_pages(report, "scan", None, {}),
    )

    content = pdf_path.read_bytes()
    assert content.startswith(b"%PDF-")
    assert len(re.findall(rb"/Type /Page(?!s)", content)) == 2
