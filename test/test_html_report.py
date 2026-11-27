"""Focused contracts for the self-contained moments HTML report."""

from __future__ import annotations

from unittest.mock import patch

import matplotlib
import numpy as np

matplotlib.use("Agg", force=True)

from input_output.output_manager import OutputManager, OutputType
from input_output.reports import write_html_report
from input_output.runs.layout import HoloRunLayout
from pipelines.pdf_report.assets import load_analysis_artifact
from pipelines.pdf_report.data import (
    ReportData,
    SystoleReportData,
    VesselGeometryReportData,
    VolumeRateReportData,
)
from pipelines.pdf_report.html import build_velocity_bvr_html


def test_html_velocity_legend_uses_existing_colorbar_artifact(tmp_path) -> None:
    image = np.zeros((4, 12, 3), dtype=np.uint8)

    with patch(
        "pipelines.pdf_report.assets.load_first_image", return_value=image,
    ) as load_first_image:
        actual = load_analysis_artifact(
            tmp_path, "scan", "velocity_colorbar", "artery"
        )

    assert actual is image
    assert load_first_image.call_args.args == (
        (tmp_path / "scan_colorbar_v_Artery.png",),
    )


def test_html_report_copies_current_content_without_pdf_dependency() -> None:
    beats = np.asarray([[1.0, 2.0, 3.0], [1.5, 2.5, 3.5]])
    systole = SystoleReportData(
        velocity=np.asarray([1.0, 3.0, 2.0, 1.0]),
        systolic_boundaries=np.asarray([0, 3]),
        maxima=np.asarray([1]),
        minima=np.asarray([2]),
        sample_interval_seconds=0.1,
        unit="mm/s",
    )
    rate = VolumeRateReportData(
        beat_curves=np.asarray([[1.0, 1.2], [2.0, 2.2]]),
        phase_curve=np.asarray([1.1, 2.1]),
        overall_value=1.6,
        unit="mm^3/s",
    )
    report = ReportData(
        parameters={
            "Min_Arterial_Velocity": {"value": 1.0, "unit": "mm/s"},
            "Average_Arterial_Velocity": {"value": 2.0, "unit": "mm/s"},
            "Max_Arterial_Velocity": {"value": 3.0, "unit": "mm/s"},
            "ARI": {"value": 0.4, "unit": ""},
            "API": {"value": 0.5, "unit": ""},
            "Min_Venous_Velocity": {"value": 1.5, "unit": "mm/s"},
            "Average_Venous_Velocity": {"value": 2.5, "unit": "mm/s"},
            "Max_Venous_Velocity": {"value": 3.5, "unit": "mm/s"},
            "VRI": {"value": 0.6, "unit": ""},
            "VPI": {"value": 0.7, "unit": ""},
            "heart_beat": {"value": 72.5, "unit": "bpm"},
        },
        velocity_per_beat={"artery": beats, "vein": beats},
        velocity_per_beat_unfiltered={"artery": beats, "vein": beats},
        velocity_units={"artery": "mm/s", "vein": "mm/s"},
        systole_signals={"artery": systole, "vein": systole},
        mask_volume_rates={"artery": rate, "vein": rate},
        vessel_geometry={
            "artery": VesselGeometryReportData(3, 12, 84.0, 80.0, 7.0),
            "vein": VesselGeometryReportData(4, 15, 96.0, 91.0, 8.0),
        },
    )
    image = np.zeros((50, 80, 3), dtype=np.uint8)

    with patch(
        "pipelines.pdf_report.html.load_velocity_artifact", return_value=image,
    ), patch(
        "pipelines.pdf_report.html.load_analysis_artifact", return_value=image,
    ):
        document = build_velocity_bvr_html(report, "scan & eye", None, {})

    assert document.startswith("<!doctype html>")
    assert "scan &amp; eye" in document
    assert 'id="velocity"' in document
    assert 'id="vessel-analysis"' in document
    assert 'class="report-section"' not in document
    assert "break-after" not in document
    assert "270mm" not in document
    assert "Artery Per-Beat Velocity" in document
    assert "Vein Mask-Derived Blood-Volume Rate" in document
    assert "Minimum filtered velocity" in document
    assert "Lumen diameter SD" in document
    assert 'class="paired-graphs"' not in document
    artery_velocity = document.split('<article class="vessel-column">', 1)[1].split(
        "</article>", 1
    )[0]
    assert artery_velocity.index("Artery Velocity Map") < artery_velocity.index(
        "Artery Velocity Distribution over Time"
    )
    assert artery_velocity.index("Artery Velocity Distribution over Time") < (
        artery_velocity.index("Artery velocity legend")
    )
    assert artery_velocity.index("Artery velocity legend") < artery_velocity.index(
        "Artery Systole Indices"
    )
    artery_analysis = document.split('id="vessel-analysis"', 1)[1].split(
        '<article class="vessel-column">', 1
    )[1].split("</article>", 1)[0]
    assert "Artery Velocity Distribution over Time" not in artery_analysis
    assert "Artery Branch Identity and Sampling Regions" in artery_analysis
    assert artery_analysis.index(
        "Artery Branch Identity and Sampling Regions"
    ) < artery_analysis.index("Artery Lumen-Diameter Distribution")
    assert "72.5 bpm" in document
    assert "mm.s⁻¹" in document
    assert "mm³.s⁻¹" in document
    assert document.count('class="graph-legend"') == 2
    assert document.count('src="data:image/png;base64,') == 16
    assert ".pdf" not in document.lower()
    assert "<iframe" not in document.lower()


def test_html_writer_uses_the_dedicated_output_folder(tmp_path) -> None:
    manager = OutputManager(HoloRunLayout.from_holo(tmp_path / "scan.holo"))
    path = manager.layout.ef_dir / OutputType.HTML.value / "scan.html"

    written = write_html_report(path, "<!doctype html><title>scan</title>")

    assert written == manager.layout.ef_dir / "html" / "scan.html"
    assert written.read_text(encoding="utf-8") == "<!doctype html><title>scan</title>"
