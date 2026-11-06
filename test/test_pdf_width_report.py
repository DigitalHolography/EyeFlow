"""Focused contracts for the shared displacement and gradient PDF."""

from pathlib import Path
from types import SimpleNamespace
import re

import h5py
import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg", force=True)

from input_output.schema import EyeFlowOutputPaths
from input_output.schema.eyeflow_output import processing_path
from input_output.writers.h5_output import PipelineH5Output
from pipeline_engine.output_model import OutputModel
from pipeline_engine.storage import ProductStorage
from pipelines.pdf_report.width_data import extract_width_report_data
from pipelines.pdf_report.width_runner import run_width_pdf_report


def _write_width_outputs(output) -> None:
    paths = EyeFlowOutputPaths.active()
    displacement = np.array(
        [
            [[[1.0], [3.0]], [[np.nan], [5.0]]],
            [[[4.0], [4.0]], [[4.0], [4.0]]],
            [[[np.nan], [np.nan]], [[np.nan], [np.nan]]],
        ],
        dtype=np.float32,
    )
    output.write(
        f"{paths.displacement_profiles_root}/LevelSetMotion/Artery/"
        "Cross_sectional_radial_movement_amplitude_profile",
        displacement,
        unit="pixels",
    )
    output.write(
        f"{paths.spatial_gradient_metrics_root}/Artery/Transverse/Masked/tk/lumen_size",
        np.array([[8.0, np.nan], [10.0, 12.0], [np.nan, np.nan]], dtype=np.float32),
        unit="pixels",
    )
    output.write(
        f"{processing_path(paths.spatial_gradient_metrics_root, 'ProcessingBandRatio')}/"
        "Vein/Transverse/Masked/tk/lumen_size",
        np.full((3, 2), 6.0, dtype=np.float32),
        unit="pixels",
    )


def test_width_summary_matches_its_plotted_phase_curves(tmp_path) -> None:
    path = tmp_path / "result.h5"
    with h5py.File(path, "w") as h5file:
        _write_width_outputs(PipelineH5Output(h5file))
        report = extract_width_report_data(PipelineH5Output(h5file))
    assert report is not None
    assert [item.label for item in report.displacement] == ["LevelSetMotion Artery"]
    displacement = report.displacement[0]
    np.testing.assert_allclose(displacement.curve[:2], [3.0, 4.0])
    assert np.isnan(displacement.curve[2])
    assert displacement.median == 3.5
    assert displacement.valid_fraction == 7 / 12
    assert (displacement.branch_count, displacement.beat_count) == (2, 2)
    assert [item.label for item in report.gradient] == ["Moments Artery", "Band ratio Vein"]
    gradient = report.gradient[0]
    np.testing.assert_allclose(gradient.curve[:2], [8.0, 11.0])
    assert gradient.median == 9.5
    assert gradient.valid_fraction == 3 / 6
    assert gradient.beat_count is None


def test_width_summary_reads_assembled_output_before_publication(tmp_path) -> None:
    with ProductStorage(directory=tmp_path / "scratch") as storage:
        model = OutputModel(storage, tmp_path / "future.h5")
        output = PipelineH5Output(model)
        _write_width_outputs(output)
        report = extract_width_report_data(output)
        assert report is not None
        assert report.displacement[0].median == 3.5
        assert report.gradient[1].median == 6.0
        assert not Path(model.filename).exists()


def test_width_pdf_is_one_shared_file_and_skips_absent_results(tmp_path) -> None:
    path = tmp_path / "result.h5"
    pdf_dir = tmp_path / "pdf"
    manager = SimpleNamespace(
        layout=SimpleNamespace(stem="scan"),
        dir_for=lambda _kind: pdf_dir,
    )
    logs = []
    with h5py.File(path, "w") as h5file:
        ctx = SimpleNamespace(
            output=SimpleNamespace(h5=PipelineH5Output(h5file), manager=manager),
            log=logs.append,
            log_warning=logs.append,
        )
        run_width_pdf_report(ctx)
        assert not pdf_dir.exists()
        _write_width_outputs(ctx.output.h5)
        run_width_pdf_report(ctx)

    pdf_path = pdf_dir / "scan_width_report.pdf"
    content = pdf_path.read_bytes()
    assert content.startswith(b"%PDF-")
    assert len(re.findall(rb"/Type /Page(?!s)", content)) == 1
    assert len(list(pdf_dir.glob("*.pdf"))) == 1


@pytest.mark.parametrize("family", ["displacement", "gradient"])
def test_width_pdf_accepts_either_analysis_alone(tmp_path, family) -> None:
    paths = EyeFlowOutputPaths.active()
    pdf_dir = tmp_path / "pdf"
    manager = SimpleNamespace(
        layout=SimpleNamespace(stem="scan"),
        dir_for=lambda _kind: pdf_dir,
    )
    with h5py.File(tmp_path / "result.h5", "w") as h5file:
        output = PipelineH5Output(h5file)
        if family == "displacement":
            output.write(
                f"{paths.displacement_profiles_root}/LevelSetMotion/Vein/"
                "Cross_sectional_radial_movement_amplitude_profile",
                np.ones((3, 2, 2, 1), dtype=np.float32),
                unit="pixels",
            )
        else:
            output.write(
                f"{paths.spatial_gradient_metrics_root}/Vein/Transverse/Masked/tk/lumen_size",
                np.ones((3, 2), dtype=np.float32),
                unit="pixels",
            )
        ctx = SimpleNamespace(
            output=SimpleNamespace(h5=output, manager=manager),
            log=lambda _message: None,
            log_warning=lambda _message: None,
        )
        run_width_pdf_report(ctx)
    assert (pdf_dir / "scan_width_report.pdf").read_bytes().startswith(b"%PDF-")
