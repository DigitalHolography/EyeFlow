"""Marker-free HD/DV folders share discovery and validation across entry points."""

from unittest.mock import patch

import h5py
import pytest
from test_input_selection import _app
from test_run_service import _descriptor, _write_input
from test_stem_list_inputs import _write_companion_h5_files

from input_output import discover_run_inputs, holo_input_status, is_run_directory
from pipeline_engine.run_service import execute_run, expand_run_inputs, resolve_run_spec
from ui.controllers.run import RunController


def _run_folder(parent, stem="scan"):
    _write_companion_h5_files(parent, stem)
    return parent / stem


@pytest.mark.parametrize("stem, select_parent", [
    ("scan", False), ("scan", True), ("scan.with.dots", True),
    ("scan.txt", False), ("scan.holo", False),
])
def test_discovery_includes_the_selected_folder_itself_or_nested_runs(tmp_path, stem, select_parent):
    folder = _run_folder(tmp_path / "nested", stem)
    assert is_run_directory(folder)
    assert discover_run_inputs(tmp_path if select_parent else folder) == [folder]
    status = holo_input_status(folder)
    assert status.hd and status.dv
    spec = resolve_run_spec(input_paths=[folder], target_names=["sample"], pipelines=[_descriptor()])
    layout = spec.requests[0].input_layout
    assert layout.root_dir == folder
    assert layout.stem == stem
    assert layout.input_path == folder
    assert layout.holo_path == folder.parent / f"{stem}.holo"
    assert spec.requests[0].output_manager.layout.ef_dir == folder / f"{stem}_EF"
    assert not layout.holo_path.exists()


def test_discovery_requires_matching_hd_and_dv_and_ignores_ef_ae_only(tmp_path):
    good = _run_folder(tmp_path / "good")
    for name, children in (
        ("outputs_only", ["outputs_only_EF", "outputs_only_AE"]),
        ("missing_dv", ["missing_dv_HD"]),
        ("wrong_stem", ["other_HD", "other_DV"]),
    ):
        root = tmp_path / name
        for child in children:
            (root / child).mkdir(parents=True)
        assert not is_run_directory(root)
    assert discover_run_inputs(tmp_path) == [good]
    (good / "scan_EF").mkdir()
    (good / "scan_AE").mkdir()
    assert discover_run_inputs(tmp_path) == [good]


def test_real_markers_and_their_run_directories_are_discovered_once(tmp_path):
    marker = _write_input(tmp_path / "with marker")
    folder = _run_folder(tmp_path / "without marker")
    assert set(discover_run_inputs(tmp_path)) == {marker, folder}
    # A selected data folder can use its real sibling marker outside the search root.
    assert discover_run_inputs(marker.parent / marker.stem) == [marker]
    app = _app()
    app.input_controller.add_holo_input_paths([marker, tmp_path])
    assert app.input_controller.selected_holo_paths() == [marker, folder]


@pytest.mark.parametrize("entry", ["native", "drop"])
def test_folder_inputs_reach_gui_status_and_execution_without_creating_markers(tmp_path, entry):
    folder = _run_folder(tmp_path / "subject")
    ae = folder / "scan_AE" / "keep.txt"
    ae.parent.mkdir()
    ae.write_text("unrelated output", encoding="utf-8")
    app = _app()
    if entry == "native":
        app.ui_services.file_dialogs.askdirectories.return_value = (str(folder.parent),)
        controller = RunController.__new__(RunController)
        controller.app = app
        controller.choose_holo_folder()
    else:
        assert app.input_controller.handle_dropped_paths([folder.parent])
    assert app.input_controller.selected_holo_paths() == [folder]
    app.input_controller.update_minimal_path_labels()
    assert app.holo_hd_status_var.get() == "HD found"
    assert app.holo_dv_status_var.get() == "DV found"
    spec = resolve_run_spec(input_paths=[folder], target_names=["sample"], pipelines=[_descriptor()])
    observed = []
    result = execute_run(spec, on_file_start=lambda path, *_: observed.append(path))
    assert result.succeeded, result.failures
    assert observed == [folder]
    assert result.outputs == (folder / "scan_EF/h5/scan_EF.h5",)
    assert ae.read_text(encoding="utf-8") == "unrelated output"
    assert not (folder.parent / "scan.holo").exists()


@pytest.mark.parametrize("select_parent", [False, True])
def test_cli_expansion_and_output_routing_do_not_add_an_extra_run_directory(tmp_path, select_parent):
    folder = _run_folder(tmp_path / "data")
    source = folder.parent if select_parent else folder
    with expand_run_inputs(source) as expanded:
        assert expanded.paths == (folder,)
        assert expanded.batch_root == folder.parent
        spec = resolve_run_spec(input_paths=expanded.paths, target_names=["sample"],
                                pipelines=[_descriptor()], batch_root=expanded.batch_root,
                                output_root=tmp_path / "results")
    assert spec.requests[0].output_manager.layout.ef_dir == tmp_path / "results/scan/scan_EF"
    assert spec.requests[0].output_manager.layout.input_path == folder


def test_cli_can_execute_a_marker_free_folder(tmp_path):
    import cli

    folder = _run_folder(tmp_path / "data")
    targets = tmp_path / "pipelines.txt"
    targets.write_text("sample\n", encoding="utf-8")
    with patch("cli._build_pipeline_registry", return_value={"sample": _descriptor()}):
        assert cli.run_cli(folder, pipelines_file=targets) == 0
    assert (folder / "scan_EF/h5/scan_EF.h5").is_file()
    assert not (folder.parent / "scan.holo").exists()


@pytest.mark.parametrize("invalid_h5", ["missing", "ambiguous"])
def test_discovery_defers_hdf5_validation_and_keeps_existing_outputs_safe(tmp_path, invalid_h5):
    folder = _run_folder(tmp_path)
    hd = folder / "scan_HD/h5/scan_HD_output.h5"
    hd.unlink()
    if invalid_h5 == "ambiguous":
        for name in ("one.h5", "two.h5"):
            with h5py.File(hd.parent / name, "w"):
                pass
    existing = folder / "scan_EF" / "previous.txt"
    existing.parent.mkdir()
    existing.write_text("previous result", encoding="utf-8")
    assert discover_run_inputs(folder) == [folder]
    assert not holo_input_status(folder).hd
    assert holo_input_status(folder).dv
    with pytest.raises(FileNotFoundError, match="HDF5"):
        resolve_run_spec(input_paths=[folder], target_names=["sample"], pipelines=[_descriptor()])
    assert existing.read_text(encoding="utf-8") == "previous result"


def test_failed_marker_free_run_reports_the_actual_folder(tmp_path):
    folder = _run_folder(tmp_path)
    spec = resolve_run_spec(input_paths=[folder], target_names=["sample"], pipelines=[_descriptor()])
    with patch("pipeline_engine.run_service.run_pipelines_to_output", side_effect=ValueError("failed")):
        result = execute_run(spec)
    assert not result.succeeded
    assert result.failures[0].input_path == folder
