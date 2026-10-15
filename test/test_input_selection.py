"""Input selection retains each run's directory across GUI entry points."""

import os
import tkinter as tk
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock, patch

import pytest
from test_run_service import _descriptor, _write_input
from test_stem_list_inputs import _Var

from pipeline_engine.run_service import (
    _batch_root,
    _relative_to_batch,
    execute_run,
    resolve_run_spec,
)
from ui.controllers.input import InputController
from ui.controllers.run import RunController
from ui.services import TkFileDialogService


@pytest.fixture(scope="module", autouse=True)
def tcl_interpreter():
    return tk.Tcl()


def _app():
    app = SimpleNamespace(
        _selected_holo_input_paths=[], _synchronizing_holo_input_var=False,
        holo_input_var=_Var(), minimal_holo_input_path_var=_Var(),
        holo_hd_status_var=_Var(), holo_dv_status_var=_Var(),
        _muted_fg="muted", _error_color="error", _success_color="success",
        progress_controller=Mock(),
        ui_services=SimpleNamespace(dialogs=Mock(), file_dialogs=Mock()),
    )
    app.input_controller = InputController(app)
    return app


def _inputs(tmp_path):
    first = _write_input(tmp_path / "drive Z" / "260618", "same_stem")
    second = _write_input(tmp_path / "drive Y" / "260129", "same_stem")
    return first, second


def _check_run_selection(app, holos):
    assert app.input_controller.selected_holo_paths() == list(holos)
    app.input_controller.update_minimal_path_labels()
    assert app.holo_hd_status_var.get() == "Found 2/2 HD"
    assert app.holo_dv_status_var.get() == "Found 2/2 DV"
    spec = resolve_run_spec(input_paths=app.input_controller.selected_holo_paths(),
                            target_names=["sample"], pipelines=[_descriptor()])
    assert [item.input_layout.holo_path for item in spec.requests] == list(holos)
    assert [item.output_manager.layout.ef_dir for item in spec.requests] == [
        holo.parent / holo.stem / f"{holo.stem}_EF" for holo in holos
    ]
    return spec


def test_popup_preserves_each_returned_file_path(tmp_path):
    holos = _inputs(tmp_path)
    app = _app()
    app.ui_services.file_dialogs.askopenfilenames.return_value = tuple(map(str, holos))
    controller = RunController.__new__(RunController)
    controller.app = app
    controller.choose_holo_file()
    _check_run_selection(app, holos)
    app.ui_services.file_dialogs.askopenfilenames.assert_called_once_with(
        filetypes=[("HOLO or path list", "*.holo *.txt"), ("All files", "*.*")],
        initialdir=os.getcwd(), title="Select .holo file(s) or one .txt path list",
    )
    app.ui_services.dialogs.showerror.assert_not_called()


def test_file_dialog_service_delegates_to_the_native_tk_dialog():
    with patch("ui.services.filedialog.askopenfilenames", return_value=("scan.holo",)) as native:
        result = TkFileDialogService().askopenfilenames(title="Select files")
    assert result == ("scan.holo",)
    native.assert_called_once_with(title="Select files")


def test_folder_dialog_service_delegates_to_the_native_tk_dialog():
    with patch("ui.services.filedialog.askdirectory", return_value="folder") as native:
        result = TkFileDialogService().askdirectory(title="Select folder", mustexist=True)
    assert result == "folder"
    native.assert_called_once_with(title="Select folder", mustexist=True)


@pytest.mark.parametrize("folder_first", [False, True])
def test_native_file_and_recursive_folder_selections_coexist(tmp_path, folder_first):
    first, second = _inputs(tmp_path)
    app = _app()
    controller = RunController.__new__(RunController)
    controller.app = app
    app.ui_services.file_dialogs.askopenfilenames.return_value = (str(first),)
    # The selected folder contains the second run one directory below it.
    app.ui_services.file_dialogs.askdirectory.return_value = str(second.parent.parent)
    if folder_first:
        controller.choose_holo_folder()
        controller.choose_holo_file()
    else:
        controller.choose_holo_file()
        controller.choose_holo_folder()
    expected = [second, first] if folder_first else [first, second]
    _check_run_selection(app, expected)
    controller.choose_holo_folder()  # Overlapping selections must not duplicate runs.
    _check_run_selection(app, expected)
    app.ui_services.file_dialogs.askdirectory.assert_called_with(
        initialdir=str(expected[0].parent),
        title="Select folder containing .holo files or HD/DV runs (including subfolders)",
        mustexist=True,
    )
    app.ui_services.dialogs.showwarning.assert_not_called()


@pytest.mark.parametrize("empty_folder", [False, True])
def test_cancelled_or_empty_folder_choice_preserves_selected_files(tmp_path, empty_folder):
    holos = _inputs(tmp_path)
    app = _app()
    app.input_controller.assign_holo_input_paths(holos)
    folder = tmp_path / "empty"
    folder.mkdir()
    app.ui_services.file_dialogs.askdirectory.return_value = str(folder) if empty_folder else ""
    controller = RunController.__new__(RunController)
    controller.app = app
    controller.choose_holo_folder()
    assert app.input_controller.selected_holo_paths() == list(holos)
    assert app.ui_services.dialogs.showwarning.call_count == int(empty_folder)


def test_clear_resets_the_combined_selection(tmp_path):
    app = _app()
    app.input_controller.assign_holo_input_paths(_inputs(tmp_path))
    controller = RunController.__new__(RunController)
    controller.app = app
    controller.clear_holo_inputs()
    app.input_controller.update_minimal_path_labels()
    assert app.input_controller.selected_holo_paths() == []
    assert app.holo_hd_status_var.get() == app.holo_dv_status_var.get() == ""


def test_text_lists_remain_a_standalone_input_when_adding_files_or_folders(tmp_path):
    first, second = _inputs(tmp_path)
    path_list = tmp_path / "inputs.txt"
    path_list.write_text(str(second), encoding="utf-8")
    app = _app()
    app.input_controller.add_holo_input_paths([first])
    app.input_controller.add_holo_input_paths([path_list])
    assert app.input_controller.selected_holo_paths() == [path_list]
    app.input_controller.add_holo_input_paths([second.parent])
    assert app.input_controller.selected_holo_paths() == [second]
    with pytest.raises(ValueError, match="by itself"):
        app.input_controller.add_holo_input_paths([first, path_list])
    assert app.input_controller.selected_holo_paths() == [second]


def test_popup_cancel_keeps_previous_selection(tmp_path):
    holos = _inputs(tmp_path)
    app = _app()
    app.input_controller.assign_holo_input_paths(holos)
    app.ui_services.file_dialogs.askopenfilenames.return_value = ()
    controller = RunController.__new__(RunController)
    controller.app = app
    controller.choose_holo_file()
    assert app.input_controller.selected_holo_paths() == list(holos)



@pytest.mark.parametrize("kind", ["files", "folders", "mixed"])
def test_drop_handles_multiple_locations_folders_and_tcl_paths_with_spaces(
    tmp_path, kind, tcl_interpreter,
):
    first, second = _inputs(tmp_path)
    dropped = {"files": [first, second], "folders": [first.parent, second.parent],
               "mixed": [first, second.parent, first.parent]}[kind]
    app = _app()
    app.tk = tcl_interpreter
    event = SimpleNamespace(data=" ".join("{" + str(path) + "}" for path in dropped))
    app.input_controller.on_input_drop(event)
    _check_run_selection(app, [first, second])
    app.ui_services.dialogs.showwarning.assert_not_called()


def test_quoted_txt_and_manual_selection_execute_the_same_independent_runs(tmp_path):
    first, second = _inputs(tmp_path)
    input_list = tmp_path / "inputs.txt"
    input_list.write_text(f'"{first}"\n"{second}"\n', encoding="utf-8")
    manual = resolve_run_spec(input_paths=[first, second], target_names=["sample"],
                              pipelines=[_descriptor()])
    listed = resolve_run_spec(input_paths=[input_list], target_names=["sample"],
                              pipelines=[_descriptor()])
    assert [request.input_layout for request in manual.requests] == [
        request.input_layout for request in listed.requests
    ]
    result = execute_run(listed)
    assert result.succeeded, result.failures
    assert len(result.outputs) == 2
    assert result.outputs[0] != result.outputs[1]


@pytest.mark.skipif(os.name != "nt", reason="Windows drive-letter semantics")
def test_output_routing_keeps_different_drive_prefixes_without_a_common_parent():
    holos = [Path("Z:/260618/same_stem.holo"), Path("Y:/260129/same_stem.holo")]
    batch_root = _batch_root(holos)
    assert [_relative_to_batch(holo, batch_root) for holo in holos] == [
        Path("Z/260618/same_stem.holo"), Path("Y/260129/same_stem.holo"),
    ]
