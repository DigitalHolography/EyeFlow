"""Marker-free HD/DV folder discovery, selection, and validation."""

import pytest
from test_run_service import _descriptor
from test_stem_list_inputs import _write_companion_h5_files

from input_output import discover_run_inputs, holo_input_status, is_run_directory
from pipeline_engine.run_service import expand_run_inputs, resolve_run_spec


def _run_folder(parent, stem="scan"):
    _write_companion_h5_files(parent, stem)
    return parent / stem


def test_discovery_accepts_the_selected_run_folder(tmp_path):
    stem = "scan"
    folder = _run_folder(tmp_path / "nested", stem)
    assert is_run_directory(folder)
    assert discover_run_inputs(folder) == [folder]
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


def test_cli_expansion_and_output_routing_do_not_add_an_extra_run_directory(tmp_path):
    folder = _run_folder(tmp_path / "data")
    source = folder.parent
    with expand_run_inputs(source) as expanded:
        assert expanded.paths == (folder,)
        assert expanded.batch_root == folder.parent
        spec = resolve_run_spec(input_paths=expanded.paths, target_names=["sample"],
                                pipelines=[_descriptor()], batch_root=expanded.batch_root,
                                output_root=tmp_path / "results")
    assert spec.requests[0].output_manager.layout.ef_dir == tmp_path / "results/scan/scan_EF"
    assert spec.requests[0].output_manager.layout.input_path == folder


def test_discovery_defers_hdf5_validation_and_keeps_existing_outputs_safe(tmp_path):
    folder = _run_folder(tmp_path)
    hd = folder / "scan_HD/h5/scan_HD_output.h5"
    hd.unlink()
    existing = folder / "scan_EF" / "previous.txt"
    existing.parent.mkdir()
    existing.write_text("previous result", encoding="utf-8")
    assert discover_run_inputs(folder) == [folder]
    assert not holo_input_status(folder).hd
    assert holo_input_status(folder).dv
    with pytest.raises(FileNotFoundError, match="HDF5"):
        resolve_run_spec(input_paths=[folder], target_names=["sample"], pipelines=[_descriptor()])
    assert existing.read_text(encoding="utf-8") == "previous result"
