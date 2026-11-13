"""Marker-free HD/DV folder discovery, selection, and validation."""

from test_run_service import _descriptor
from test_stem_list_inputs import _write_companion_h5_files

from input_output import discover_run_inputs, is_run_directory
from pipeline_engine.run_service import expand_run_inputs, resolve_run_spec


def _run_folder(parent, stem="scan"):
    _write_companion_h5_files(parent, stem)
    return parent / stem


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
