"""JSON workflow selection, validation, and downstream execution contracts."""

import json
from unittest.mock import patch

import h5py
import pytest

from app_settings import AppSettingsStore, selected_velocity_methods
from pipeline_engine.run_service import execute_run, resolve_run_spec
from pipelines import load_pipeline_catalog
from pipelines.velocity.sources import load_velocity_inputs
from test_dual_velocity_workflows import _dual_input


@pytest.mark.parametrize("selection, methods", [
    ("doppler_moments", ("doppler_moments",)),
    ("bandratio", ("frequency_bands",)),
    ("both", ("doppler_moments", "frequency_bands")),
    ("frequency_bands", ("frequency_bands",)),
])
def test_settings_selection_survives_import_save_and_reload(tmp_path, selection, methods):
    source = tmp_path / "import.json"
    source.write_text(json.dumps({"velocity_estimation_method": selection}))
    store = AppSettingsStore(path=tmp_path / "settings.json", default_template_path=source)
    assert selected_velocity_methods(store.load_velocity_estimation_method()) == methods
    store.import_file(source)
    store.save_ui_mode("advanced")
    reopened = AppSettingsStore(path=store.path)
    assert reopened.load()["velocity_estimation_method"] == selection
    assert selected_velocity_methods(reopened.load_velocity_estimation_method()) == methods


def test_missing_selection_defaults_to_both(tmp_path):
    store = AppSettingsStore(path=tmp_path / "settings.json")
    assert store.load_velocity_estimation_method() == "both"
    store.save({"ui_mode": "minimal"})
    assert store.load_velocity_estimation_method() == "both"


@pytest.mark.parametrize("value", [None, "unknown", "", [], {}, True, 1])
def test_invalid_selection_preserves_active_settings_and_rejects_run(tmp_path, value):
    store = AppSettingsStore(path=tmp_path / "settings.json")
    store.save({"velocity_estimation_method": "both"})
    source = tmp_path / "invalid.json"
    source.write_text(json.dumps({"velocity_estimation_method": value}))
    with pytest.raises(ValueError, match="velocity_estimation_method"):
        store.import_file(source)
    assert store.load_velocity_estimation_method() == "both"
    with pytest.raises(ValueError, match="velocity_estimation_method"):
        resolve_run_spec(input_paths=[], target_names=[], pipelines=[],
                         velocity_estimation_method=value)


@pytest.mark.parametrize("selection, roots", [
    ("doppler_moments", {"ProcessingMoments"}),
    ("bandratio", {"ProcessingBandRatio"}),
    ("both", {"ProcessingMoments", "ProcessingBandRatio"}),
])
def test_only_selected_workflows_load_inputs_and_execute_downstream(tmp_path, selection, roots):
    holo, _ = _dual_input(tmp_path)
    available, _ = load_pipeline_catalog()
    spec = resolve_run_spec(
        input_paths=[holo], target_names=["velocity_analysis"], pipelines=available,
        pipeline_options={"velocity_analysis": ()}, velocity_estimation_method=selection,
    )
    with patch("pipelines.velocity.runner.load_velocity_inputs", wraps=load_velocity_inputs) as load:
        result = execute_run(spec)
    assert result.succeeded, result.failures
    methods = selected_velocity_methods(selection)
    assert [call.args[0].velocity_estimation_method for call in load.call_args_list] == list(methods)
    with h5py.File(result.outputs[0], "r") as output:
        assert set(output) & {"ProcessingMoments", "ProcessingBandRatio"} == roots
        assert list(output.attrs["velocity_estimation_methods"]) == list(methods)
        statuses = json.loads(output.attrs["pipeline_status"])
        for method in methods:
            assert statuses["variants"][method]["velocity_analysis"] == "succeeded"
        assert set(statuses["variants"]) == set(methods)
        assert json.loads(output.attrs["execution_variant_failures"]) == {}
