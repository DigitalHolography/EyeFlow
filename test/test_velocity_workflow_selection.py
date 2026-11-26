"""JSON workflow selection, validation, and downstream execution contracts."""

import json
import pytest

from app_settings import UserSettings, selected_velocity_methods
from pipeline_engine.execution.run_service import resolve_run_spec


@pytest.mark.parametrize("selection, methods", [
    ("doppler_moments", ("doppler_moments",)),
    ("bandratio", ("frequency_bands",)),
    ("both", ("doppler_moments", "frequency_bands")),
    ("frequency_bands", ("frequency_bands",)),
])
def test_settings_selection_survives_import_save_and_reload(tmp_path, selection, methods):
    source = tmp_path / "import.json"
    source.write_text(json.dumps({"velocity_estimation_method": selection}))
    store = UserSettings(path=tmp_path / "settings.json", default_template_path=source)
    assert selected_velocity_methods(store.velocity_estimation_method) == methods
    store.import_file(source)
    store.ui_mode = "advanced"
    store.save()
    reopened = UserSettings(path=store.path)
    assert reopened.load()["velocity_estimation_method"] == selection
    assert selected_velocity_methods(reopened.velocity_estimation_method) == methods


def test_missing_selection_defaults_to_both(tmp_path):
    store = UserSettings(path=tmp_path / "settings.json")
    assert store.velocity_estimation_method == "both"
    store.save({"ui_mode": "minimal"})
    assert store.velocity_estimation_method == "both"


@pytest.mark.parametrize("value", [None, "unknown", "", [], {}, True, 1])
def test_invalid_selection_preserves_active_settings_and_rejects_run(tmp_path, value):
    store = UserSettings(path=tmp_path / "settings.json")
    store.save({"velocity_estimation_method": "both"})
    source = tmp_path / "invalid.json"
    source.write_text(json.dumps({"velocity_estimation_method": value}))
    with pytest.raises(ValueError, match="velocity_estimation_method"):
        store.import_file(source)
    assert store.velocity_estimation_method == "both"
    with pytest.raises(ValueError, match="velocity_estimation_method"):
        resolve_run_spec(input_paths=[], target_names=[], pipelines=[],
                         settings=UserSettings(values={"velocity_estimation_method": value}))
