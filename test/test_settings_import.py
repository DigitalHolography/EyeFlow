"""Tests for importing EyeFlow settings through the Tools menu action."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock, patch
import pytest

from app_settings import (
    AppSettingsStore,
    BooleanSettings,
    PipelineOptionsSettings,
    PipelineVisibilitySettings,
    UserSettings,
    user_settings,
)
from ui.controllers.settings import SettingsController


@pytest.mark.parametrize("settings,count", [
    ({}, 42),
    ({"number_of_radii_in_half_diagonal": 16}, 16),
])
def test_saved_radius_count_and_default(tmp_path, settings, count):
    path = tmp_path / "import.json"
    path.write_text(json.dumps(settings), encoding="utf-8")
    store = UserSettings(path=tmp_path / "active.json")
    store.import_file(path)
    assert store.number_of_radii_in_half_diagonal == count
    assert store.load() == settings


@pytest.mark.parametrize("value", [0, -1, 2.5, None, float("nan"), float("inf"), True, "4"])
def test_invalid_radius_count_does_not_replace_saved_configuration(tmp_path, value):
    store = UserSettings(path=tmp_path / "active.json")
    original = {"number_of_radii_in_half_diagonal": 5}
    store.save(original)
    source = tmp_path / "invalid.json"
    source.write_text(json.dumps({"number_of_radii_in_half_diagonal": value}), encoding="utf-8")
    with pytest.raises(ValueError, match="number_of_radii_in_half_diagonal"):
        store.import_file(source)
    assert store.load() == original


def test_gui_and_cli_import_the_same_settings_dataset():
    import cli
    from ui import app

    assert cli.user_settings is user_settings
    assert app.user_settings is user_settings


def test_user_settings_load_once_and_return_independent_snapshots(tmp_path):
    with patch("app_settings._load_settings_file", return_value={"ui_mode": "advanced"}) as read:
        settings = UserSettings(path=tmp_path / "settings.json")
        settings.path.touch()
        read.assert_not_called()
        assert settings.ui_mode == "advanced"
        snapshot = settings.load()
        snapshot["ui_mode"] = "minimal"
        assert settings.ui_mode == "advanced"
        assert read.call_count == 1


def test_user_settings_initialization_preserves_shared_edits(tmp_path):
    template = tmp_path / "defaults.json"
    template.write_text(json.dumps({
        "ui_mode": "minimal",
        "pipeline_visibility": {"velocity_analysis": True},
    }), encoding="utf-8")
    settings = UserSettings(path=tmp_path / "active.json", default_template_path=template)
    assert settings.ui_mode == "minimal"
    settings.ui_mode = "advanced"
    assert settings.initialize_from_defaults() is True
    assert settings.initialize_from_defaults() is False
    assert settings.ui_mode == "advanced"
    assert UserSettings(settings.path).ui_mode == "advanced"
    assert settings.pipeline_visibility.selected_names() == ("velocity_analysis",)


def test_user_settings_share_live_sections_and_persist_explicitly(tmp_path):
    path = tmp_path / "settings.json"
    settings = UserSettings(path=path)
    visibility = settings.pipeline_visibility
    options = settings.pipeline_options
    assert isinstance(visibility, PipelineVisibilitySettings)
    assert isinstance(options, PipelineOptionsSettings)
    visibility["velocity_analysis"] = True
    segment_options = options.setdefault("velocity_analysis", {})
    assert isinstance(segment_options, BooleanSettings)
    segment_options["segments"] = True
    settings.band_limited_harmonics = 5
    settings.number_of_radii_in_half_diagonal = 16
    settings["custom_setting"] = {"preserved": [1, 2]}
    assert not path.exists()

    # Replacing with a live view must snapshot it before clearing the group.
    options.replace(options)
    settings.save(settings)
    disk = json.loads(path.read_text(encoding="utf-8"))
    assert disk["pipeline_options"] == {"velocity_analysis": {"segments": True}}
    assert disk["pipeline_visibility"] == {"velocity_analysis": True}
    assert disk["band_limited_harmonics"] == 5
    assert disk["number_of_radii_in_half_diagonal"] == 16
    assert disk["custom_setting"] == {"preserved": [1, 2]}
    assert visibility is settings.pipeline_visibility
    assert segment_options["segments"] is True
    assert UserSettings(path).load() == disk


def test_user_settings_import_and_reload_refresh_existing_section_references(tmp_path):
    settings = UserSettings(path=tmp_path / "active.json")
    settings.save({
        "pipeline_visibility": {"velocity_analysis": True},
        "pipeline_options": {"velocity_analysis": {"segments": True}},
    })
    visibility = settings.pipeline_visibility
    options = settings.pipeline_options["velocity_analysis"]
    source = tmp_path / "import.json"
    source.write_text(json.dumps({
        "pipeline_visibility": {"velocity_analysis": False},
        "pipeline_options": {"velocity_analysis": {"segments": False}},
        "ui_mode": "advanced",
    }), encoding="utf-8")
    settings.import_file(source)
    assert visibility["velocity_analysis"] is False
    assert options["segments"] is False
    assert settings.ui_mode == "advanced"

    UserSettings(settings.path).save({"ui_mode": "minimal"})
    assert settings.ui_mode == "advanced"
    settings.reload()
    assert settings.ui_mode == "minimal"
    assert not visibility
    assert not options


def test_user_settings_reject_failed_import_and_failed_save_without_replacing_cache(tmp_path):
    settings = UserSettings(path=tmp_path / "active.json")
    original = {"ui_mode": "advanced", "band_ratio_frequency_scale_hz": 3.0}
    settings.save(original)
    source = tmp_path / "import.json"
    source.write_text('{"band_ratio_frequency_scale_hz": 0}', encoding="utf-8")
    with pytest.raises(ValueError, match="band_ratio_frequency_scale_hz"):
        settings.import_file(source)
    with patch.object(AppSettingsStore, "save", side_effect=OSError("read only")):
        with pytest.raises(OSError, match="read only"):
            settings.save({"ui_mode": "minimal"})
    assert settings.load() == original
    assert json.loads(settings.path.read_text(encoding="utf-8")) == original


def test_user_settings_validated_fields_and_sections_reject_invalid_edits(tmp_path):
    settings = UserSettings(path=tmp_path / "settings.json")
    with pytest.raises(ValueError, match="band_limited_harmonics"):
        settings.band_limited_harmonics = -1
    with pytest.raises(ValueError, match="Boolean"):
        settings.pipeline_visibility["velocity_analysis"] = "false"
    with pytest.raises(ValueError, match="Boolean"):
        settings.pipeline_options["velocity_analysis"] = {"segments": 1}
    assert not settings.path.exists()


def test_gui_controller_uses_live_settings_views(tmp_path):
    from pipeline_engine import PipelineDescriptor, PipelineOption
    from ui.controllers.pipeline_library import PipelineLibraryController

    app = SimpleNamespace(settings=UserSettings(path=tmp_path / "settings.json"))
    app.settings.pipeline_visibility.replace({"velocity_analysis": True})
    app.settings.pipeline_options.replace({"velocity_analysis": {"segments": False}})
    options = app.settings.pipeline_options["velocity_analysis"]
    options["segments"] = True
    app.settings.ui_mode = "advanced"
    app.settings.save()
    assert options["segments"] is True
    app.settings_controller = SettingsController(app)

    pipeline = PipelineDescriptor(
        name="velocity_analysis", description="Velocity", available=True,
        options=(PipelineOption("segments", "Segments", default_enabled=False),),
    )
    app.pipeline_catalog = {pipeline.name: pipeline}
    app.pipeline_option_vars = {}
    controller = PipelineLibraryController(app)
    controller.update_summary = Mock()
    controller.sync_visibility([pipeline])
    controller.sync_options([pipeline])
    controller.set_option_visibility("velocity_analysis", "segments", False)
    assert app.settings.pipeline_options["velocity_analysis"]["segments"] is False
    assert UserSettings(app.settings.path).pipeline_options == {
        "velocity_analysis": {"segments": False},
    }


class _FileDialogs:
    def __init__(self, selected_path: str) -> None:
        self.selected_path = selected_path
        self.askopenfilename_calls = 0

    def askopenfilename(self, **_options) -> str:
        self.askopenfilename_calls += 1
        return self.selected_path


class _Variable:
    def __init__(self, value: str) -> None:
        self.value = value

    def get(self) -> str:
        return self.value

    def set(self, value: str) -> None:
        self.value = value


class SettingsImportTests(unittest.TestCase):
    def test_last_input_directory_is_persisted_without_changing_other_settings(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            directory = Path(temp_dir)
            settings_path = directory / "settings.json"
            store = UserSettings(path=settings_path)
            original = {"ui_mode": "advanced", "pipeline_visibility": {"pdf_report": False}}
            store.save(original)
            store.last_input_directory = directory
            store.save()
            reopened = UserSettings(path=settings_path)
            self.assertEqual(directory, reopened.last_input_directory)
            self.assertEqual({**original, "last_input_directory": str(directory)}, reopened.load())

    def test_invalid_last_input_directories_are_ignored(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            store = UserSettings(path=root / "settings.json")
            file_path = root / "input.holo"
            file_path.touch()
            for value in (None, 123, [], "", ".", str(root / "missing"), str(file_path), "\0"):
                with self.subTest(value=value):
                    store.save({"last_input_directory": value})
                    self.assertIsNone(store.last_input_directory)

    def test_import_file_replaces_active_settings(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            source_path = root / "selected.json"
            source = {
                "pipeline_visibility": {"pdf_report": False},
                "ui_mode": "advanced",
            }
            source_path.write_text(json.dumps(source), encoding="utf-8")
            store = UserSettings(
                path=root / "active" / "settings.json",
                default_template_path=None,
            )

            store.import_file(source_path)

            self.assertEqual(source, store.load())

    def test_import_file_rejects_invalid_json_without_replacing_settings(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            source_path = root / "invalid.json"
            source_path.write_text('{"ui_mode":', encoding="utf-8")
            store = UserSettings(
                path=root / "settings.json",
                default_template_path=None,
            )
            original = {"ui_mode": "minimal"}
            store.save(original)

            with self.assertRaisesRegex(ValueError, "Invalid JSON"):
                store.import_file(source_path)

            self.assertEqual(original, store.load())

    def test_import_file_requires_a_json_object(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            source_path = root / "list.json"
            source_path.write_text("[]", encoding="utf-8")
            store = UserSettings(
                path=root / "settings.json",
                default_template_path=None,
            )

            with self.assertRaisesRegex(TypeError, "JSON object"):
                store.import_file(source_path)

            self.assertFalse(store.path.exists())

    def test_existing_settings_use_default_band_ratio_scale(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            store = UserSettings(
                path=Path(temp_dir) / "settings.json",
                default_template_path=None,
            )
            store.save({"ui_mode": "minimal"})

            self.assertEqual(1.0, store.band_ratio_frequency_scale_hz)
            self.assertFalse(store.copy_hd_passthrough)

    def test_import_rejects_invalid_hd_passthrough_setting(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            source = root / "invalid.json"
            source.write_text(json.dumps({"copy_hd_passthrough": "false"}), encoding="utf-8")
            store = UserSettings(path=root / "settings.json", default_template_path=None)
            with self.assertRaisesRegex(ValueError, "copy_hd_passthrough"):
                store.import_file(source)
            self.assertFalse(store.path.exists())

    def test_import_rejects_invalid_band_ratio_scale_without_replacing_settings(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            source_path = root / "invalid_band_scale.json"
            source_path.write_text(
                json.dumps({"band_ratio_frequency_scale_hz": 0.0}),
                encoding="utf-8",
            )
            store = UserSettings(
                path=root / "settings.json",
                default_template_path=None,
            )
            original = {"band_ratio_frequency_scale_hz": 1.0}
            store.save(original)

            with self.assertRaisesRegex(
                ValueError,
                "band_ratio_frequency_scale_hz",
            ):
                store.import_file(source_path)

            self.assertEqual(original, store.load())

    def test_choose_config_file_imports_and_refreshes_the_ui(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            source_path = root / "selected.json"
            source_path.write_text(
                json.dumps({"ui_mode": "advanced"}),
                encoding="utf-8",
            )
            store = UserSettings(
                path=root / "settings.json",
                default_template_path=None,
            )
            dialogs = SimpleNamespace(
                showerror=Mock(),
                showinfo=Mock(),
                showwarning=Mock(),
            )
            pipeline_controller = SimpleNamespace(register=Mock())
            app = SimpleNamespace(
                settings=store,
                pipeline_library_controller=pipeline_controller,
                ui_services=SimpleNamespace(
                    file_dialogs=_FileDialogs(str(source_path)),
                    dialogs=dialogs,
                ),
            )
            controller = SettingsController(app)
            controller.apply_ui_mode = Mock()

            controller.choose_config_file()

            self.assertEqual({"ui_mode": "advanced"}, store.load())
            pipeline_controller.register.assert_called_once_with()
            controller.apply_ui_mode.assert_called_once_with(
                "advanced",
                persist=False,
            )
            dialogs.showinfo.assert_called_once()
            dialogs.showerror.assert_not_called()

    def test_choose_config_file_is_blocked_during_a_pipeline_run(self) -> None:
        file_dialogs = _FileDialogs("selected.json")
        dialogs = SimpleNamespace(showwarning=Mock())
        app = SimpleNamespace(
            _pipeline_run_active=True,
            ui_services=SimpleNamespace(
                file_dialogs=file_dialogs,
                dialogs=dialogs,
            ),
        )

        SettingsController(app).choose_config_file()

        dialogs.showwarning.assert_called_once()
        self.assertEqual(0, file_dialogs.askopenfilename_calls)


if __name__ == "__main__":
    unittest.main()
