"""Tests for importing EyeFlow settings through the Tools menu action."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

from app_settings import AppSettingsStore
from ui.controllers.settings import SettingsController


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
            store = AppSettingsStore(path=settings_path)
            original = {"ui_mode": "advanced", "pipeline_visibility": {"pdf_report": False}}
            store.save(original)
            store.save_last_input_directory(directory)

            reopened = AppSettingsStore(path=settings_path)
            self.assertEqual(directory, reopened.load_last_input_directory())
            self.assertEqual({**original, "last_input_directory": str(directory)}, reopened.load())

    def test_invalid_last_input_directories_are_ignored(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            store = AppSettingsStore(path=root / "settings.json")
            file_path = root / "input.holo"
            file_path.touch()
            for value in (None, 123, [], "", ".", str(root / "missing"), str(file_path), "\0"):
                with self.subTest(value=value):
                    store.save({"last_input_directory": value})
                    self.assertIsNone(store.load_last_input_directory())

    def test_import_file_replaces_active_settings(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            source_path = root / "selected.json"
            source = {
                "pipeline_visibility": {"pdf_report": False},
                "ui_mode": "advanced",
            }
            source_path.write_text(json.dumps(source), encoding="utf-8")
            store = AppSettingsStore(
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
            store = AppSettingsStore(
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
            store = AppSettingsStore(
                path=root / "settings.json",
                default_template_path=None,
            )

            with self.assertRaisesRegex(TypeError, "JSON object"):
                store.import_file(source_path)

            self.assertFalse(store.path.exists())

    def test_existing_settings_use_default_band_ratio_scale(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            store = AppSettingsStore(
                path=Path(temp_dir) / "settings.json",
                default_template_path=None,
            )
            store.save({"ui_mode": "minimal"})

            self.assertEqual(1.0, store.load_band_ratio_frequency_scale_hz())
            self.assertFalse(store.load_copy_hd_passthrough())

    def test_import_rejects_invalid_hd_passthrough_setting(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            source = root / "invalid.json"
            source.write_text(json.dumps({"copy_hd_passthrough": "false"}), encoding="utf-8")
            store = AppSettingsStore(path=root / "settings.json", default_template_path=None)
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
            store = AppSettingsStore(
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
            store = AppSettingsStore(
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
                settings_store=store,
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
