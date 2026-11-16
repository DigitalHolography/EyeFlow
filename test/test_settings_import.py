"""Tests for importing EyeFlow settings through the Tools menu action."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock
import pytest

from app_settings import AppSettingsStore
from ui.controllers.settings import SettingsController
from app_settings import validate_annulus_settings
from calculations.topology import OpticDisc
from pipelines.topology_core.annulus_sampling import resolve_annulus_geometry


@pytest.mark.parametrize("settings,size,stride", [
    ({}, 4.0, 2.0),
    ({"annulus_size": 7.5}, 7.5, 3.75),
    ({"annulus_stride": 1.25}, 4.0, 1.25),
    ({"number_of_radii_in_FOV": 10}, 10.0, 10.0),
    ({"NumberOfRadiiOverFOV": 20}, 5.0, 5.0),
    ({"annulus_size": 7.5, "number_of_radii_in_FOV": -1}, 7.5, 3.75),
    ({"annulus_stride": 1.25, "number_of_radii_in_FOV": 10}, 4.0, 1.25),
])
def test_saved_annulus_settings_and_count_migration(tmp_path, settings, size, stride):
    path = tmp_path / "import.json"
    path.write_text(json.dumps(settings), encoding="utf-8")
    store = AppSettingsStore(path=tmp_path / "active.json")
    store.import_file(path)
    ctx = SimpleNamespace(annulus_settings=store.load_annulus_settings())
    # Acquisition counts are not intentional user settings and must be ignored.
    ctx.attrs = {"number_of_radii_in_FOV": 2}
    geometry = resolve_annulus_geometry(ctx, OpticDisc(None, (50, 40), 4, 4), (80, 100))
    metadata = geometry.metadata((80, 100))
    assert metadata["annulus_size_pixels"] == pytest.approx(size)
    assert metadata["annulus_stride_pixels"] == pytest.approx(stride)


@pytest.mark.parametrize("value", [0, -0.5, float("nan"), float("inf"), True, "4"])
@pytest.mark.parametrize("name", ["annulus_size", "annulus_stride"])
def test_invalid_annulus_settings_do_not_replace_saved_configuration(tmp_path, name, value):
    store = AppSettingsStore(path=tmp_path / "active.json")
    store.save({"annulus_size": 5})
    source = tmp_path / "invalid.json"
    source.write_text(json.dumps({name: value}), encoding="utf-8")
    with pytest.raises(ValueError, match=name):
        store.import_file(source)
    assert store.load() == {"annulus_size": 5}


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
