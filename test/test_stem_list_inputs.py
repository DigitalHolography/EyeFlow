"""Tests for resolving direct stem-list inputs."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

import h5py

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from input_output import read_holo_input_list, resolve_selected_run_layouts  # noqa: E402
from ui.controllers.input import InputController  # noqa: E402


def _write_holo_file(root: Path, stem: str) -> Path:
    holo_path = root / f"{stem}.holo"
    holo_path.write_text("holo", encoding="utf-8")
    return holo_path


def _write_companion_h5_files(root: Path, stem: str) -> tuple[Path, Path]:
    hd_h5 = root / stem / f"{stem}_HD" / "h5" / f"{stem}_HD_output.h5"
    dv_h5 = root / stem / f"{stem}_DV" / "h5" / f"{stem}_DV.h5"
    hd_h5.parent.mkdir(parents=True)
    dv_h5.parent.mkdir(parents=True)
    with h5py.File(hd_h5, "w"):
        pass
    with h5py.File(dv_h5, "w"):
        pass
    return hd_h5, dv_h5


def _write_hd_h5_file(root: Path, stem: str) -> Path:
    hd_h5 = root / stem / f"{stem}_HD" / "h5" / f"{stem}_HD_output.h5"
    hd_h5.parent.mkdir(parents=True)
    with h5py.File(hd_h5, "w"):
        pass
    return hd_h5


class _Var:
    def __init__(self, value: str = "") -> None:
        self.value = value

    def get(self) -> str:
        return self.value

    def set(self, value: str) -> None:
        self.value = value


class _InputStateHarness:
    def __init__(self) -> None:
        self.holo_hd_status_var = _Var()
        self.holo_dv_status_var = _Var()
        self._muted_fg = "muted"
        self._error_color = "error"
        self._success_color = "success"
        self.input_controller = InputController(self)


class HoloInputListTests(unittest.TestCase):
    def test_reads_holo_paths_as_root_dir_stem_pairs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir)
            first_dir = root / "first"
            second_dir = root / "second"
            first_dir.mkdir()
            second_dir.mkdir()
            first_holo = _write_holo_file(first_dir, "subject_a")
            second_holo = _write_holo_file(second_dir, "subject_b")
            input_list = root / "inputs.txt"
            input_list.write_text(f"{first_holo}\n{second_holo}\n", encoding="utf-8")

            parsed = read_holo_input_list(input_list)

            self.assertEqual(
                ((first_dir, "subject_a"), (second_dir, "subject_b")),
                parsed.path_stem_pairs,
            )

    def test_resolves_txt_holo_paths_directly(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir)
            first_dir = root / "first"
            second_dir = root / "second"
            first_dir.mkdir()
            second_dir.mkdir()
            first_holo = _write_holo_file(first_dir, "subject_a")
            second_holo = _write_holo_file(second_dir, "subject_b")
            first_hd, first_dv = _write_companion_h5_files(first_dir, "subject_a")
            second_hd, second_dv = _write_companion_h5_files(second_dir, "subject_b")
            input_list = root / "inputs.txt"
            input_list.write_text(f"{first_holo}\n{second_holo}\n", encoding="utf-8")

            resolved = resolve_selected_run_layouts([input_list])

            self.assertEqual(
                ["subject_a", "subject_b"],
                [item.stem for item in resolved],
            )
            self.assertEqual(
                [first_dir / "subject_a", second_dir / "subject_b"],
                [item.root_dir for item in resolved],
            )
            self.assertEqual(
                [first_holo, second_holo],
                [item.holo_path for item in resolved],
            )
            self.assertEqual([first_hd, second_hd], [item.hd_h5 for item in resolved])
            self.assertEqual([first_dv, second_dv], [item.dv_h5 for item in resolved])

    def test_ignores_blank_txt_lines(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir)
            holo_path = _write_holo_file(root, "subject_a")
            hd_h5, dv_h5 = _write_companion_h5_files(root, "subject_a")
            input_list = root / "inputs.txt"
            input_list.write_text(f"\n{holo_path}\n\n", encoding="utf-8")

            resolved = resolve_selected_run_layouts([input_list])

            self.assertEqual(["subject_a"], [item.stem for item in resolved])
            self.assertEqual([hd_h5], [item.hd_h5 for item in resolved])
            self.assertEqual([dv_h5], [item.dv_h5 for item in resolved])

    def test_gui_status_counts_txt_holo_path_list_hd_dv_inputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir)
            first_holo = _write_holo_file(root, "subject_a")
            second_holo = _write_holo_file(root, "subject_b")
            _write_companion_h5_files(root, "subject_a")
            _write_hd_h5_file(root, "subject_b")
            input_list = root / "inputs.txt"
            input_list.write_text(f"{first_holo}\n{second_holo}\n", encoding="utf-8")
            app = _InputStateHarness()

            app.input_controller.update_holo_input_found_statuses([input_list])

            self.assertEqual("Found 2/2 HD", app.holo_hd_status_var.get())
            self.assertEqual(
                "Found 1/2 DV: missing subject_b",
                app.holo_dv_status_var.get(),
            )

    def test_empty_txt_raises_value_error(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            input_list = Path(tmp_dir) / "inputs.txt"
            input_list.write_text("\n\n", encoding="utf-8")

            with self.assertRaises(ValueError):
                resolve_selected_run_layouts([input_list])

    def test_mixed_txt_and_holo_selection_raises_value_error(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir)
            input_list = root / "inputs.txt"
            holo_path = root / "subject_a.holo"
            input_list.write_text(f"{holo_path}\n", encoding="utf-8")
            holo_path.write_text("holo", encoding="utf-8")

            with self.assertRaises(ValueError):
                resolve_selected_run_layouts([input_list, holo_path])

    def test_existing_holo_selection_still_resolves(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir)
            hd_h5, dv_h5 = _write_companion_h5_files(root, "subject_a")
            holo_path = root / "subject_a.holo"
            holo_path.write_text("holo", encoding="utf-8")

            resolved = resolve_selected_run_layouts([holo_path])

            self.assertEqual(["subject_a"], [item.stem for item in resolved])
            self.assertEqual([hd_h5], [item.hd_h5 for item in resolved])
            self.assertEqual([dv_h5], [item.dv_h5 for item in resolved])


if __name__ == "__main__":
    unittest.main()
