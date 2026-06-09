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

from input_output import resolve_selected_run_layouts  # noqa: E402
from ui.input_state import InputStateMixin  # noqa: E402


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


class _InputStateHarness(InputStateMixin):
    def __init__(self) -> None:
        self.holo_hd_status_var = _Var()
        self.holo_dv_status_var = _Var()
        self._muted_fg = "muted"
        self._error_color = "error"
        self._success_color = "success"


class StemListInputTests(unittest.TestCase):
    def test_resolves_txt_stems_directly(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir)
            first_hd, first_dv = _write_companion_h5_files(root, "subject_a")
            second_hd, second_dv = _write_companion_h5_files(root, "subject_b")
            input_list = root / "inputs.txt"
            input_list.write_text("subject_a\nsubject_b\n", encoding="utf-8")

            resolved = resolve_selected_run_layouts([input_list])

            self.assertEqual(
                ["subject_a", "subject_b"],
                [item.stem for item in resolved],
            )
            self.assertEqual(
                [root / "subject_a", root / "subject_b"],
                [item.root_dir for item in resolved],
            )
            self.assertEqual([first_hd, second_hd], [item.hd_h5 for item in resolved])
            self.assertEqual([first_dv, second_dv], [item.dv_h5 for item in resolved])

    def test_ignores_blank_txt_lines(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir)
            hd_h5, dv_h5 = _write_companion_h5_files(root, "subject_a")
            input_list = root / "inputs.txt"
            input_list.write_text("\nsubject_a\n\n", encoding="utf-8")

            resolved = resolve_selected_run_layouts([input_list])

            self.assertEqual(["subject_a"], [item.stem for item in resolved])
            self.assertEqual([hd_h5], [item.hd_h5 for item in resolved])
            self.assertEqual([dv_h5], [item.dv_h5 for item in resolved])

    def test_gui_status_counts_txt_stem_list_hd_dv_inputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir)
            _write_companion_h5_files(root, "subject_a")
            _write_hd_h5_file(root, "subject_b")
            input_list = root / "inputs.txt"
            input_list.write_text("subject_a\nsubject_b\n", encoding="utf-8")
            app = _InputStateHarness()

            app._update_minimal_found_statuses([input_list])

            self.assertEqual("HD 2/2 found", app.holo_hd_status_var.get())
            self.assertEqual(
                "DV 1/2 found: missing subject_b",
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
            input_list.write_text("subject_a\n", encoding="utf-8")
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
