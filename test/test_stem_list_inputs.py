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
from input_output.runs.inputs import holo_input_status, stem_input_status  # noqa: E402


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


class HoloInputListTests(unittest.TestCase):
    def test_quoted_list_resolves_each_companion_under_its_own_directory(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir)
            first_dir, second_dir = root / "first drive", root / "second drive"
            first_dir.mkdir()
            second_dir.mkdir()
            holos = [_write_holo_file(parent, "same_stem") for parent in (first_dir, second_dir)]
            companions = [_write_companion_h5_files(parent, "same_stem")
                          for parent in (first_dir, second_dir)]
            input_list = root / "inputs.txt"
            input_list.write_text(f'  "{holos[0]}"  \n\n\'{holos[1]}\'\n', encoding="utf-8")
            resolved = resolve_selected_run_layouts([input_list])
            self.assertEqual(holos, [item.holo_path for item in resolved])
            self.assertEqual([pair[0] for pair in companions], [item.hd_h5 for item in resolved])
            self.assertEqual([pair[1] for pair in companions], [item.dv_h5 for item in resolved])

    def test_status_requires_an_unambiguous_h5_selection(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir)
            stem = "subject_a"
            holo = _write_holo_file(root, stem)
            hd_dir = root / stem / f"{stem}_HD" / "h5"
            hd_dir.mkdir(parents=True)
            for name in ("one.h5", "two.h5"):
                with h5py.File(hd_dir / name, "w"):
                    pass

            self.assertFalse(holo_input_status(holo).hd)
            self.assertFalse(stem_input_status(stem, root).hd)

            with h5py.File(hd_dir / f"{stem}_HD_output.h5", "w"):
                pass
            self.assertTrue(holo_input_status(holo).hd)
            self.assertTrue(stem_input_status(stem, root).hd)

    def test_mixed_txt_and_holo_selection_raises_value_error(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir)
            input_list = root / "inputs.txt"
            holo_path = root / "subject_a.holo"
            input_list.write_text(f"{holo_path}\n", encoding="utf-8")
            holo_path.write_text("holo", encoding="utf-8")

            with self.assertRaises(ValueError):
                resolve_selected_run_layouts([input_list, holo_path])


if __name__ == "__main__":
    unittest.main()
