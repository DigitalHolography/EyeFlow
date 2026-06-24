"""Controller for selected input paths and drag-and-drop."""

from __future__ import annotations

import tkinter as tk
from collections.abc import Sequence
from pathlib import Path

from input_output import (
    INPUT_LIST_SUFFIX,
    holo_input_status,
    read_holo_input_list,
    stem_input_status,
)
from input_output.output_manager import OutputManager

from ..services import services_for

try:
    from tkinterdnd2 import DND_FILES
except ImportError:  # optional dependency
    DND_FILES = None

HOLO_SUFFIX = ".holo"


class InputController:
    def __init__(self, app) -> None:
        self.app = app

    def on_holo_input_changed(self, *_args) -> None:
        if not self.app._synchronizing_holo_input_var:
            self.app._selected_holo_input_paths = []
        self.update_minimal_path_labels()

    def path_from_var(self, raw_value: str) -> Path | None:
        value = raw_value.strip()
        if not value:
            return None
        path = Path(value).expanduser()
        return path if path.is_absolute() else Path.cwd() / path

    def set_holo_input_var(self, value: str) -> None:
        self.app._synchronizing_holo_input_var = True
        try:
            self.app.holo_input_var.set(value)
        finally:
            self.app._synchronizing_holo_input_var = False

    def selected_holo_paths(self) -> list[Path]:
        if self.app._selected_holo_input_paths:
            return list(self.app._selected_holo_input_paths)
        selected_path = self.selected_holo_path()
        return [selected_path] if selected_path is not None else []

    def selected_holo_path(self) -> Path | None:
        if self.app._selected_holo_input_paths:
            return self.app._selected_holo_input_paths[0]
        return self.path_from_var(self.app.holo_input_var.get() or "")

    def assign_holo_input_path(self, input_path: Path) -> None:
        normalized_path = self._normalize_path(input_path)
        self.app._selected_holo_input_paths = [normalized_path]
        self.set_holo_input_var(str(normalized_path))
        self.apply_input_defaults(normalized_path)

    def assign_holo_input_paths(self, input_paths: Sequence[Path]) -> None:
        normalized_paths = [self._normalize_path(path) for path in input_paths]
        if not normalized_paths:
            self.app._selected_holo_input_paths = []
            self.set_holo_input_var("")
            self.apply_input_defaults(None)
            return

        self.app._selected_holo_input_paths = normalized_paths
        display_value = str(normalized_paths[0])
        if len(normalized_paths) > 1:
            display_value += f" (+{len(normalized_paths) - 1} more)"
        self.set_holo_input_var(display_value)
        self.apply_input_defaults(normalized_paths[0])

    def reference_holo_tooltip_text(self) -> str:
        return "Pick one or more reference .holo files, or a .txt path list."

    def install_drop_targets(self) -> None:
        if DND_FILES is None:
            return
        for widget in (
            self.app,
            self.app.minimal_view,
            self.app.advanced_view,
            self.app.run_tab,
            self.app.pipeline_library_tab,
        ):
            self.register_drop_target(widget)

        slot_widgets = (
            getattr(self.app, "minimal_holo_browse_button", None),
            getattr(self.app, "minimal_holo_input_path_label", None),
            getattr(self.app, "holo_input_entry", None),
            getattr(self.app, "holo_browse_button", None),
        )
        for widget in slot_widgets:
            if widget is not None:
                self.register_drop_target(widget, "holo")

    def register_drop_target(
        self,
        widget: tk.Misc,
        slot_hint: str | None = None,
    ) -> None:
        if DND_FILES is None:
            return
        try:
            widget.drop_target_register(DND_FILES)
            widget.dnd_bind(
                "<<Drop>>",
                lambda event, target_slot=slot_hint: self.on_input_drop(
                    event, slot_hint=target_slot
                ),
            )
        except (AttributeError, tk.TclError):
            pass

    def handle_dropped_paths(
        self,
        dropped_paths: Sequence[Path],
        *,
        slot_hint: str | None = None,
    ) -> bool:
        del slot_hint
        valid_paths = [path for path in self._clean_dropped_paths(dropped_paths)]
        if not valid_paths:
            return False

        if len(valid_paths) == 1:
            self.assign_holo_input_path(valid_paths[0])
            self.app.progress_controller.log_run(
                f"[INPUT] Drag and drop HOLO -> {valid_paths[0]}"
            )
            return True

        self.assign_holo_input_paths(valid_paths)
        self.app.progress_controller.log_run(
            f"[INPUT] Drag and drop HOLO -> {valid_paths[0]} "
            f"(+{len(valid_paths) - 1} more)"
        )
        return True

    def on_input_drop(self, event, *, slot_hint: str | None = None) -> None:
        raw_data = getattr(event, "data", "")
        try:
            dropped_values = self.app.tk.splitlist(raw_data)
        except tk.TclError:
            dropped_values = (raw_data,)

        dropped_paths = [Path(value) for value in dropped_values if value]
        if self.handle_dropped_paths(dropped_paths, slot_hint=slot_hint):
            return

        services_for(self.app).dialogs.showwarning(
            "Unsupported drop",
            "Drop one or more .holo files, or one .txt path list.",
        )

    def update_minimal_path_labels(self) -> None:
        holo_paths = self.selected_holo_paths()
        if not holo_paths:
            self.app.minimal_holo_input_path_var.set("No input selected")
        elif len(holo_paths) == 1:
            self.app.minimal_holo_input_path_var.set(str(holo_paths[0]))
        else:
            self.app.minimal_holo_input_path_var.set(
                f"{holo_paths[0]} (+{len(holo_paths) - 1} more)"
            )
        self.update_holo_input_found_statuses(holo_paths)

    def update_holo_input_found_statuses(self, holo_paths: Sequence[Path]) -> None:
        if not holo_paths:
            self.set_shared_holo_status(
                hd_text="",
                hd_color=self.app._muted_fg,
                dv_text="",
                dv_color=self.app._muted_fg,
            )
            return
        if len(holo_paths) == 1 and holo_paths[0].suffix.lower() == INPUT_LIST_SUFFIX:
            self.update_input_list_found_statuses(holo_paths[0])
            return

        hd_found_count, dv_found_count, missing_hd, missing_dv = (
            self._count_holo_statuses(holo_paths)
        )
        total_files = len(holo_paths)
        self.set_shared_holo_status(
            hd_text=found_status_text(
                "HD", hd_found_count, total_files, missing_hd
            ),
            hd_color=self.app._error_color if missing_hd else self.app._success_color,
            dv_text=found_status_text(
                "DV", dv_found_count, total_files, missing_dv
            ),
            dv_color=self.app._error_color if missing_dv else self.app._success_color,
        )

    def update_input_list_found_statuses(self, input_list_path: Path) -> None:
        try:
            input_list = read_holo_input_list(input_list_path)
        except (OSError, ValueError) as exc:
            self.set_shared_holo_status(
                hd_text=str(exc),
                hd_color=self.app._error_color,
                dv_text="DV unavailable",
                dv_color=self.app._error_color,
            )
            return

        inputs = list(input_list.path_stem_pairs)
        stems = [stem for _root_dir, stem in inputs]
        statuses = [stem_input_status(stem, root_dir) for root_dir, stem in inputs]
        missing_hd = [stem for stem, status in zip(stems, statuses) if not status.hd]
        missing_dv = [stem for stem, status in zip(stems, statuses) if not status.dv]
        self.set_shared_holo_status(
            hd_text=found_status_text(
                "HD", len(stems) - len(missing_hd), len(stems), missing_hd
            ),
            hd_color=self.app._error_color if missing_hd else self.app._success_color,
            dv_text=found_status_text(
                "DV", len(stems) - len(missing_dv), len(stems), missing_dv
            ),
            dv_color=self.app._error_color if missing_dv else self.app._success_color,
        )

    def set_shared_holo_status(
        self,
        *,
        hd_text: str,
        hd_color: str,
        dv_text: str,
        dv_color: str,
    ) -> None:
        self.app.holo_hd_status_var.set(hd_text)
        self.app.holo_dv_status_var.set(dv_text)
        self.set_holo_status_rows_visible(bool(hd_text or dv_text))
        for label_name, color in (
            ("minimal_holo_hd_status_label", hd_color),
            ("minimal_holo_dv_status_label", dv_color),
            ("holo_hd_status_label", hd_color),
            ("holo_dv_status_label", dv_color),
        ):
            label = getattr(self.app, label_name, None)
            if label is not None:
                label.configure(fg=color)

    def set_holo_status_rows_visible(self, visible: bool) -> None:
        for frame_name in ("minimal_holo_status_frame", "holo_status_frame"):
            frame = getattr(self.app, frame_name, None)
            if frame is None:
                continue
            if visible:
                frame.grid()
            else:
                frame.grid_remove()

    def prepare_output_manager_for_input(self, input_path: Path) -> OutputManager:
        output_manager = OutputManager.from_holo(input_path)
        output_manager.prepare(replace=True)
        return output_manager

    def apply_input_defaults(self, input_path: Path | None) -> None:
        del input_path
        self.app.progress_controller.reset_progress()
        self.app.progress_controller.set_minimal_status("Ready.")

    def _normalize_path(self, input_path: Path) -> Path:
        path = input_path.expanduser()
        return path if path.is_absolute() else Path.cwd() / path

    def _clean_dropped_paths(self, dropped_paths: Sequence[Path]) -> list[Path]:
        valid_paths: list[Path] = []
        for dropped_path in dropped_paths:
            cleaned = str(dropped_path).strip().strip("{}")
            if not cleaned:
                continue
            candidate = Path(cleaned).expanduser()
            if candidate.is_file() and candidate.suffix.lower() in {
                HOLO_SUFFIX,
                INPUT_LIST_SUFFIX,
            }:
                valid_paths.append(candidate)
        return valid_paths

    def _count_holo_statuses(
        self,
        holo_paths: Sequence[Path],
    ) -> tuple[int, int, list[str], list[str]]:
        hd_found_count = 0
        dv_found_count = 0
        missing_hd: list[str] = []
        missing_dv: list[str] = []
        for holo_path in holo_paths:
            normalized_holo = self._normalize_path(holo_path)
            status = holo_input_status(normalized_holo, require_holo_file=True)
            if status.hd:
                hd_found_count += 1
            else:
                missing_hd.append(normalized_holo.stem)
            if status.dv:
                dv_found_count += 1
            else:
                missing_dv.append(normalized_holo.stem)
        return hd_found_count, dv_found_count, missing_hd, missing_dv


def found_status_text(
    label: str,
    found_count: int,
    total_count: int,
    missing_stems: Sequence[str],
) -> str:
    if total_count == 1:
        return f"{label} found" if found_count else f"{label} not found"
    text = f"Found {found_count}/{total_count} {label}"
    if missing_stems:
        text += ": missing " + ", ".join(missing_stems)
    return text
