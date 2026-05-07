from __future__ import annotations

import tkinter as tk
from collections.abc import Sequence
from pathlib import Path
from tkinter import messagebox

from input_output import HOLO_SUFFIX

from .compat import DND_FILES


class DragDropMixin:
    def _install_drop_targets(self) -> None:
        if DND_FILES is None:
            return
        for widget in (
            self,
            self.minimal_view,
            self.advanced_view,
            self.run_tab,
            self.pipeline_library_tab,
        ):
            self._register_drop_target(widget)

        slot_widgets = {
            "holo": (
                getattr(self, "minimal_holo_browse_button", None),
                getattr(self, "minimal_holo_input_path_label", None),
                getattr(self, "holo_input_entry", None),
                getattr(self, "holo_browse_button", None),
            ),
        }
        for slot, widgets in slot_widgets.items():
            for widget in widgets:
                if widget is not None:
                    self._register_drop_target(widget, slot)

    def _register_drop_target(
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
                lambda event, target_slot=slot_hint: self._on_input_drop(
                    event, slot_hint=target_slot
                ),
            )
        except (AttributeError, tk.TclError):
            pass

    def _handle_dropped_paths(
        self,
        dropped_paths: Sequence[Path],
        *,
        slot_hint: str | None = None,
    ) -> bool:
        del slot_hint
        valid_paths: list[Path] = []
        for dropped_path in dropped_paths:
            cleaned = str(dropped_path).strip().strip("{}")
            if not cleaned:
                continue
            candidate = Path(cleaned).expanduser()
            if candidate.is_file() and candidate.suffix.lower() == HOLO_SUFFIX:
                valid_paths.append(candidate)

        if not valid_paths:
            return False

        if len(valid_paths) != 1:
            return False

        self._assign_holo_input_path(valid_paths[0])
        self._log_run(f"[INPUT] Drag and drop HOLO -> {valid_paths[0]}")
        return True

    def _on_input_drop(self, event, *, slot_hint: str | None = None) -> None:
        raw_data = getattr(event, "data", "")
        try:
            dropped_values = self.tk.splitlist(raw_data)
        except tk.TclError:
            dropped_values = (raw_data,)

        dropped_paths = [Path(value) for value in dropped_values if value]
        if self._handle_dropped_paths(dropped_paths, slot_hint=slot_hint):
            return

        messagebox.showwarning(
            "Unsupported drop",
            "Drop one .holo file into the input area.",
        )
