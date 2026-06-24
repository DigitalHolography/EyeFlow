"""Controller for settings persistence and window mode behavior."""

from __future__ import annotations

import tkinter as tk

from ..services import services_for


class SettingsController:
    def __init__(self, app) -> None:
        self.app = app

    def ensure_default_settings(self) -> None:
        try:
            self.app.settings_store.initialize_from_defaults()
        except OSError as exc:
            self.show_settings_warning(
                "Settings not initialized",
                f"Could not create default settings file:\n{exc}",
            )

    def show_settings_warning(self, title: str, details: str) -> None:
        if self.app._settings_warning_shown:
            return
        self.app._settings_warning_shown = True
        services_for(self.app).dialogs.showwarning(title, details)

    def persist_ui_mode(self) -> None:
        try:
            self.app.settings_store.save_ui_mode(self.app.ui_mode)
        except OSError as exc:
            self.show_settings_warning(
                "Settings not saved",
                f"Could not save UI mode preference:\n{exc}",
            )

    def persist_trim_h5source(self) -> None:
        try:
            self.app.settings_store.save_trim_h5source(
                self.app._trim_h5source.get()
            )
        except OSError as exc:
            self.show_settings_warning(
                "Settings not saved",
                f"Could not save trim preference:\n{exc}",
            )

    def window_size_for_mode(self, mode: str) -> tuple[int, int, int, int]:
        screen_width = self.app.winfo_screenwidth()
        screen_height = self.app.winfo_screenheight()
        if mode == "advanced":
            width = min(900, max(760, screen_width - 240), screen_width)
            height = min(520, max(520, screen_height - 240), screen_height)
            min_width = min(620, width)
            min_height = min(420, height)
        else:
            width = max(560, min(660, screen_width - 260))
            height = max(460, min(500, screen_height - 260))
            min_width = min(500, width)
            min_height = min(460, height)
        return width, height, min_width, min_height

    def ensure_window_size_for_mode(
        self,
        mode: str,
        *,
        force_target_size: bool = False,
    ) -> None:
        target_width, target_height, min_width, min_height = (
            self.window_size_for_mode(mode)
        )
        if mode == "minimal":
            self.app.minimal_content.update_idletasks()
            min_width = max(min_width, self.app.minimal_content.winfo_reqwidth())
            min_height = max(min_height, self.app.minimal_content.winfo_reqheight())
        self.app.minsize(min_width, min_height)

        width, height = self._target_window_size(
            mode,
            target_width,
            target_height,
            min_width,
            min_height,
            force_target_size=force_target_size,
        )
        if width is None or height is None:
            return

        screen_width = self.app.winfo_screenwidth()
        screen_height = self.app.winfo_screenheight()
        x = max(min(self.app.winfo_x(), screen_width - width), 0)
        y = max(min(self.app.winfo_y(), screen_height - height), 0)
        self.app.geometry(f"{width}x{height}+{x}+{y}")

    def apply_ui_mode(self, mode: str, *, persist: bool = True) -> None:
        normalized_mode = "advanced" if mode == "advanced" else "minimal"
        self.app.ui_mode = normalized_mode
        self.app.ui_mode_var.set(normalized_mode)

        self.app.minimal_view.pack_forget()
        self.app.advanced_view.pack_forget()
        if normalized_mode == "advanced":
            self.app.advanced_view.pack(fill="both", expand=True)
        else:
            self.app.minimal_view.pack(fill="both", expand=True)

        self.app.update_idletasks()
        self.ensure_window_size_for_mode(
            normalized_mode,
            force_target_size=(normalized_mode == "minimal"),
        )
        if persist:
            self.persist_ui_mode()

    def on_close(self) -> None:
        if getattr(self.app, "_pipeline_run_active", False):
            services_for(self.app).dialogs.showwarning(
                "Run in progress",
                "Wait for the current pipeline run to finish before closing EyeFlow.",
            )
            return

        self.persist_ui_mode()
        self.persist_trim_h5source()
        self.app.destroy()

    def _target_window_size(
        self,
        mode: str,
        target_width: int,
        target_height: int,
        min_width: int,
        min_height: int,
        *,
        force_target_size: bool,
    ) -> tuple[int | None, int | None]:
        current_width = max(self.app.winfo_width(), 1)
        current_height = max(self.app.winfo_height(), 1)
        if (
            not force_target_size
            and current_width >= min_width
            and current_height >= min_height
        ):
            return None, None

        screen_width = self.app.winfo_screenwidth()
        screen_height = self.app.winfo_screenheight()
        if not force_target_size:
            return (
                min(max(current_width, min_width), screen_width),
                min(max(current_height, min_height), screen_height),
            )
        if mode == "minimal":
            try:
                if self.app.state() != "normal":
                    self.app.state("normal")
            except tk.TclError:
                pass
            self.app.minimal_content.update_idletasks()
            return min(min_width, screen_width), min(min_height, screen_height)
        return min(target_width, screen_width), min(target_height, screen_height)
