"""Build the main Tk views used by EyeFlow."""

from __future__ import annotations

import tkinter as tk
from tkinter import ttk

from .advanced import build_advanced_view
from .minimal import build_minimal_view


class AppViewBuilder:
    def __init__(self, app) -> None:
        self.app = app

    def build_ui(self) -> None:
        self._build_menu()

        container = ttk.Frame(self.app, padding=0)
        container.pack(fill="both", expand=True)
        self.app.main_container = container

        self.app.minimal_view = ttk.Frame(container, padding=0)
        self.app.advanced_view = ttk.Frame(container, padding=(20, 20, 20, 20))

        build_minimal_view(self.app, self.app.minimal_view)
        build_advanced_view(self.app, self.app.advanced_view)

    def _build_menu(self) -> None:
        self.app.ui_mode_var = tk.StringVar(value=self.app.ui_mode)
        menu_bar = tk.Menu(self.app)
        view_menu = tk.Menu(menu_bar, tearoff=False)
        view_menu.add_radiobutton(
            label="Minimal UI",
            value="minimal",
            variable=self.app.ui_mode_var,
            command=lambda: self.app.settings_controller.apply_ui_mode(
                self.app.ui_mode_var.get()
            ),
        )
        view_menu.add_radiobutton(
            label="Advanced UI",
            value="advanced",
            variable=self.app.ui_mode_var,
            command=lambda: self.app.settings_controller.apply_ui_mode(
                self.app.ui_mode_var.get()
            ),
        )
        menu_bar.add_cascade(label="View", menu=view_menu)
        help_menu = tk.Menu(menu_bar, tearoff=False)
        help_menu.add_command(
            label="Open Log File",
            command=self.app.progress_controller.open_run_log_file,
        )
        menu_bar.add_cascade(label="Help", menu=help_menu)
        self.app.configure(menu=menu_bar)
