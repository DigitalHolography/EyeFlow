"""Controller for window setup, theme, and view construction."""

from __future__ import annotations

import tkinter as tk
import tkinter.font as tkfont
from tkinter import ttk

from ..views import AppViewBuilder

try:
    import sv_ttk
except ImportError:  # optional dependency
    sv_ttk = None


class ViewController:
    def __init__(self, app) -> None:
        self.app = app
        self.builder = AppViewBuilder(app)

    def set_initial_window_size(self) -> None:
        width, height, min_width, min_height = (
            self.app.settings_controller.window_size_for_mode(self.app.ui_mode)
        )
        screen_width = self.app.winfo_screenwidth()
        screen_height = self.app.winfo_screenheight()
        width = min(width, screen_width)
        height = min(height, screen_height)
        x = max((screen_width - width) // 2, 0)
        y = max((screen_height - height) // 2, 0)
        self.app.geometry(f"{width}x{height}+{x}+{y}")
        self.app.minsize(min_width, min_height)

    def apply_theme(self) -> None:
        style = ttk.Style(self.app)
        self.app._style = style
        if sv_ttk:
            try:
                sv_ttk.set_theme("dark")
            except Exception:
                pass

        fallback_bg = "#0f1116"
        fallback_surface = "#1b1f27"
        fallback_fg = "#e8eef5"
        fallback_muted = "#9aa6b5"
        fallback_accent = "#4f9dff"

        bg = style.lookup("TFrame", "background") or fallback_bg
        fg = style.lookup("TLabel", "foreground") or fallback_fg
        surface = (
            style.lookup("TEntry", "fieldbackground")
            or style.lookup("TEntry", "background")
            or fallback_surface
        )
        muted = (
            style.lookup("TLabel", "foreground", state=("disabled",))
            or fallback_muted
        )
        accent = (
            style.lookup("TButton", "bordercolor")
            or style.lookup("TNotebook", "foreground")
            or fallback_accent
        )

        self.app.configure(bg=bg)
        self.app._text_bg = surface
        self.app._text_fg = fg
        self.app._muted_fg = muted
        self.app._bg_color = bg
        self.app._surface_color = surface
        self.app._accent_color = accent
        self.app._success_color = "#3fb37f"
        self.app._error_color = "#ff6b6b"
        self.configure_progress_styles()
        self.configure_velocity_estimator_styles()

    def configure_progress_styles(self) -> None:
        progress_colors = {
            self.app._progress_primary_style: self.app._accent_color,
            self.app._progress_final_style: "#3fb37f",
        }
        for style_name, color in progress_colors.items():
            try:
                self.app._style.configure(
                    style_name,
                    troughcolor=self.app._surface_color,
                    background=color,
                    bordercolor=color,
                    lightcolor=color,
                    darkcolor=color,
                )
            except tk.TclError:
                self.app._style.configure(style_name, background=color)

    def configure_velocity_estimator_styles(self) -> None:
        heading_font = tkfont.nametofont("TkDefaultFont").copy()
        heading_font.configure(weight="bold")
        self.app._velocity_estimator_heading_font = heading_font
        self.app._style.configure(
            "EstimatorHeading.TLabel",
            font=heading_font,
            foreground=self.app._accent_color,
        )
        self.app._style.configure(
            "EstimatorToggle.TFrame",
            background=self.app._surface_color,
            bordercolor=self.app._accent_color,
            borderwidth=1,
            relief="solid",
        )
        self.app._style.configure(
            "EstimatorToggle.Toolbutton",
            background=self.app._surface_color,
            foreground=self.app._text_fg,
            bordercolor=self.app._surface_color,
            lightcolor=self.app._surface_color,
            darkcolor=self.app._surface_color,
            focuscolor=self.app._surface_color,
            padding=(10, 5),
            relief="flat",
        )
        self.app._style.map(
            "EstimatorToggle.Toolbutton",
            background=[
                ("selected", self.app._accent_color),
                ("pressed", self.app._accent_color),
                ("active", self.app._accent_color),
                ("!selected", self.app._surface_color),
            ],
            foreground=[
                ("selected", self.app._bg_color),
                ("pressed", self.app._bg_color),
                ("active", self.app._bg_color),
                ("!selected", self.app._text_fg),
            ],
            bordercolor=[
                ("selected", self.app._accent_color),
                ("pressed", self.app._accent_color),
                ("!selected", self.app._surface_color),
            ],
            lightcolor=[
                ("selected", self.app._accent_color),
                ("pressed", self.app._accent_color),
                ("!selected", self.app._surface_color),
            ],
            darkcolor=[
                ("selected", self.app._accent_color),
                ("pressed", self.app._accent_color),
                ("!selected", self.app._surface_color),
            ],
            focuscolor=[
                ("selected", self.app._accent_color),
                ("pressed", self.app._accent_color),
                ("!selected", self.app._surface_color),
            ],
            relief=[
                ("selected", "flat"),
                ("pressed", "flat"),
            ],
        )

    def build_ui(self) -> None:
        self.builder.build_ui()
