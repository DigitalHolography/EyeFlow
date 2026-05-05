from __future__ import annotations

import tkinter as tk
import tkinter.font as tkfont
from pathlib import Path
from tkinter import ttk

from runtime_limits import configure_numeric_threads

configure_numeric_threads()

from app_settings import AppSettingsStore
from pipelines import PipelineDescriptor
from pipeline_engine import PipelineDAG

from .compat import BaseAppTk, sv_ttk
from .drag_drop import DragDropMixin
from .input_state import InputStateMixin
from .pipeline_library import PipelineLibraryMixin
from .progress import ProgressMixin
from .resources import ResourceMixin
from .run import RunMixin
from .settings import SettingsMixin
from .views import ViewBuilderMixin


class ProcessApp(
    ViewBuilderMixin,
    DragDropMixin,
    ResourceMixin,
    SettingsMixin,
    InputStateMixin,
    ProgressMixin,
    PipelineLibraryMixin,
    RunMixin,
    BaseAppTk,
):
    def __init__(self) -> None:
        super().__init__()
        self.title("EyeFlow")
        self.settings_store = AppSettingsStore()
        self._settings_warning_shown = False
        self._ensure_default_settings()
        self.ui_mode = self.settings_store.load_ui_mode()
        self.pipeline_registry: dict[str, PipelineDescriptor] = {}
        self.pipeline_catalog: dict[str, PipelineDescriptor] = {}
        self.pipeline_dag: PipelineDAG | None = None
        self.pipeline_rows: list[PipelineDescriptor] = []
        self.pipeline_visibility: dict[str, bool] = {}
        self.pipeline_visibility_vars: dict[str, tk.BooleanVar] = {}
        self.pipeline_row_widgets: dict[str, tk.Widget] = {}
        self.holo_input_var = tk.StringVar()
        self.run_progress_var = tk.DoubleVar(value=0.0)
        self._selected_holo_input_paths: list[Path] = []
        self._synchronizing_holo_input_var = False
        self.minimal_status_var = tk.StringVar(value="Ready.")
        self.pipeline_library_summary_var = tk.StringVar(value="")
        self.minimal_holo_input_path_var = tk.StringVar(value="No input selected")
        self.holo_hd_status_var = tk.StringVar(value="HD waiting")
        self.holo_dv_status_var = tk.StringVar(value="DV waiting")
        self._progress_total_units = 1.0
        self._progress_completed_units = 0.0
        self._last_saved_run_log_path: Path | None = None
        self._progress_primary_style = "MinimalPrimary.Horizontal.TProgressbar"
        self._progress_final_style = "MinimalFinal.Horizontal.TProgressbar"
        self._window_icon_image: tk.PhotoImage | None = None
        self._minimal_logo_image: tk.PhotoImage | None = None
        self._minimal_title_font: tkfont.Font | None = None
        self._trim_h5source = tk.BooleanVar(
            value=self.settings_store.load_trim_h5source()
        )

        self._set_initial_window_size()
        self._apply_theme()
        self._set_window_icon()
        self._build_ui()
        self._install_drop_targets()
        self.holo_input_var.trace_add("write", self._on_holo_input_changed)
        self.protocol("WM_DELETE_WINDOW", self._on_close)
        self._register_pipelines()
        self._reset_run_log()
        self._update_minimal_path_labels()
        self._apply_ui_mode(self.ui_mode, persist=False)

    def _set_initial_window_size(self) -> None:
        width, height, min_width, min_height = self._window_size_for_mode(self.ui_mode)
        screen_width = self.winfo_screenwidth()
        screen_height = self.winfo_screenheight()
        width = min(width, screen_width)
        height = min(height, screen_height)
        x = max((screen_width - width) // 2, 0)
        y = max((screen_height - height) // 2, 0)
        self.geometry(f"{width}x{height}+{x}+{y}")
        self.minsize(min_width, min_height)

    def _apply_theme(self) -> None:
        """
        Apply the Sun Valley ttk theme when available; otherwise fall back to a simple dark palette.
        """
        style = ttk.Style(self)
        self._style = style
        if sv_ttk:
            try:
                sv_ttk.set_theme("dark")
            except Exception:
                pass

        # Fallback palette aligned with Sun Valley dark.
        fallback_bg = "#0f1116"
        fallback_surface = "#1b1f27"
        fallback_fg = "#e8eef5"
        fallback_muted = "#9aa6b5"
        fallback_accent = "#4f9dff"

        # Derive colors from the active theme when possible to keep consistency.
        bg = style.lookup("TFrame", "background") or fallback_bg
        fg = style.lookup("TLabel", "foreground") or fallback_fg
        surface = (
            style.lookup("TEntry", "fieldbackground")
            or style.lookup("TEntry", "background")
            or fallback_surface
        )
        muted = (
            style.lookup("TLabel", "foreground", state=("disabled",)) or fallback_muted
        )
        accent = (
            style.lookup("TButton", "bordercolor")
            or style.lookup("TNotebook", "foreground")
            or fallback_accent
        )

        self.configure(bg=bg)
        # set texts colors when created.
        self._text_bg = surface
        self._text_fg = fg
        self._muted_fg = muted
        self._bg_color = bg
        self._surface_color = surface
        self._accent_color = accent
        self._success_color = "#3fb37f"
        self._error_color = "#ff6b6b"
        self._configure_progress_styles()

    def _configure_progress_styles(self) -> None:
        progress_colors = {
            self._progress_primary_style: self._accent_color,
            self._progress_final_style: "#3fb37f",
        }
        for style_name, color in progress_colors.items():
            try:
                self._style.configure(
                    style_name,
                    troughcolor=self._surface_color,
                    background=color,
                    bordercolor=color,
                    lightcolor=color,
                    darkcolor=color,
                )
            except tk.TclError:
                self._style.configure(style_name, background=color)


def main() -> None:
    app = ProcessApp()
    app.mainloop()
