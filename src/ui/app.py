"""Tk application shell for EyeFlow."""

from __future__ import annotations

import tkinter as tk
import tkinter.font as tkfont
from pathlib import Path

from runtime_limits import configure_numeric_threads

configure_numeric_threads()

from app_settings import AppSettingsStore, app_display_name
from pipelines import PipelineDescriptor
from pipeline_engine import PipelineDAG

from .controllers import (
    InputController,
    PipelineLibraryController,
    ProgressController,
    ResourceController,
    RunController,
    SettingsController,
    ViewController,
)
from .services import UiServices
from utils.logger import configure_logger

try:
    from tkinterdnd2 import TkinterDnD
except ImportError:  # optional dependency
    TkinterDnD = None

BaseAppTk = TkinterDnD.Tk if TkinterDnD is not None else tk.Tk


class ProcessApp(BaseAppTk):
    def __init__(self) -> None:
        super().__init__()
        self.title(app_display_name())
        self.ui_services = UiServices()
        self.settings_store = AppSettingsStore()
        configure_logger(self.settings_store.path)
        self._settings_warning_shown = False
        self.settings_controller = SettingsController(self)
        self.settings_controller.ensure_default_settings()
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
        self.holo_hd_status_var = tk.StringVar(value="")
        self.holo_dv_status_var = tk.StringVar(value="")
        self._progress_total_units = 1.0
        self._progress_completed_units = 0.0
        self._pipeline_run_active = False
        self._pipeline_run_thread = None
        self._pipeline_ui_events = None
        self._last_saved_run_log_path: Path | None = None
        self._progress_primary_style = "MinimalPrimary.Horizontal.TProgressbar"
        self._progress_final_style = "MinimalFinal.Horizontal.TProgressbar"
        self._window_icon_image: tk.PhotoImage | None = None
        self._minimal_logo_image: tk.PhotoImage | None = None
        self._minimal_title_font: tkfont.Font | None = None
        self._trim_h5source = tk.BooleanVar(
            value=self.settings_store.load_trim_h5source()
        )
        self.resource_controller = ResourceController(self)
        self.view_controller = ViewController(self)
        self.progress_controller = ProgressController(self)
        self.input_controller = InputController(self)
        self.pipeline_library_controller = PipelineLibraryController(self)
        self.run_controller = RunController(self)

        self.view_controller.set_initial_window_size()
        self.view_controller.apply_theme()
        self.resource_controller.set_window_icon()
        self.view_controller.build_ui()
        self.input_controller.install_drop_targets()
        self.holo_input_var.trace_add(
            "write",
            self.input_controller.on_holo_input_changed,
        )
        self.protocol("WM_DELETE_WINDOW", self.settings_controller.on_close)
        self.pipeline_library_controller.register()
        self.progress_controller.reset_run_log()
        self.input_controller.update_minimal_path_labels()
        self.settings_controller.apply_ui_mode(self.ui_mode, persist=False)


def main() -> None:
    app = ProcessApp()
    app.mainloop()
