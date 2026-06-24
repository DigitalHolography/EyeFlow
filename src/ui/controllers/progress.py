"""Controller for progress state and run-log UI updates."""

from __future__ import annotations

from pathlib import Path

from utils.logger import current_logger

from ..services import services_for


class ProgressController:
    def __init__(self, app) -> None:
        self.app = app

    def set_minimal_status(self, text: str) -> None:
        self.app.minimal_status_var.set(text)

    def run_log_path(self) -> Path:
        return current_logger().path

    def persist_run_log_snapshot(self) -> None:
        try:
            logger = current_logger()
            logger.write_snapshot(self.app.run_log.get("1.0", "end-1c"))
        except OSError:
            self.app._last_saved_run_log_path = None
            return
        self.app._last_saved_run_log_path = logger.last_saved_path

    def open_run_log_file(self) -> None:
        log_path = self.run_log_path()
        try:
            self.persist_run_log_snapshot()
            current_logger().ensure_file()
            services_for(self.app).path_opener.open_path(log_path)
        except OSError as exc:
            services_for(self.app).dialogs.showerror(
                "Log file unavailable",
                f"Could not open the log file:\n{log_path}\n\n{exc}",
            )

    def set_progress_style(self, style_name: str) -> None:
        if hasattr(self.app, "minimal_progress"):
            self.app.minimal_progress.configure(style=style_name)

    def reset_progress(self) -> None:
        self.app._progress_total_units = 1.0
        self.app._progress_completed_units = 0.0
        self.set_progress_style(self.app._progress_primary_style)
        self.app.run_progress_var.set(0.0)

    def start_progress(
        self,
        total_units: float,
        *,
        style_name: str | None = None,
        status_text: str | None = None,
    ) -> None:
        self.app._progress_total_units = max(float(total_units), 1.0)
        self.app._progress_completed_units = 0.0
        self.set_progress_style(style_name or self.app._progress_primary_style)
        self.app.run_progress_var.set(0.0)
        if status_text is not None:
            self.set_minimal_status(status_text)

    def set_progress_units(self, completed_units: float) -> None:
        clamped_units = min(
            max(float(completed_units), 0.0),
            max(self.app._progress_total_units, 1.0),
        )
        self.app._progress_completed_units = clamped_units
        self.app.run_progress_var.set(
            (clamped_units / max(self.app._progress_total_units, 1.0)) * 100.0
        )

    def advance_progress(self, units: float = 1.0) -> None:
        self.set_progress_units(self.app._progress_completed_units + units)

    def reset_run_log(
        self,
        message: str = "Select one .holo file, then run.",
    ) -> None:
        self.app.run_log.configure(state="normal")
        self.app.run_log.delete("1.0", "end")
        self.app.run_log.insert("end", message)
        self.app.run_log.configure(state="disabled")
        self.persist_run_log_snapshot()

    def log_run(self, text: str) -> None:
        self.app.run_log.configure(state="normal")
        self.app.run_log.insert("end", f"{text}\n")
        self.app.run_log.see("end")
        self.app.run_log.configure(state="disabled")
        self.persist_run_log_snapshot()

    def show_run_error_dialog(self, message: str) -> None:
        self.app.bell()
        if self.app._last_saved_run_log_path is not None:
            message = (
                f"{message}\n\nLatest log saved to:\n"
                f"{self.app._last_saved_run_log_path}"
            )
        services_for(self.app).dialogs.showwarning(
            "Run completed with errors",
            message,
        )
