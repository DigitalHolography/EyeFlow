"""Controller for Run-tab user actions."""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from queue import Empty, Queue
from threading import Event, Thread

from pipeline_engine import RunResult, RunSpec, execute_run, resolve_run_spec
from utils.logger import Logger, LogLevel

from ..services import services_for

_PipelineUiEvent = tuple[str, object]


@dataclass(frozen=True)
class _FileStarted:
    path: Path
    index: int
    total: int
    pipeline_count: int


@dataclass(frozen=True)
class _PipelineStarted:
    name: str
    index: int
    total: int


class RunController:
    def __init__(self, app) -> None:
        self.app = app
        self._stop_after_current_file = Event()
        self._current_file: _FileStarted | None = None
        self._current_pipeline: _PipelineStarted | None = None
        Logger.configure(on_log=self._queue_log)

    def _queue_log(self, message: str) -> None:
        self._queue_pipeline_ui_event("log", message)

    def choose_holo_file(self) -> None:
        selected_holo = self.app.input_controller.selected_holo_path()
        initial_dir = (
            str(selected_holo.parent)
            if selected_holo is not None
            else os.path.abspath("example_file")
        )
        paths = services_for(self.app).file_dialogs.askopenfilenames(
            filetypes=[("HOLO or path list", "*.holo *.txt"), ("All files", "*.*")],
            initialdir=initial_dir,
            title="Select .holo file(s) or one .txt path list",
        )
        if paths:
            self.app.input_controller.assign_holo_input_paths(
                [Path(path) for path in paths]
            )

    def run_process(self) -> None:
        if getattr(self.app, "_pipeline_run_active", False):
            services_for(self.app).dialogs.showwarning(
                "Run in progress",
                "Wait for the current run to finish.",
            )
            return

        self.app.progress_controller.reset_progress()
        spec = self._build_run_spec()
        if spec is None:
            return

        self.app.progress_controller.start_progress(
            len(spec.plan.descriptors),
            style_name=self.app._progress_primary_style,
            status_text=f"Starting batch of {len(spec.requests)} file(s)...",
        )
        self.start_pipeline_thread(spec)

    def _build_run_spec(self) -> RunSpec | None:
        input_paths = self.app.input_controller.selected_holo_paths()
        if not input_paths:
            services_for(self.app).dialogs.showwarning(
                "Missing input",
                "Select one or more .holo files or one .txt path list.",
            )
            return None
        target_names = (
            self.app.pipeline_library_controller.selected_target_pipeline_names()
        )
        if not target_names:
            services_for(self.app).dialogs.showwarning(
                "No pipelines",
                "Select at least one pipeline in Pipeline Library.",
            )
            return None

        try:
            spec = resolve_run_spec(
                input_paths=input_paths,
                target_names=target_names,
                pipelines=self.app.pipeline_catalog.values(),
                pipeline_options=(
                    self.app.pipeline_library_controller.selected_pipeline_options()
                ),
            )
        except (FileNotFoundError, OSError, RuntimeError, ValueError) as exc:
            services_for(self.app).dialogs.showerror(
                "Run configuration error",
                str(exc),
            )
            return None
        self.app.progress_controller.reset_run_log("Starting pipeline run...\n")
        return spec

    def start_pipeline_thread(self, spec: RunSpec) -> None:
        self._stop_after_current_file.clear()
        self._current_file = None
        self._current_pipeline = None
        self.app._pipeline_ui_events = Queue()
        self._set_pipeline_run_active(True)
        thread = Thread(
            target=self._run_pipeline_worker,
            args=(spec,),
            name="EyeFlowPipelineWorker",
            daemon=True,
        )
        self.app._pipeline_run_thread = thread
        thread.start()
        self.app.after(50, self.drain_pipeline_ui_events)

    def stop_process(self) -> None:
        if not getattr(self.app, "_pipeline_run_active", False):
            return
        if self._stop_after_current_file.is_set():
            return

        self._stop_after_current_file.set()
        self._set_stop_controls_enabled(False)
        self.app.progress_controller.set_minimal_status(
            self._current_progress_status_text()
        )
        self.app.progress_controller.log_run(
            "[STOP] Stop requested; the current input file will finish."
        )

    def _run_pipeline_worker(self, spec: RunSpec) -> None:
        try:
            summary = self._execute_run(spec)
        except Exception as exc:  # noqa: BLE001
            self._queue_pipeline_ui_event("failure", str(exc))
            return

        self._queue_pipeline_ui_event("success", summary)

    def _set_pipeline_run_active(self, active: bool) -> None:
        self.app._pipeline_run_active = active
        self._set_run_controls_enabled(not active)
        self._set_stop_controls_enabled(active)
        if not active:
            self.app._pipeline_run_thread = None

    def _set_run_controls_enabled(self, enabled: bool) -> None:
        state = "normal" if enabled else "disabled"
        for attr in ("minimal_run_button", "advanced_run_button"):
            button = getattr(self.app, attr, None)
            if button is not None:
                button.configure(state=state)

    def _set_stop_controls_enabled(self, enabled: bool) -> None:
        state = "normal" if enabled else "disabled"
        for attr in ("minimal_stop_button", "advanced_stop_button"):
            button = getattr(self.app, attr, None)
            if button is not None:
                button.configure(state=state)

    def _queue_pipeline_ui_event(self, event_type: str, payload: object) -> None:
        events: Queue[_PipelineUiEvent] | None = getattr(
            self.app,
            "_pipeline_ui_events",
            None,
        )
        if events is not None:
            events.put((event_type, payload))

    def drain_pipeline_ui_events(self) -> None:
        events: Queue[_PipelineUiEvent] | None = getattr(
            self.app,
            "_pipeline_ui_events",
            None,
        )
        if events is None:
            return

        self._handle_pending_pipeline_events(events)
        if getattr(self.app, "_pipeline_run_active", False):
            self.app.after(50, self.drain_pipeline_ui_events)
        else:
            self.app._pipeline_ui_events = None

    def _handle_pending_pipeline_events(
        self,
        events: Queue[_PipelineUiEvent],
    ) -> None:
        while True:
            try:
                event_type, payload = events.get_nowait()
            except Empty:
                return
            self._handle_pipeline_ui_event(event_type, payload)

    def _handle_pipeline_ui_event(self, event_type: str, payload: object) -> None:
        if event_type == "log":
            self.app.progress_controller.log_run(str(payload))
        elif event_type == "file_start":
            assert isinstance(payload, _FileStarted)
            self._handle_file_started(payload)
        elif event_type == "pipeline_start":
            assert isinstance(payload, _PipelineStarted)
            self._handle_pipeline_started(payload)
        elif event_type == "progress":
            self.app.progress_controller.advance_progress()
        elif event_type == "success":
            self._finish_pipeline_run_success(payload)
        elif event_type == "failure":
            self._finish_pipeline_run_failure(str(payload))

    def _handle_file_started(self, event: _FileStarted) -> None:
        self._current_file = event
        self._current_pipeline = None
        self.app.progress_controller.start_progress(
            event.pipeline_count,
            style_name=self.app._progress_primary_style,
            status_text=self._current_progress_status_text(),
        )

    def _handle_pipeline_started(self, event: _PipelineStarted) -> None:
        self._current_pipeline = event
        self.app.progress_controller.set_minimal_status(
            self._current_progress_status_text()
        )

    def _current_progress_status_text(self) -> str:
        current_file = self._current_file
        current_pipeline = self._current_pipeline
        if current_file is None:
            if self._stop_after_current_file.is_set():
                return "Stop requested. Finishing the current file..."
            return "Processing input..."

        file_status = (
            f"File {current_file.index}/{current_file.total}: {current_file.path}"
        )
        detail_status = (
            f"Pipeline {current_pipeline.index}/{current_pipeline.total}: "
            f"{current_pipeline.name}"
            if current_pipeline is not None
            else "Preparing input..."
        )
        if self._stop_after_current_file.is_set():
            detail_status += " — Stop requested; finishing current file..."
        return f"{file_status}\n{detail_status}"

    def _finish_pipeline_run_success(self, summary: object) -> None:
        assert isinstance(summary, RunResult)
        if summary.last_output_path is not None:
            self.app.progress_controller.log_run(
                f"Completed. Output file: {summary.last_output_path}"
            )
        if summary.failures:
            self.app.progress_controller.log_run(
                f"{len(summary.failures)} failure(s):"
            )
            for failure in summary.failures:
                self.app.progress_controller.log_run(f" - {failure}")
        if summary.stopped:
            self.app.progress_controller.set_minimal_status(
                "Batch stopped after the current file."
            )
        elif summary.failures:
            self.app.progress_controller.set_minimal_status(
                f"Process ended with {len(summary.failures)} failure(s)."
            )
        else:
            self.app.progress_controller.set_minimal_status("Process ended.")
        self._set_pipeline_run_active(False)

    def _finish_pipeline_run_failure(self, failure_message: str) -> None:
        self.app.progress_controller.log_run(
            failure_message,
            level=LogLevel.ERROR,
        )
        self.app.progress_controller.set_minimal_status("Run failed.")
        self._set_pipeline_run_active(False)
        services_for(self.app).dialogs.showerror("Run failed", failure_message)

    def _execute_run(self, spec: RunSpec) -> RunResult:
        pipeline_count = len(spec.plan.descriptors)
        return execute_run(
            spec,
            on_file_start=lambda path, index, total: self._queue_pipeline_ui_event(
                "file_start",
                _FileStarted(path, index, total, pipeline_count),
            ),
            on_pipeline_start=lambda name, index, total: (
                self._queue_pipeline_ui_event(
                    "pipeline_start",
                    _PipelineStarted(name, index, total),
                )
            ),
            on_progress=lambda: self._queue_pipeline_ui_event("progress", None),
            should_stop=self._stop_after_current_file.is_set,
        )
