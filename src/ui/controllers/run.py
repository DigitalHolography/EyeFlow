"""Controller for Run-tab user actions."""

from __future__ import annotations

import os
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from queue import Empty, Queue
from threading import Thread

from input_output import resolve_selected_run_layouts
from input_output.holo_run_layout import HoloRunLayout
from input_output.output_manager import OutputManager, OutputType
from pipelines import PipelineDescriptor
from pipeline_engine import PipelineExecutionPlan, run_pipelines_to_output

from ..services import services_for


@dataclass(frozen=True)
class _PipelineRunRequest:
    output_manager: OutputManager
    pipelines: Sequence[PipelineDescriptor]
    target_names: Sequence[str]
    holodoppler_h5: Path | None
    doppler_vision_h5: Path | None


@dataclass(frozen=True)
class _PipelineRunSummary:
    output_h5_path: Path | None
    failures: tuple[str, ...]


_PipelineUiEvent = tuple[str, object]


class RunController:
    def __init__(self, app) -> None:
        self.app = app

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
        requests = self._build_pipeline_run_requests()
        if requests is None:
            return

        total_units = sum(len(request.pipelines) for request in requests)
        self.app.progress_controller.start_progress(
            total_units,
            style_name=self.app._progress_primary_style,
            status_text="Running pipelines...",
        )
        self.start_pipeline_thread(requests)

    def _build_pipeline_run_requests(self) -> list[_PipelineRunRequest] | None:
        plan = self._resolve_selected_run_plan()
        if plan is None:
            return None

        pipelines = list(plan.descriptors)
        if not self._reject_unavailable_pipelines(pipelines):
            return None

        run_layouts = self._validate_selected_inputs(
            self.app.input_controller.selected_holo_paths()
        )
        if run_layouts is None:
            return None

        self.app.progress_controller.reset_run_log("Starting pipeline run...\n")
        self.app.progress_controller.log_run(
            f"[DAG] Targets -> {', '.join(plan.targets)}"
        )
        self.app.progress_controller.log_run(
            f"[DAG] Execution order -> {', '.join(plan.names)}"
        )
        for run_layout in run_layouts:
            self._log_run_layout(run_layout)

        return [
            self._create_pipeline_run_request(
                run_layout=run_layout,
                pipelines=pipelines,
                target_names=plan.targets,
            )
            for run_layout in run_layouts
        ]

    def _validate_selected_inputs(
        self,
        holo_paths: Sequence[Path],
    ) -> list[HoloRunLayout] | None:
        if not holo_paths:
            services_for(self.app).dialogs.showwarning(
                "Missing input",
                "Select one or more .holo files or one .txt path list.",
            )
            return None

        try:
            return resolve_selected_run_layouts(holo_paths)
        except ValueError as exc:
            services_for(self.app).dialogs.showerror("Invalid input", str(exc))
            return None
        except FileNotFoundError as exc:
            services_for(self.app).dialogs.showerror("Missing data", str(exc))
            return None

    def _resolve_selected_run_plan(self) -> PipelineExecutionPlan | None:
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
            return self.app.pipeline_library_controller.resolve_plan(target_names)
        except (RuntimeError, ValueError) as exc:
            services_for(self.app).dialogs.showerror(
                "Pipeline DAG error",
                str(exc),
            )
            return None

    def _reject_unavailable_pipelines(
        self,
        pipelines: Sequence[PipelineDescriptor],
    ) -> bool:
        unavailable = [pipeline for pipeline in pipelines if not pipeline.available]
        if not unavailable:
            return True

        details = []
        for pipeline in unavailable:
            reason = ", ".join(pipeline.missing_deps or pipeline.requires)
            details.append(
                f"{pipeline.name}" + (f" ({reason})" if reason else "")
            )
        services_for(self.app).dialogs.showerror(
            "Pipeline unavailable",
            "The DAG requires unavailable pipeline(s):\n" + "\n".join(details),
        )
        return False

    def _log_run_layout(self, run_layout: HoloRunLayout) -> None:
        self.app.progress_controller.log_run(f"[INPUT] HOLO -> {run_layout.holo_path}")
        self.app.progress_controller.log_run(
            f"[INPUT] DATA DIR -> {run_layout.root_dir}"
        )
        self.app.progress_controller.log_run(f"[RESOLVED] HD -> {run_layout.hd_h5}")
        self.app.progress_controller.log_run(f"[RESOLVED] DV -> {run_layout.dv_h5}")

    def _create_pipeline_run_request(
        self,
        *,
        run_layout: HoloRunLayout,
        pipelines: Sequence[PipelineDescriptor],
        target_names: Sequence[str],
    ) -> _PipelineRunRequest | None:
        try:
            output_manager = self.app.input_controller.prepare_output_manager_for_input(
                run_layout.holo_path
            )
        except (OSError, RuntimeError) as exc:
            self.app.progress_controller.log_run(f"[FAIL] {exc}")
            self.app.progress_controller.set_minimal_status("Run failed.")
            services_for(self.app).dialogs.showerror(
                "Output folder unavailable",
                str(exc),
            )
            return None

        output_h5_path = output_manager.path_for(OutputType.H5)
        self.app.progress_controller.log_run(f"[OUTPUT] {output_h5_path}")
        return _PipelineRunRequest(
            output_manager=output_manager,
            pipelines=pipelines,
            target_names=target_names,
            holodoppler_h5=run_layout.hd_h5,
            doppler_vision_h5=run_layout.dv_h5,
        )

    def start_pipeline_thread(self, requests: list[_PipelineRunRequest]) -> None:
        self.app._pipeline_ui_events = Queue()
        self._set_pipeline_run_active(True)
        thread = Thread(
            target=self._run_pipeline_worker,
            args=(requests,),
            name="EyeFlowPipelineWorker",
            daemon=True,
        )
        self.app._pipeline_run_thread = thread
        thread.start()
        self.app.after(50, self.drain_pipeline_ui_events)

    def _run_pipeline_worker(self, requests: list[_PipelineRunRequest]) -> None:
        try:
            summary = self._run_pipelines_to_output(requests)
        except Exception as exc:  # noqa: BLE001
            self._queue_pipeline_ui_event("failure", str(exc))
            return

        self._queue_pipeline_ui_event("success", summary)

    def _set_pipeline_run_active(self, active: bool) -> None:
        self.app._pipeline_run_active = active
        self._set_run_controls_enabled(not active)
        if not active:
            self.app._pipeline_run_thread = None

    def _set_run_controls_enabled(self, enabled: bool) -> None:
        state = "normal" if enabled else "disabled"
        for attr in ("minimal_run_button", "advanced_run_button"):
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
        elif event_type == "progress":
            self.app.progress_controller.advance_progress()
        elif event_type == "success":
            self._finish_pipeline_run_success(payload)
        elif event_type == "failure":
            self._finish_pipeline_run_failure(str(payload))

    def _finish_pipeline_run_success(self, summary: object) -> None:
        assert isinstance(summary, _PipelineRunSummary)
        self.app.progress_controller.set_progress_units(self.app._progress_total_units)
        if summary.output_h5_path is not None:
            self.app.progress_controller.log_run(
                f"Completed. Output file: {summary.output_h5_path}"
            )
        if summary.failures:
            self.app.progress_controller.log_run(
                f"{len(summary.failures)} failure(s):"
            )
            for failure in summary.failures:
                self.app.progress_controller.log_run(f" - {failure}")
            self.app.progress_controller.set_minimal_status(
                f"Process ended with {len(summary.failures)} failure(s)."
            )
        else:
            self.app.progress_controller.set_minimal_status("Process ended.")
        self._set_pipeline_run_active(False)

    def _finish_pipeline_run_failure(self, failure_message: str) -> None:
        self.app.progress_controller.log_run(f"[FAIL] {failure_message}")
        self.app.progress_controller.set_minimal_status("Run failed.")
        self._set_pipeline_run_active(False)
        services_for(self.app).dialogs.showerror("Run failed", failure_message)

    def _run_pipelines_to_output(
        self,
        requests: list[_PipelineRunRequest],
    ) -> _PipelineRunSummary:
        last_output_path: Path | None = None
        failures: list[str] = []
        for request in requests:
            run_name = request.output_manager.layout.holo_path.name
            try:
                last_output_path = run_pipelines_to_output(
                    output_manager=request.output_manager,
                    pipelines=request.pipelines,
                    target_names=request.target_names,
                    holodoppler_h5=request.holodoppler_h5,
                    doppler_vision_h5=request.doppler_vision_h5,
                    on_log=lambda message: self._queue_pipeline_ui_event(
                        "log",
                        message,
                    ),
                    on_progress=lambda: self._queue_pipeline_ui_event(
                        "progress",
                        None,
                    ),
                )
            except Exception as exc:  # noqa: BLE001
                failure = f"{run_name}: {exc}"
                failures.append(failure)
                self._queue_pipeline_ui_event("log", f"[FAIL] {failure}")
                continue

            self._queue_pipeline_ui_event(
                "log",
                f"Completed run for {run_name}: {last_output_path}",
            )

        return _PipelineRunSummary(last_output_path, tuple(failures))
