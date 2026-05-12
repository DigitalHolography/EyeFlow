"""Run-tab actions and background pipeline execution."""

import os
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from queue import Empty, Queue
from threading import Thread
from tkinter import filedialog, messagebox

from input_output import resolve_holo_input
from input_output.holo_run_layout import HoloRunLayout
from input_output.output_manager import OutputManager, OutputType
from pipelines import PipelineDescriptor
from pipeline_engine import PipelineExecutionPlan, run_pipelines_to_output


@dataclass(frozen=True)
class _PipelineRunRequest:
    output_manager: OutputManager
    pipelines: Sequence[PipelineDescriptor]
    target_names: Sequence[str]
    holodoppler_h5: Path | None
    doppler_vision_h5: Path | None


_PipelineUiEvent = tuple[str, object]


class RunMixin:
    def choose_holo_file(self) -> None:
        selected_holo = self._selected_holo_path()
        initial_dir = (
            str(selected_holo.parent)
            if selected_holo is not None
            else os.path.abspath("example_file")
        )
        path = filedialog.askopenfilename(
            filetypes=[("HOLO", "*.holo"), ("All files", "*.*")],
            initialdir=initial_dir,
            title="Select .holo file",
        )
        if path:
            self._assign_holo_input_path(Path(path))

    def _validate_selected_input(
        self,
        holo_path: Path | None,
    ) -> HoloRunLayout | None:
        if holo_path is None:
            messagebox.showwarning(
                "Missing input",
                "Select one .holo file.",
            )
            return None

        try:
            run_layout = resolve_holo_input(holo_path)
        except ValueError as exc:
            messagebox.showerror("Invalid input", str(exc))
            return None
        except FileNotFoundError as exc:
            messagebox.showerror("Missing data", str(exc))
            return None

        return run_layout

    def run_process(self) -> None:
        if getattr(self, "_pipeline_run_active", False):
            messagebox.showwarning(
                "Run in progress",
                "Wait for the current run to finish.",
            )
            return

        self._reset_progress()
        request = self._build_pipeline_run_request()
        if request is None:
            return

        self._start_progress(
            len(request.pipelines),
            style_name=self._progress_primary_style,
            status_text="Running pipelines...",
        )
        self._start_pipeline_thread(request)

    def _build_pipeline_run_request(self) -> _PipelineRunRequest | None:
        plan = self._resolve_selected_run_plan()
        if plan is None:
            return None

        pipelines = list(plan.descriptors)
        if not self._reject_unavailable_pipelines(pipelines):
            return None

        run_layout = self._validate_selected_input(self._selected_holo_path())
        if run_layout is None:
            return None

        self._reset_run_log("Starting pipeline run...\n")
        self._log_run(f"[DAG] Targets -> {', '.join(plan.targets)}")
        self._log_run(f"[DAG] Execution order -> {', '.join(plan.names)}")
        self._log_run_layout(run_layout)

        return self._create_pipeline_run_request(
            run_layout=run_layout,
            pipelines=pipelines,
            target_names=plan.targets,
        )

    def _resolve_selected_run_plan(self) -> PipelineExecutionPlan | None:
        target_names = self._selected_target_pipeline_names()
        if not target_names:
            messagebox.showwarning(
                "No pipelines",
                "Select at least one pipeline in Pipeline Library.",
            )
            return None

        try:
            return self._resolve_pipeline_plan(target_names)
        except (RuntimeError, ValueError) as exc:
            messagebox.showerror(
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
        messagebox.showerror(
            "Pipeline unavailable",
            "The DAG requires unavailable pipeline(s):\n" + "\n".join(details),
        )
        return False

    def _log_run_layout(self, run_layout: HoloRunLayout) -> None:
        self._log_run(f"[INPUT] HOLO -> {run_layout.holo_path}")
        self._log_run(f"[INPUT] DATA DIR -> {run_layout.root_dir}")
        self._log_run(f"[RESOLVED] HD -> {run_layout.hd_h5}")
        self._log_run(f"[RESOLVED] DV -> {run_layout.dv_h5}")

    def _create_pipeline_run_request(
        self,
        *,
        run_layout: HoloRunLayout,
        pipelines: Sequence[PipelineDescriptor],
        target_names: Sequence[str],
    ) -> _PipelineRunRequest | None:
        try:
            output_manager = self._prepare_output_manager_for_input(
                run_layout.holo_path
            )
        except (OSError, RuntimeError) as exc:
            self._log_run(f"[FAIL] {exc}")
            self._set_minimal_status("Run failed.")
            messagebox.showerror("Output folder unavailable", str(exc))
            return None

        output_h5_path = output_manager.path_for(OutputType.H5)
        self._log_run(f"[OUTPUT] {output_h5_path}")
        return _PipelineRunRequest(
            output_manager=output_manager,
            pipelines=pipelines,
            target_names=target_names,
            holodoppler_h5=run_layout.hd_h5,
            doppler_vision_h5=run_layout.dv_h5,
        )

    def _start_pipeline_thread(self, request: _PipelineRunRequest) -> None:
        self._pipeline_ui_events = Queue()
        self._set_pipeline_run_active(True)
        thread = Thread(
            target=self._run_pipeline_worker,
            args=(request,),
            name="EyeFlowPipelineWorker",
            daemon=True,
        )
        self._pipeline_run_thread = thread
        thread.start()
        self.after(50, self._drain_pipeline_ui_events)

    def _run_pipeline_worker(self, request: _PipelineRunRequest) -> None:
        try:
            output_h5_path = self._run_pipelines_to_output(request)
        except Exception as exc:  # noqa: BLE001
            self._queue_pipeline_ui_event("failure", str(exc))
            return

        self._queue_pipeline_ui_event("success", output_h5_path)

    def _set_pipeline_run_active(self, active: bool) -> None:
        self._pipeline_run_active = active
        self._set_run_controls_enabled(not active)
        if not active:
            self._pipeline_run_thread = None

    def _set_run_controls_enabled(self, enabled: bool) -> None:
        state = "normal" if enabled else "disabled"
        for attr in ("minimal_run_button", "advanced_run_button"):
            button = getattr(self, attr, None)
            if button is not None:
                button.configure(state=state)

    def _queue_pipeline_ui_event(self, event_type: str, payload: object) -> None:
        events: Queue[_PipelineUiEvent] | None = getattr(
            self,
            "_pipeline_ui_events",
            None,
        )
        if events is not None:
            events.put((event_type, payload))

    def _drain_pipeline_ui_events(self) -> None:
        events: Queue[_PipelineUiEvent] | None = getattr(
            self,
            "_pipeline_ui_events",
            None,
        )
        if events is None:
            return

        self._handle_pending_pipeline_events(events)
        if getattr(self, "_pipeline_run_active", False):
            self.after(50, self._drain_pipeline_ui_events)
        else:
            self._pipeline_ui_events = None

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
            self._log_run(str(payload))
        elif event_type == "progress":
            self._advance_progress()
        elif event_type == "success":
            self._finish_pipeline_run_success(Path(payload))
        elif event_type == "failure":
            self._finish_pipeline_run_failure(str(payload))

    def _finish_pipeline_run_success(self, output_h5_path: Path) -> None:
        self._set_progress_units(self._progress_total_units)
        self._log_run(f"Completed. Output file: {output_h5_path}")
        self._set_minimal_status("Process ended.")
        self._set_pipeline_run_active(False)

    def _finish_pipeline_run_failure(self, failure_message: str) -> None:
        self._log_run(f"[FAIL] {failure_message}")
        self._set_minimal_status("Run failed.")
        self._set_pipeline_run_active(False)
        messagebox.showerror("Run failed", failure_message)

    def _run_pipelines_to_output(
        self,
        request: _PipelineRunRequest,
    ) -> Path:
        return run_pipelines_to_output(
            output_manager=request.output_manager,
            pipelines=request.pipelines,
            target_names=request.target_names,
            holodoppler_h5=request.holodoppler_h5,
            doppler_vision_h5=request.doppler_vision_h5,
            on_pipeline_success=lambda name: self._queue_pipeline_ui_event(
                "log",
                f"[OK] {name}",
            ),
            on_progress=lambda: self._queue_pipeline_ui_event("progress", None),
        )
