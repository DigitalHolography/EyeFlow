from __future__ import annotations

import os
from collections.abc import Sequence
from pathlib import Path
from tkinter import filedialog, messagebox

from input_output import (
    ResolvedHoloInput,
    resolve_selected_holo_inputs,
    run_pipelines_to_output_h5,
)
from pipelines import PipelineDescriptor


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
    ) -> ResolvedHoloInput | None:
        if holo_path is None:
            messagebox.showwarning(
                "Missing input",
                "Select one .holo file.",
            )
            return None

        try:
            resolved_inputs = resolve_selected_holo_inputs([holo_path])
        except ValueError as exc:
            messagebox.showerror("Invalid input", str(exc))
            return None
        except FileNotFoundError as exc:
            messagebox.showerror("Missing data", str(exc))
            return None

        return resolved_inputs[0]

    def run_process(self) -> None:
        self._reset_progress()
        holo_path = self._selected_holo_path()

        target_names = self._selected_target_pipeline_names()
        if not target_names:
            messagebox.showwarning(
                "No pipelines",
                "Select at least one pipeline in Pipeline Library.",
            )
            return

        try:
            plan = self._resolve_pipeline_plan(target_names)
        except (RuntimeError, ValueError) as exc:
            messagebox.showerror(
                "Pipeline DAG error",
                str(exc),
            )
            return

        pipelines = list(plan.descriptors)
        unavailable = [pipeline for pipeline in pipelines if not pipeline.available]
        if unavailable:
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
            return

        resolved_input = self._validate_selected_input(holo_path)
        if resolved_input is None:
            return

        self._reset_run_log("Starting pipeline run...\n")
        self._log_run(f"[DAG] Targets -> {', '.join(plan.targets)}")
        self._log_run(f"[DAG] Execution order -> {', '.join(plan.names)}")
        self._log_run(f"[INPUT] HOLO -> {resolved_input.holo_path}")
        self._log_run(f"[INPUT] DATA DIR -> {resolved_input.data_dir}")
        self._log_run(f"[RESOLVED] HD -> {resolved_input.hd_h5}")
        self._log_run(f"[RESOLVED] DV -> {resolved_input.dv_h5}")

        try:
            output_dir = self._prepare_default_output_dir(resolved_input.holo_path)
        except (OSError, RuntimeError) as exc:
            self._log_run(f"[FAIL] {exc}")
            self._set_minimal_status("Run failed.")
            messagebox.showerror("Output folder unavailable", str(exc))
            return

        output_h5_path = output_dir / self._default_work_h5_name_for_input(
            resolved_input.holo_path
        )
        self._log_run(f"[OUTPUT] {output_h5_path}")

        self._start_progress(
            len(pipelines),
            style_name=self._progress_primary_style,
            status_text="Running pipelines...",
        )

        try:
            self._run_pipelines_to_output(
                output_h5_path=output_h5_path,
                pipelines=pipelines,
                target_names=plan.targets,
                holodoppler_h5=resolved_input.hd_h5,
                doppler_vision_h5=resolved_input.dv_h5,
            )
        except Exception as exc:  # noqa: BLE001
            failure_message = str(exc)
            self._log_run(f"[FAIL] {failure_message}")
            self._set_minimal_status("Run failed.")
            messagebox.showerror("Run failed", failure_message)
            return

        self._set_progress_units(self._progress_total_units)
        self._log_run(f"Completed. Output file: {output_h5_path}")
        self._set_minimal_status("Process ended.")

    def _run_pipelines_to_output(
        self,
        *,
        output_h5_path: Path,
        pipelines: Sequence[PipelineDescriptor],
        target_names: Sequence[str] = (),
        holodoppler_h5: Path | None,
        doppler_vision_h5: Path | None,
    ) -> Path:
        return run_pipelines_to_output_h5(
            output_h5_path=output_h5_path,
            pipelines=pipelines,
            target_names=target_names,
            holodoppler_h5=holodoppler_h5,
            doppler_vision_h5=doppler_vision_h5,
            on_pipeline_success=lambda name: self._log_run(f"[OK] {name}"),
            on_progress=self._advance_progress,
        )
