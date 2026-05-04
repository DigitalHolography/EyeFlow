from __future__ import annotations

import os
from collections.abc import Sequence
from contextlib import ExitStack
from pathlib import Path
from tkinter import filedialog, messagebox

from input_output import (
    EyeFlowOutputManager,
    PipelineInputView,
    ResolvedHoloInput,
    open_h5,
    resolve_selected_holo_inputs,
)
from pipelines import PipelineDescriptor
from pipelines.core.errors import format_pipeline_exception


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

        selected_names = [
            pipeline.name
            for pipeline in self.pipeline_rows
            if pipeline.available and self.pipeline_visibility.get(pipeline.name, False)
        ]
        if not selected_names:
            messagebox.showwarning(
                "No pipelines",
                "Select at least one pipeline in Pipeline Library.",
            )
            return

        pipelines: list[PipelineDescriptor] = []
        missing: list[str] = []
        for name in selected_names:
            pipeline = self.pipeline_registry.get(name)
            if pipeline is None:
                missing.append(name)
            else:
                pipelines.append(pipeline)
        if missing:
            messagebox.showerror(
                "Pipeline missing", f"Pipeline(s) not registered: {', '.join(missing)}"
            )
            return

        resolved_input = self._validate_selected_input(holo_path)
        if resolved_input is None:
            return

        self._reset_run_log("Starting pipeline run...\n")
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
        holodoppler_h5: Path | None,
        doppler_vision_h5: Path | None,
    ) -> Path:
        output_h5_path.parent.mkdir(parents=True, exist_ok=True)
        with ExitStack() as stack:
            work_h5 = stack.enter_context(open_h5(output_h5_path, "w"))
            hd_h5 = (
                stack.enter_context(open_h5(holodoppler_h5, "r"))
                if holodoppler_h5 is not None
                else None
            )
            dv_h5 = (
                stack.enter_context(open_h5(doppler_vision_h5, "r"))
                if doppler_vision_h5 is not None
                else None
            )
            output_manager = EyeFlowOutputManager(work_h5)
            output_manager.initialize(
                holodoppler_source_file=(
                    str(holodoppler_h5) if holodoppler_h5 is not None else None
                ),
                doppler_vision_source_file=(
                    str(doppler_vision_h5) if doppler_vision_h5 is not None else None
                ),
            )
            work_h5.attrs["trim_h5source"] = True
            work_h5.attrs["pipeline_order"] = [pipeline.name for pipeline in pipelines]

            for pipeline_desc in pipelines:
                pipeline = pipeline_desc.instantiate()
                pipeline_input = PipelineInputView(
                    work_h5=work_h5,
                    holodoppler_h5=hd_h5,
                    doppler_vision_h5=dv_h5,
                )
                try:
                    result = pipeline.run(pipeline_input)
                except Exception as exc:  # noqa: BLE001
                    raise RuntimeError(
                        format_pipeline_exception(exc, pipeline)
                    ) from exc
                output_manager.append_pipeline_result(pipeline.name, result)
                result.output_h5_path = str(output_h5_path)
                self._log_run(f"[OK] {pipeline.name}")
                self._advance_progress()
        return output_h5_path
