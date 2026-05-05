from __future__ import annotations

import tkinter as tk
from tkinter import ttk

from app_settings import normalize_pipeline_visibility
from pipelines import PipelineDescriptor, load_pipeline_catalog
from pipeline_engine import PipelineDAG, PipelineExecutionPlan

from .widgets import _Tooltip


class PipelineLibraryMixin:
    def _build_pipeline_library_tab(self, parent: ttk.Frame) -> None:
        parent.columnconfigure(0, weight=1)
        parent.rowconfigure(2, weight=1)

        ttk.Label(parent, text="Select pipeline targets.").grid(
            row=0, column=0, sticky="w"
        )

        controls = ttk.Frame(parent)
        controls.grid(row=1, column=0, sticky="ew", pady=(8, 4))
        controls.columnconfigure(2, weight=1)
        ttk.Button(
            controls,
            text="Select all",
            command=self.select_all_pipelines,
        ).grid(row=0, column=0, sticky="w")
        ttk.Button(
            controls,
            text="Deselect all",
            command=self.deselect_all_pipelines,
        ).grid(row=0, column=1, sticky="w", padx=(4, 0))
        ttk.Label(controls, textvariable=self.pipeline_library_summary_var).grid(
            row=0, column=2, sticky="e"
        )

        library_container = ttk.Frame(parent)
        library_container.grid(row=2, column=0, sticky="nsew")
        library_container.columnconfigure(0, weight=1)
        library_container.rowconfigure(0, weight=1)

        self.pipeline_library_canvas = tk.Canvas(
            library_container, highlightthickness=0, bg=self._bg_color
        )
        self.pipeline_library_canvas.grid(row=0, column=0, sticky="nsew")
        library_scroll = ttk.Scrollbar(
            library_container,
            orient="vertical",
            command=self.pipeline_library_canvas.yview,
        )
        library_scroll.grid(row=0, column=1, sticky="ns")
        self.pipeline_library_canvas.configure(yscrollcommand=library_scroll.set)
        self.pipeline_library_inner = ttk.Frame(self.pipeline_library_canvas)
        self.pipeline_library_window = self.pipeline_library_canvas.create_window(
            (0, 0), window=self.pipeline_library_inner, anchor="nw"
        )
        self.pipeline_library_inner.bind(
            "<Configure>",
            lambda _evt: self.pipeline_library_canvas.configure(
                scrollregion=self.pipeline_library_canvas.bbox("all")
            ),
        )
        self.pipeline_library_canvas.bind(
            "<Configure>",
            lambda evt: self.pipeline_library_canvas.itemconfigure(
                self.pipeline_library_window, width=evt.width
            ),
        )
        self._bind_vertical_mousewheel(
            self.pipeline_library_canvas, self.pipeline_library_canvas
        )
        self._bind_vertical_mousewheel(
            self.pipeline_library_inner, self.pipeline_library_canvas
        )
        self._bind_vertical_mousewheel(library_scroll, self.pipeline_library_canvas)

    def _bind_vertical_mousewheel(self, widget: tk.Misc, canvas: tk.Canvas) -> None:
        for sequence in ("<MouseWheel>", "<Button-4>", "<Button-5>"):
            widget.bind(
                sequence,
                lambda event, target_canvas=canvas: self._on_vertical_mousewheel(
                    event, target_canvas
                ),
                add="+",
            )

    @staticmethod
    def _mousewheel_scroll_units(event: tk.Event) -> int:
        delta = int(getattr(event, "delta", 0) or 0)
        if delta:
            steps = max(1, abs(delta) // 120) if abs(delta) >= 120 else 1
            return -steps if delta > 0 else steps

        button = getattr(event, "num", None)
        if button == 4:
            return -1
        if button == 5:
            return 1
        return 0

    def _on_vertical_mousewheel(self, event: tk.Event, canvas: tk.Canvas) -> str | None:
        scroll_units = self._mousewheel_scroll_units(event)
        if not scroll_units:
            return None
        canvas.yview_scroll(scroll_units, "units")
        return "break"

    def _register_pipelines(self) -> None:
        available, missing = load_pipeline_catalog()
        rows = sorted(
            [*available, *missing], key=lambda pipeline: pipeline.name.lower()
        )
        self.pipeline_registry = {p.name: p for p in available}
        self.pipeline_catalog = {p.name: p for p in rows}

        try:
            self.pipeline_dag = PipelineDAG(rows)
        except (RuntimeError, ValueError) as exc:
            self.pipeline_dag = None
            self._show_settings_warning(
                "Pipeline DAG error",
                f"Pipeline dependency graph is invalid:\n{exc}",
            )
        else:
            rows = list(self.pipeline_dag.ordered_descriptors)

        self.pipeline_rows = rows
        self._sync_pipeline_visibility(rows)
        self._populate_pipeline_library(rows)
        self._install_drop_targets()

    def _descriptor_tooltip_text(self, descriptor: PipelineDescriptor) -> str:
        parts: list[str] = []
        description = getattr(descriptor, "description", "")
        if description:
            parts.append(description)

        dag_requires = getattr(descriptor, "dag_requires", ())
        if dag_requires:
            parts.append(f"DAG requires: {', '.join(dag_requires)}")

        dag_produces = getattr(descriptor, "dag_produces", ())
        if dag_produces:
            parts.append(f"DAG produces: {', '.join(dag_produces)}")

        missing_deps = getattr(descriptor, "missing_deps", []) or getattr(
            descriptor, "requires", []
        )
        if missing_deps:
            parts.append(f"Install: {', '.join(missing_deps)}")
        return "\n".join(parts)

    def _pipeline_status_text(self, pipeline: PipelineDescriptor) -> str:
        if pipeline.available:
            return "Available"
        if pipeline.missing_deps:
            return f"Missing deps: {', '.join(pipeline.missing_deps)}"
        return "Unavailable"

    def _populate_pipeline_library(self, rows: list[PipelineDescriptor]) -> None:
        for child in self.pipeline_library_inner.winfo_children():
            child.destroy()
        self.pipeline_visibility_vars = {}
        self.pipeline_row_widgets = {}
        self.pipeline_library_inner.columnconfigure(0, weight=1)

        selected_header = ttk.Label(self.pipeline_library_inner, text="Target")
        selected_header.grid(row=0, column=0, sticky="w", pady=(0, 6))
        order_header = ttk.Label(self.pipeline_library_inner, text="Pipeline")
        order_header.grid(row=0, column=1, sticky="w", padx=(12, 0), pady=(0, 6))
        status_header = ttk.Label(self.pipeline_library_inner, text="Status")
        status_header.grid(row=0, column=2, sticky="w", padx=(12, 0), pady=(0, 6))
        self._bind_vertical_mousewheel(selected_header, self.pipeline_library_canvas)
        self._bind_vertical_mousewheel(order_header, self.pipeline_library_canvas)
        self._bind_vertical_mousewheel(status_header, self.pipeline_library_canvas)

        for idx, pipeline in enumerate(rows, start=1):
            is_available = getattr(pipeline, "available", True)
            var = tk.BooleanVar(
                value=self.pipeline_visibility.get(pipeline.name, False)
                and is_available
            )
            row_frame = ttk.Frame(self.pipeline_library_inner)
            row_frame.grid(
                row=idx,
                column=0,
                columnspan=3,
                sticky="ew",
                pady=(0, 6),
            )
            row_frame.columnconfigure(1, weight=1)

            check = ttk.Checkbutton(
                row_frame,
                text="",
                variable=var,
                state="normal" if is_available else "disabled",
                command=lambda name=pipeline.name, visible_var=var: (
                    self._set_pipeline_visibility(name, visible_var.get())
                ),
            )
            check.grid(row=0, column=0, sticky="w")
            name_label = ttk.Label(row_frame, text=pipeline.name)
            name_label.grid(row=0, column=1, sticky="w", padx=(12, 0))
            self._bind_vertical_mousewheel(check, self.pipeline_library_canvas)
            self._bind_vertical_mousewheel(name_label, self.pipeline_library_canvas)
            self._bind_vertical_mousewheel(row_frame, self.pipeline_library_canvas)
            self.pipeline_row_widgets[pipeline.name] = row_frame

            status_text = self._pipeline_status_text(pipeline)
            status = ttk.Label(row_frame, text=status_text)
            status.grid(row=0, column=2, sticky="w", padx=(12, 0))
            self._bind_vertical_mousewheel(status, self.pipeline_library_canvas)

            tip_text = self._descriptor_tooltip_text(pipeline)
            if tip_text:
                _Tooltip(check, tip_text, bg=self._surface_color, fg=self._text_fg)
                _Tooltip(name_label, tip_text, bg=self._surface_color, fg=self._text_fg)
                _Tooltip(status, tip_text, bg=self._surface_color, fg=self._text_fg)

            self.pipeline_visibility_vars[pipeline.name] = var

        self._update_pipeline_library_summary()

    def _sync_pipeline_visibility(self, rows: list[PipelineDescriptor]) -> None:
        visibility, changed = normalize_pipeline_visibility(
            (pipeline.name for pipeline in rows),
            self.settings_store.load_pipeline_visibility(),
        )
        for pipeline in rows:
            if not pipeline.available and visibility.get(pipeline.name, False):
                visibility[pipeline.name] = False
                changed = True
        self.pipeline_visibility = visibility
        if changed:
            self._persist_pipeline_visibility()

    def _persist_pipeline_visibility(self) -> None:
        try:
            self.settings_store.save_pipeline_visibility(self.pipeline_visibility)
        except OSError as exc:
            self._show_settings_warning(
                "Settings not saved",
                f"Could not save pipeline selection preferences:\n{exc}",
            )

    def _set_pipeline_visibility(self, name: str, visible: bool) -> None:
        pipeline = self.pipeline_catalog.get(name)
        if pipeline is not None and not pipeline.available:
            visible = False
        if self.pipeline_visibility.get(name) == visible:
            return
        self.pipeline_visibility[name] = visible
        self._persist_pipeline_visibility()
        self._update_pipeline_library_summary()

    def _set_all_pipeline_visibility(self, visible: bool) -> None:
        changed = False
        target_values = {
            pipeline.name: visible and pipeline.available
            for pipeline in self.pipeline_rows
        }
        for name, target_value in target_values.items():
            if self.pipeline_visibility.get(name) != target_value:
                self.pipeline_visibility[name] = target_value
                changed = True
        if not changed:
            return
        for name, var in self.pipeline_visibility_vars.items():
            var.set(self.pipeline_visibility.get(name, False))
        self._persist_pipeline_visibility()
        self._update_pipeline_library_summary()

    def _selected_target_pipeline_names(self) -> list[str]:
        return [
            pipeline.name
            for pipeline in self.pipeline_rows
            if pipeline.available and self.pipeline_visibility.get(pipeline.name, False)
        ]

    def _resolve_pipeline_plan(
        self,
        target_names: list[str],
    ) -> PipelineExecutionPlan:
        dag = getattr(self, "pipeline_dag", None)
        if dag is None:
            dag = PipelineDAG(self.pipeline_rows)
            self.pipeline_dag = dag
        return dag.resolve_targets(target_names)

    def _update_pipeline_library_summary(self) -> None:
        target_names = self._selected_target_pipeline_names()
        available_count = sum(
            1 for pipeline in self.pipeline_rows if pipeline.available
        )
        if not target_names:
            self.pipeline_library_summary_var.set(f"Targets: 0/{available_count}")
            return

        try:
            plan = self._resolve_pipeline_plan(target_names)
        except (RuntimeError, ValueError):
            self.pipeline_library_summary_var.set("DAG error")
            return

        self.pipeline_library_summary_var.set(
            f"Targets: {len(target_names)}/{available_count} | "
            f"Run steps: {len(plan.descriptors)}"
        )

    def select_all_pipelines(self) -> None:
        self._set_all_pipeline_visibility(True)

    def deselect_all_pipelines(self) -> None:
        self._set_all_pipeline_visibility(False)
