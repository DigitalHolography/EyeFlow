"""Controller for pipeline target selection and DAG summaries."""

from __future__ import annotations

import tkinter as tk
from tkinter import ttk

from app_settings import normalize_pipeline_visibility
from pipelines import PipelineDescriptor, load_pipeline_catalog
from pipeline_engine import PipelineDAG, PipelineExecutionPlan

from ..widgets import Tooltip


class PipelineLibraryController:
    def __init__(self, app) -> None:
        self.app = app

    def build_tab(self, parent: ttk.Frame) -> None:
        parent.columnconfigure(0, weight=1)
        parent.rowconfigure(2, weight=1)

        ttk.Label(parent, text="Select pipeline targets.").grid(
            row=0, column=0, sticky="w"
        )
        self._build_controls(parent)
        self._build_library_container(parent)

    def register(self) -> None:
        available, missing = load_pipeline_catalog()
        rows = sorted(
            [
                pipeline
                for pipeline in [*available, *missing]
                if getattr(pipeline, "visibility", "visible") != "hidden"
            ],
            key=lambda pipeline: pipeline.name.lower(),
        )
        self.app.pipeline_registry = {p.name: p for p in available}
        self.app.pipeline_catalog = {p.name: p for p in rows}

        try:
            self.app.pipeline_dag = PipelineDAG(rows)
        except (RuntimeError, ValueError) as exc:
            self.app.pipeline_dag = None
            self.app.settings_controller.show_settings_warning(
                "Pipeline DAG error",
                f"Pipeline dependency graph is invalid:\n{exc}",
            )
        else:
            rows = list(self.app.pipeline_dag.ordered_descriptors)

        self.app.pipeline_rows = rows
        self.sync_visibility(rows)
        self.populate(rows)

    def selected_target_pipeline_names(self) -> list[str]:
        return [
            pipeline.name
            for pipeline in self.app.pipeline_rows
            if pipeline.available
            and self.app.pipeline_visibility.get(pipeline.name, False)
        ]

    def resolve_plan(self, target_names: list[str]) -> PipelineExecutionPlan:
        dag = self.app.pipeline_dag
        if dag is None:
            dag = PipelineDAG(self.app.pipeline_rows)
            self.app.pipeline_dag = dag
        return dag.resolve_targets(target_names)

    def select_all(self) -> None:
        self.set_all_visibility(True)

    def deselect_all(self) -> None:
        self.set_all_visibility(False)

    def bind_vertical_mousewheel(self, widget: tk.Misc, canvas: tk.Canvas) -> None:
        for sequence in ("<MouseWheel>", "<Button-4>", "<Button-5>"):
            widget.bind(
                sequence,
                lambda event, target_canvas=canvas: self.on_vertical_mousewheel(
                    event, target_canvas
                ),
                add="+",
            )

    def on_vertical_mousewheel(self, event: tk.Event, canvas: tk.Canvas) -> str | None:
        scroll_units = mousewheel_scroll_units(event)
        if not scroll_units:
            return None
        canvas.yview_scroll(scroll_units, "units")
        return "break"

    def populate(self, rows: list[PipelineDescriptor]) -> None:
        for child in self.app.pipeline_library_inner.winfo_children():
            child.destroy()
        self.app.pipeline_visibility_vars = {}
        self.app.pipeline_row_widgets = {}
        self.app.pipeline_library_inner.columnconfigure(0, weight=0)
        self.app.pipeline_library_inner.columnconfigure(1, weight=1)
        self.app.pipeline_library_inner.columnconfigure(2, weight=0)
        self._build_header()

        for idx, pipeline in enumerate(rows, start=1):
            self._build_pipeline_row(idx, pipeline)

        self.update_summary()

    def sync_visibility(self, rows: list[PipelineDescriptor]) -> None:
        visibility, changed = normalize_pipeline_visibility(
            (pipeline.name for pipeline in rows),
            self.app.settings_store.load_pipeline_visibility(),
        )
        for pipeline in rows:
            if not pipeline.available and visibility.get(pipeline.name, False):
                visibility[pipeline.name] = False
                changed = True
        self.app.pipeline_visibility = visibility
        if changed:
            self.persist_visibility()

    def persist_visibility(self) -> None:
        try:
            self.app.settings_store.save_pipeline_visibility(
                self.app.pipeline_visibility
            )
        except OSError as exc:
            self.app.settings_controller.show_settings_warning(
                "Settings not saved",
                f"Could not save pipeline selection preferences:\n{exc}",
            )

    def set_visibility(self, name: str, visible: bool) -> None:
        pipeline = self.app.pipeline_catalog.get(name)
        if pipeline is not None and not pipeline.available:
            visible = False
        if self.app.pipeline_visibility.get(name) == visible:
            return
        self.app.pipeline_visibility[name] = visible
        var = self.app.pipeline_visibility_vars.get(name)
        if var is not None and var.get() != visible:
            var.set(visible)
        self.persist_visibility()
        self.update_summary()

    def set_all_visibility(self, visible: bool) -> None:
        changed = False
        target_values = {
            pipeline.name: visible and pipeline.available
            for pipeline in self.app.pipeline_rows
        }
        for name, target_value in target_values.items():
            if self.app.pipeline_visibility.get(name) != target_value:
                self.app.pipeline_visibility[name] = target_value
                changed = True
        if not changed:
            return
        for name, var in self.app.pipeline_visibility_vars.items():
            var.set(self.app.pipeline_visibility.get(name, False))
        self.persist_visibility()
        self.update_summary()

    def update_summary(self) -> None:
        target_names = self.selected_target_pipeline_names()
        available_count = sum(
            1 for pipeline in self.app.pipeline_rows if pipeline.available
        )
        if not target_names:
            self.app.pipeline_library_summary_var.set(f"Targets: 0/{available_count}")
            return

        try:
            plan = self.resolve_plan(target_names)
        except (RuntimeError, ValueError):
            self.app.pipeline_library_summary_var.set("DAG error")
            return

        self.app.pipeline_library_summary_var.set(
            f"Targets: {len(target_names)}/{available_count} | "
            f"Run steps: {len(plan.descriptors)}"
        )

    def _build_controls(self, parent: ttk.Frame) -> None:
        controls = ttk.Frame(parent)
        controls.grid(row=1, column=0, sticky="ew", pady=(8, 4))
        controls.columnconfigure(2, weight=1)
        ttk.Button(
            controls,
            text="Select all",
            command=self.select_all,
        ).grid(row=0, column=0, sticky="w")
        ttk.Button(
            controls,
            text="Deselect all",
            command=self.deselect_all,
        ).grid(row=0, column=1, sticky="w", padx=(4, 0))
        ttk.Label(
            controls,
            textvariable=self.app.pipeline_library_summary_var,
        ).grid(row=0, column=2, sticky="e")

    def _build_library_container(self, parent: ttk.Frame) -> None:
        library_container = ttk.Frame(parent)
        library_container.grid(row=2, column=0, sticky="nsew", padx=(10, 10))
        library_container.columnconfigure(0, weight=1)
        library_container.rowconfigure(0, weight=1)

        library_bg = self.app._style.lookup("TFrame", "background") or self.app._bg_color
        self.app.pipeline_library_canvas = tk.Canvas(
            library_container, highlightthickness=0, bg=library_bg
        )
        self.app.pipeline_library_canvas.grid(row=0, column=0, sticky="nsew")
        library_scroll = ttk.Scrollbar(
            library_container,
            orient="vertical",
            command=self.app.pipeline_library_canvas.yview,
        )
        library_padding = tk.Frame(library_container, width=10, bg=library_bg)
        library_padding.grid(row=0, column=1, sticky="ns")
        library_scroll.grid(row=0, column=2, sticky="ns")
        self.app.pipeline_library_canvas.configure(yscrollcommand=library_scroll.set)
        self.app.pipeline_library_inner = ttk.Frame(self.app.pipeline_library_canvas)
        self.app.pipeline_library_window = (
            self.app.pipeline_library_canvas.create_window(
                (0, 0), window=self.app.pipeline_library_inner, anchor="nw"
            )
        )
        self._bind_library_canvas(library_scroll)

    def _bind_library_canvas(self, library_scroll: ttk.Scrollbar) -> None:
        self.app.pipeline_library_inner.bind(
            "<Configure>",
            lambda _evt: self.app.pipeline_library_canvas.configure(
                scrollregion=self.app.pipeline_library_canvas.bbox("all")
            ),
        )
        self.app.pipeline_library_canvas.bind(
            "<Configure>",
            lambda evt: self.app.pipeline_library_canvas.itemconfigure(
                self.app.pipeline_library_window, width=evt.width
            ),
        )
        for widget in (
            self.app.pipeline_library_canvas,
            self.app.pipeline_library_inner,
            library_scroll,
        ):
            self.bind_vertical_mousewheel(widget, self.app.pipeline_library_canvas)

    def _build_header(self) -> None:
        selected_header = ttk.Label(self.app.pipeline_library_inner, text="Target")
        selected_header.grid(row=0, column=0, sticky="w", pady=(0, 6))
        order_header = ttk.Label(self.app.pipeline_library_inner, text="Pipeline")
        order_header.grid(row=0, column=1, sticky="w", padx=(12, 0), pady=(0, 6))
        status_header = ttk.Label(self.app.pipeline_library_inner, text="Status")
        status_header.grid(row=0, column=2, sticky="e", padx=(12, 0), pady=(0, 6))
        for widget in (selected_header, order_header, status_header):
            self.bind_vertical_mousewheel(widget, self.app.pipeline_library_canvas)

    def _build_pipeline_row(self, idx: int, pipeline: PipelineDescriptor) -> None:
        is_available = getattr(pipeline, "available", True)
        var = tk.BooleanVar(
            value=self.app.pipeline_visibility.get(pipeline.name, False)
            and is_available
        )
        check = ttk.Checkbutton(
            self.app.pipeline_library_inner,
            text="",
            variable=var,
            state="normal" if is_available else "disabled",
            command=lambda name=pipeline.name, visible_var=var: (
                self.set_visibility(name, visible_var.get())
            ),
        )
        check.grid(row=idx, column=0, sticky="w", pady=(0, 6))
        name_label = ttk.Label(self.app.pipeline_library_inner, text=pipeline.name)
        name_label.grid(row=idx, column=1, sticky="w", padx=(12, 0), pady=(0, 6))
        status = ttk.Label(
            self.app.pipeline_library_inner,
            text=pipeline_status_text(pipeline),
        )
        status.grid(row=idx, column=2, sticky="e", padx=(12, 0), pady=(0, 6))
        self._bind_row_widgets(check, name_label, status)
        self._bind_pipeline_row_toggle(
            pipeline.name,
            var,
            name_label,
            status,
        )
        self._add_tooltips(pipeline, check, name_label, status)
        self.app.pipeline_row_widgets[pipeline.name] = name_label
        self.app.pipeline_visibility_vars[pipeline.name] = var

    def _bind_row_widgets(self, *widgets: tk.Misc) -> None:
        for widget in widgets:
            self.bind_vertical_mousewheel(widget, self.app.pipeline_library_canvas)

    def _bind_pipeline_row_toggle(
        self,
        name: str,
        var: tk.BooleanVar,
        *widgets: tk.Misc,
    ) -> None:
        def toggle(_event: tk.Event) -> str:
            pipeline = self.app.pipeline_catalog.get(name)
            if pipeline is not None and not pipeline.available:
                return "break"
            var.set(not var.get())
            self.set_visibility(name, var.get())
            return "break"

        for widget in widgets:
            widget.bind("<Button-1>", toggle, add="+")

    def _add_tooltips(
        self,
        pipeline: PipelineDescriptor,
        *widgets: tk.Misc,
    ) -> None:
        tip_text = descriptor_tooltip_text(pipeline)
        if not tip_text:
            return
        for widget in widgets:
            Tooltip(widget, tip_text, bg=self.app._surface_color, fg=self.app._text_fg)


def mousewheel_scroll_units(event: tk.Event) -> int:
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


def descriptor_tooltip_text(descriptor: PipelineDescriptor) -> str:
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


def pipeline_status_text(pipeline: PipelineDescriptor) -> str:
    if pipeline.available:
        return "Available"
    if pipeline.missing_deps:
        return f"Missing deps: {', '.join(pipeline.missing_deps)}"
    return "Unavailable"
