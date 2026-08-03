"""Controller for pipeline target selection and DAG summaries."""

from __future__ import annotations

import tkinter as tk
from tkinter import ttk

from app_settings import (
    normalize_pipeline_options,
    normalize_pipeline_visibility,
    runtime_pipelines_path,
)
from pipelines import PipelineDescriptor, load_pipeline_catalog
from pipeline_engine import PipelineDAG, PipelineExecutionPlan

from ..services import services_for
from ..widgets import Tooltip


class PipelineLibraryController:
    _STATUS_COLUMN_PADDING = (24, 24)

    def __init__(self, app) -> None:
        self.app = app

    def configure_library_columns(self, inner, *, row_count: int = 1) -> None:
        inner.columnconfigure(0, weight=1, minsize=180)
        inner.columnconfigure(1, weight=0, minsize=6)
        inner.columnconfigure(2, weight=0, minsize=0)
        divider_slot = ttk.Frame(inner, cursor="sb_h_double_arrow")
        divider_slot.grid(
            row=0,
            column=1,
            rowspan=max(1, row_count),
            sticky="nsew",
        )
        divider_slot.grid_propagate(False)
        divider_slot._library_resize_inner = inner
        for relx in (0.35, 0.65):
            handle_line = ttk.Separator(divider_slot, orient="vertical")
            handle_line._library_resize_inner = inner
            handle_line.place(
                relx=relx,
                rely=0.5,
                anchor="center",
                width=1,
                height=22,
            )
            handle_line.bind(
                "<ButtonPress-1>", self._start_library_column_resize
            )
            handle_line.bind("<B1-Motion>", self._resize_library_columns)
            handle_line.bind(
                "<ButtonRelease-1>", self._finish_library_column_resize
            )
        divider_slot.bind(
            "<ButtonPress-1>", self._start_library_column_resize
        )
        divider_slot.bind("<B1-Motion>", self._resize_library_columns)
        divider_slot.bind(
            "<ButtonRelease-1>", self._finish_library_column_resize
        )

    def _start_library_column_resize(self, event) -> None:
        inner = getattr(event.widget, "_library_resize_inner", event.widget.master)
        inner.update_idletasks()
        name_width = inner.grid_bbox(0, 0)[2]
        status_width = max(1, self._library_column_requested_width(inner, 2))
        self._library_resize_state = (
            inner,
            event.x_root,
            name_width,
            status_width,
        )

    def _resize_library_columns(self, event) -> None:
        state = getattr(self, "_library_resize_state", None)
        if state is None:
            return
        inner, start_x, start_width, status_min_width = state
        divider_width = 6
        new_width = start_width + event.x_root - start_x
        max_width = max(180, inner.winfo_width() - status_min_width - divider_width)
        new_width = min(max(180, new_width), max_width)
        inner.columnconfigure(0, weight=0, minsize=new_width)
        inner.columnconfigure(2, weight=1, minsize=status_min_width)

    @staticmethod
    def _library_column_requested_width(inner, column: int) -> int:
        widths = []
        for widget in inner.grid_slaves(column=column):
            grid_info = widget.grid_info()
            padx = grid_info.get("padx", 0)
            if isinstance(padx, (tuple, list)):
                padding = sum(int(value) for value in padx)
            elif isinstance(padx, str) and " " in padx:
                padding = sum(int(value) for value in padx.split())
            else:
                padding = int(padx or 0) * 2
            widths.append(widget.winfo_reqwidth() + padding)
        return max(widths, default=0)

    def _finish_library_column_resize(self, _event) -> None:
        self._library_resize_state = None

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
        catalog = sorted(
            [*available, *missing],
            key=lambda pipeline: pipeline.name.lower(),
        )
        self.app.pipeline_registry = {p.name: p for p in available}
        self.app.pipeline_catalog = {p.name: p for p in catalog}

        try:
            self.app.pipeline_dag = PipelineDAG(catalog)
        except (RuntimeError, ValueError) as exc:
            self.app.pipeline_dag = None
            rows = [
                pipeline
                for pipeline in catalog
                if getattr(pipeline, "visibility", "visible") != "hidden"
            ]
            self.app.settings_controller.show_settings_warning(
                "Pipeline DAG error",
                f"Pipeline dependency graph is invalid:\n{exc}",
            )
        else:
            rows = [
                pipeline
                for pipeline in self.app.pipeline_dag.ordered_descriptors
                if getattr(pipeline, "visibility", "visible") != "hidden"
            ]

        self.app.pipeline_rows = rows
        self.sync_visibility(rows)
        self.sync_options(rows)
        self.populate(rows)

    def selected_target_pipeline_names(self) -> list[str]:
        return [
            pipeline.name
            for pipeline in self.app.pipeline_rows
            if pipeline.available
            and self.app.pipeline_visibility.get(pipeline.name, False)
        ]

    def selected_pipeline_options(self) -> dict[str, tuple[str, ...]]:
        return {
            pipeline.name: tuple(
                option.name
                for option in pipeline.options
                if self.app.pipeline_option_visibility
                .get(pipeline.name, {})
                .get(option.name, option.default_enabled)
            )
            for pipeline in self.app.pipeline_rows
            if pipeline.options
        }

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

    def open_pipelines_folder(self) -> None:
        try:
            pipelines_folder = runtime_pipelines_path()
            pipelines_folder.mkdir(parents=True, exist_ok=True)
            services_for(self.app).path_opener.open_path(pipelines_folder)
        except OSError as exc:
            services_for(self.app).dialogs.showerror(
                "Folder unavailable",
                f"Could not open pipelines folder:\n{exc}",
            )

    def reload_pipelines(self) -> None:
        if getattr(self.app, "_pipeline_run_active", False):
            services_for(self.app).dialogs.showwarning(
                "Pipeline run active",
                "Wait for the current pipeline run to finish before reloading pipelines.",
            )
            return
        self.register()

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
        self.app.pipeline_option_vars = {}
        self.app.pipeline_option_widgets = {}
        self.app.pipeline_expand_buttons = {}
        row_count = 1 + sum(1 + len(pipeline.options) for pipeline in rows)
        self.configure_library_columns(
            self.app.pipeline_library_inner,
            row_count=row_count,
        )
        self._build_header()

        row_index = 1
        for pipeline in rows:
            row_index += self._build_pipeline_row(row_index, pipeline)

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

    def sync_options(self, rows: list[PipelineDescriptor]) -> None:
        selections, changed = normalize_pipeline_options(
            {
                pipeline.name: pipeline.options
                for pipeline in rows
                if pipeline.options
            },
            self.app.settings_store.load_pipeline_options(),
        )
        self.app.pipeline_option_visibility = selections
        if changed:
            self.persist_options()

    def persist_options(self) -> None:
        try:
            self.app.settings_store.save_pipeline_options(
                self.app.pipeline_option_visibility
            )
        except OSError as exc:
            self.app.settings_controller.show_settings_warning(
                "Settings not saved",
                f"Could not save pipeline option preferences:\n{exc}",
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
        self._update_option_widget_states(name)
        self.update_summary()

    def set_option_visibility(
        self,
        pipeline_name: str,
        option_name: str,
        enabled: bool,
    ) -> None:
        pipeline_values = self.app.pipeline_option_visibility.setdefault(
            pipeline_name,
            {},
        )
        if pipeline_values.get(option_name) == enabled:
            return
        pipeline_values[option_name] = enabled
        var = self.app.pipeline_option_vars.get(pipeline_name, {}).get(option_name)
        if var is not None and var.get() != enabled:
            var.set(enabled)
        self.persist_options()
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
        for name, var in self.app.pipeline_visibility_vars.items():
            var.set(self.app.pipeline_visibility.get(name, False))
            self._update_option_widget_states(name)
        options_changed = False
        for pipeline in self.app.pipeline_rows:
            values = self.app.pipeline_option_visibility.setdefault(
                pipeline.name,
                {},
            )
            for option in pipeline.options:
                if values.get(option.name) != visible:
                    values[option.name] = visible
                    options_changed = True
                option_var = self.app.pipeline_option_vars.get(
                    pipeline.name,
                    {},
                ).get(option.name)
                if option_var is not None:
                    option_var.set(visible)
        if changed:
            self.persist_visibility()
        if options_changed:
            self.persist_options()
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
        controls.columnconfigure(4, weight=1)
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
        ttk.Button(
            controls,
            text="Open folder",
            command=self.open_pipelines_folder,
        ).grid(row=0, column=2, sticky="w", padx=(4, 0))
        ttk.Button(
            controls,
            text="Reload",
            command=self.reload_pipelines,
        ).grid(row=0, column=3, sticky="w", padx=(4, 0))
        ttk.Label(
            controls,
            textvariable=self.app.pipeline_library_summary_var,
        ).grid(row=0, column=4, sticky="e")

    def _build_library_container(self, parent: ttk.Frame) -> None:
        library_container = ttk.Frame(parent)
        library_container.grid(row=2, column=0, sticky="nsew")
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
        library_scroll.grid(row=0, column=1, sticky="ns")
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
        status_header = ttk.Label(self.app.pipeline_library_inner, text="Status")
        status_header.grid(
            row=0,
            column=2,
            sticky="w",
            padx=self._STATUS_COLUMN_PADDING,
            pady=(0, 6),
        )
        for widget in (selected_header, status_header):
            self.bind_vertical_mousewheel(widget, self.app.pipeline_library_canvas)

    def _build_pipeline_row(self, idx: int, pipeline: PipelineDescriptor) -> int:
        is_available = getattr(pipeline, "available", True)
        var = tk.BooleanVar(
            value=self.app.pipeline_visibility.get(pipeline.name, False)
            and is_available
        )
        target_frame = ttk.Frame(self.app.pipeline_library_inner)
        target_frame.grid(row=idx, column=0, sticky="w", pady=(0, 6))
        target_column = 0
        if pipeline.options:
            expanded = self.app.pipeline_expanded.setdefault(pipeline.name, True)
            expand_button = ttk.Button(
                target_frame,
                text="-" if expanded else "+",
                width=2,
                command=lambda name=pipeline.name: (
                    self._toggle_pipeline_options(name)
                ),
            )
            expand_button.grid(row=0, column=0, sticky="w", padx=(0, 3))
            self.bind_vertical_mousewheel(
                expand_button,
                self.app.pipeline_library_canvas,
            )
            self.app.pipeline_expand_buttons[pipeline.name] = expand_button
            target_column = 1
        check = ttk.Checkbutton(
            target_frame,
            text=pipeline.name,
            variable=var,
            state="normal" if is_available else "disabled",
            command=lambda name=pipeline.name, visible_var=var: (
                self.set_visibility(name, visible_var.get())
            ),
        )
        check.grid(row=0, column=target_column, sticky="w")
        status = ttk.Label(
            self.app.pipeline_library_inner,
            text=pipeline_status_text(pipeline),
        )
        status.grid(
            row=idx,
            column=2,
            sticky="w",
            padx=self._STATUS_COLUMN_PADDING,
            pady=(0, 6),
        )
        self._bind_row_widgets(target_frame, check, status)
        self._bind_pipeline_row_toggle(
            pipeline.name,
            var,
            status,
        )
        self._add_tooltips(pipeline, check, status)
        self.app.pipeline_row_widgets[pipeline.name] = check
        self.app.pipeline_visibility_vars[pipeline.name] = var
        self.app.pipeline_option_vars[pipeline.name] = {}
        option_widgets: list[tk.Widget] = []
        expanded = self.app.pipeline_expanded.get(pipeline.name, True)
        for offset, option in enumerate(pipeline.options, start=1):
            option_var = tk.BooleanVar(
                value=self.app.pipeline_option_visibility
                .get(pipeline.name, {})
                .get(option.name, option.default_enabled)
            )
            option_check = ttk.Checkbutton(
                self.app.pipeline_library_inner,
                text=option.label,
                variable=option_var,
                command=lambda pipeline_name=pipeline.name,
                option_name=option.name,
                selected_var=option_var: self.set_option_visibility(
                    pipeline_name,
                    option_name,
                    selected_var.get(),
                ),
            )
            option_check.grid(
                row=idx + offset,
                column=0,
                sticky="w",
                padx=(34, 0),
                pady=(0, 4),
            )
            option_status = ttk.Label(
                self.app.pipeline_library_inner,
                text=option.description,
            )
            option_status.grid(
                row=idx + offset,
                column=2,
                sticky="w",
                padx=self._STATUS_COLUMN_PADDING,
                pady=(0, 4),
            )
            self._bind_row_widgets(option_check, option_status)
            if option.description:
                Tooltip(
                    option_check,
                    option.description,
                    bg=self.app._surface_color,
                    fg=self.app._text_fg,
                )
            self.app.pipeline_option_vars[pipeline.name][option.name] = option_var
            option_widgets.extend((option_check, option_status))
            if not expanded:
                option_check.grid_remove()
                option_status.grid_remove()
        self.app.pipeline_option_widgets[pipeline.name] = option_widgets
        self._update_option_widget_states(pipeline.name)
        return 1 + len(pipeline.options)

    def _toggle_pipeline_options(self, pipeline_name: str) -> None:
        expanded = not self.app.pipeline_expanded.get(pipeline_name, True)
        self.app.pipeline_expanded[pipeline_name] = expanded
        button = self.app.pipeline_expand_buttons.get(pipeline_name)
        if button is not None:
            button.configure(text="-" if expanded else "+")
        for widget in self.app.pipeline_option_widgets.get(pipeline_name, []):
            if expanded:
                widget.grid()
            else:
                widget.grid_remove()

    def _update_option_widget_states(self, pipeline_name: str) -> None:
        pipeline = self.app.pipeline_catalog.get(pipeline_name)
        enabled = bool(
            pipeline is not None
            and pipeline.available
            and self.app.pipeline_visibility.get(pipeline_name, False)
        )
        state = "normal" if enabled else "disabled"
        for widget in self.app.pipeline_option_widgets.get(pipeline_name, []):
            if isinstance(widget, ttk.Checkbutton):
                widget.configure(state=state)

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
