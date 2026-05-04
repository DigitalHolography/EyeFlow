from __future__ import annotations

import tkinter as tk
from collections.abc import Sequence
from tkinter import ttk

from app_settings import normalize_pipeline_order, normalize_pipeline_visibility
from pipelines import PipelineDescriptor, load_pipeline_catalog

from .widgets import _Tooltip


class PipelineLibraryMixin:
    def _build_pipeline_library_tab(self, parent: ttk.Frame) -> None:
        parent.columnconfigure(0, weight=1)
        parent.rowconfigure(2, weight=1)

        ttk.Label(
            parent,
            text="Select the pipelines to run. "
            "This preference is saved between app launches.",
        ).grid(row=0, column=0, sticky="w")

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
        rows = self._sync_pipeline_order(rows)
        self.pipeline_rows = rows
        self._sync_pipeline_visibility(rows)
        self._populate_pipeline_library(rows)
        self._install_drop_targets()

    def _descriptor_tooltip_text(self, descriptor) -> str:
        parts: list[str] = []
        description = getattr(descriptor, "description", "")
        if description:
            parts.append(description)
        required_pipelines = getattr(descriptor, "required_pipelines", [])
        if required_pipelines:
            parts.append(f"Requires pipelines: {', '.join(required_pipelines)}")
        missing_pipelines = getattr(descriptor, "missing_pipelines", [])
        if missing_pipelines:
            parts.append(
                "Unavailable until these pipelines are available: "
                f"{', '.join(missing_pipelines)}"
            )
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

        selected_header = ttk.Label(self.pipeline_library_inner, text="Selected")
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
            name_label = ttk.Label(
                row_frame,
                text=pipeline.name,
                cursor="fleur" if is_available else "",
            )
            name_label.grid(row=0, column=1, sticky="w", padx=(12, 0))
            self._bind_vertical_mousewheel(check, self.pipeline_library_canvas)
            self._bind_vertical_mousewheel(name_label, self.pipeline_library_canvas)
            self._bind_vertical_mousewheel(row_frame, self.pipeline_library_canvas)
            self.pipeline_row_widgets[pipeline.name] = row_frame

            if is_available:
                for widget in (row_frame, name_label):
                    self._bind_pipeline_drag(widget, pipeline.name)

            status_text = self._pipeline_status_text(pipeline)
            status = ttk.Label(row_frame, text=status_text)
            status.grid(row=0, column=2, sticky="w", padx=(12, 0))
            self._bind_vertical_mousewheel(status, self.pipeline_library_canvas)
            if is_available:
                self._bind_pipeline_drag(status, pipeline.name)

            tip_text = self._descriptor_tooltip_text(pipeline)
            if tip_text:
                _Tooltip(check, tip_text, bg=self._surface_color, fg=self._text_fg)
                _Tooltip(name_label, tip_text, bg=self._surface_color, fg=self._text_fg)
                _Tooltip(status, tip_text, bg=self._surface_color, fg=self._text_fg)

            self.pipeline_visibility_vars[pipeline.name] = var

        if (
            self._pipeline_drop_indicator is None
            or not self._pipeline_drop_indicator.winfo_exists()
        ):
            self._pipeline_drop_indicator = tk.Frame(
                self.pipeline_library_inner,
                bg=self._accent_color,
                height=2,
                highlightthickness=0,
                bd=0,
            )
        self._hide_pipeline_drop_indicator()
        self._update_pipeline_library_summary()

    def _bind_pipeline_drag(self, widget: tk.Widget, name: str) -> None:
        widget.bind(
            "<ButtonPress-1>",
            lambda event, pipeline_name=name: self._start_pipeline_drag(
                event,
                pipeline_name,
            ),
            add="+",
        )
        widget.bind(
            "<B1-Motion>",
            self._on_pipeline_drag_motion,
            add="+",
        )
        widget.bind(
            "<ButtonRelease-1>",
            self._finish_pipeline_drag,
            add="+",
        )

    def _sync_pipeline_order(
        self,
        rows: list[PipelineDescriptor],
    ) -> list[PipelineDescriptor]:
        ordered_names, changed = normalize_pipeline_order(
            (pipeline.name for pipeline in rows),
            self.settings_store.load_pipeline_order(),
        )
        rows_by_name = {pipeline.name: pipeline for pipeline in rows}
        ordered_rows = [
            rows_by_name[name] for name in ordered_names if name in rows_by_name
        ]
        if changed:
            self._persist_pipeline_order(ordered_rows)
        return ordered_rows

    def _persist_pipeline_order(
        self,
        rows: Sequence[PipelineDescriptor] | None = None,
    ) -> None:
        order = [pipeline.name for pipeline in (rows or self.pipeline_rows)]
        try:
            self.settings_store.save_pipeline_order(order)
        except OSError as exc:
            self._show_settings_warning(
                "Settings not saved",
                f"Could not save pipeline order preference:\n{exc}",
            )

    def _pipeline_index(self, name: str) -> int | None:
        return next(
            (
                idx
                for idx, pipeline in enumerate(self.pipeline_rows)
                if pipeline.name == name
            ),
            None,
        )

    def _move_pipeline_to_index(
        self,
        name: str,
        target_index: int,
        *,
        persist: bool = True,
        refresh: bool = True,
    ) -> bool:
        current_index = self._pipeline_index(name)
        if current_index is None:
            return False

        rows = list(self.pipeline_rows)
        pipeline = rows.pop(current_index)
        target_index = max(0, min(int(target_index), len(rows)))
        if target_index > current_index:
            target_index -= 1
        rows.insert(target_index, pipeline)
        if rows == self.pipeline_rows:
            return False

        self.pipeline_rows = rows
        if persist:
            self._persist_pipeline_order()
        if refresh:
            self._populate_pipeline_library(self.pipeline_rows)
        return True

    def _move_pipeline_to_top(
        self,
        name: str,
        *,
        persist: bool = True,
        refresh: bool = True,
    ) -> bool:
        return self._move_pipeline_to_index(
            name,
            0,
            persist=persist,
            refresh=refresh,
        )

    def _pipeline_drop_index(self, root_y: int) -> int:
        if not self.pipeline_rows:
            return 0
        self.pipeline_library_inner.update_idletasks()
        for idx, pipeline in enumerate(self.pipeline_rows):
            widget = self.pipeline_row_widgets.get(pipeline.name)
            if widget is None or not widget.winfo_exists():
                continue
            midpoint = widget.winfo_rooty() + (widget.winfo_height() / 2.0)
            if root_y < midpoint:
                return idx
        return len(self.pipeline_rows)

    def _pipeline_drop_indicator_y(self, drop_index: int) -> int:
        if not self.pipeline_rows:
            return 0
        clamped_index = max(0, min(drop_index, len(self.pipeline_rows)))
        if clamped_index >= len(self.pipeline_rows):
            last_widget = self.pipeline_row_widgets.get(self.pipeline_rows[-1].name)
            if last_widget is None or not last_widget.winfo_exists():
                return 0
            return int(last_widget.winfo_y() + last_widget.winfo_height())

        widget = self.pipeline_row_widgets.get(self.pipeline_rows[clamped_index].name)
        if widget is None or not widget.winfo_exists():
            return 0
        return int(widget.winfo_y())

    def _show_pipeline_drop_indicator(self, drop_index: int) -> None:
        indicator = getattr(self, "_pipeline_drop_indicator", None)
        if indicator is None:
            return
        self.pipeline_library_inner.update_idletasks()
        indicator_y = self._pipeline_drop_indicator_y(drop_index)
        indicator_width = max(self.pipeline_library_inner.winfo_width(), 1)
        indicator.place(
            x=0,
            y=max(indicator_y - 1, 0),
            width=indicator_width,
            height=2,
        )
        indicator.lift()

    def _hide_pipeline_drop_indicator(self) -> None:
        indicator = getattr(self, "_pipeline_drop_indicator", None)
        if indicator is not None:
            indicator.place_forget()

    def _start_pipeline_drag(self, event, name: str) -> str:
        self._dragging_pipeline_name = name
        self._dragging_pipeline_active = False
        self._drag_start_root_y = int(event.y_root)
        try:
            event.widget.grab_set()
        except tk.TclError:
            pass
        return "break"

    def _on_pipeline_drag_motion(self, event) -> str:
        if getattr(self, "_dragging_pipeline_name", None) is None:
            return "break"
        if not getattr(self, "_dragging_pipeline_active", False):
            if abs(int(event.y_root) - self._drag_start_root_y) < 4:
                return "break"
            self._dragging_pipeline_active = True
        drop_index = self._pipeline_drop_index(int(event.y_root))
        self._show_pipeline_drop_indicator(drop_index)
        return "break"

    def _finish_pipeline_drag(self, event) -> str:
        name = getattr(self, "_dragging_pipeline_name", None)
        self._dragging_pipeline_name = None
        was_active = getattr(self, "_dragging_pipeline_active", False)
        self._dragging_pipeline_active = False
        try:
            event.widget.grab_release()
        except tk.TclError:
            pass
        self._hide_pipeline_drop_indicator()
        if not name:
            return "break"
        if not was_active:
            return "break"
        target_index = self._pipeline_drop_index(int(event.y_root))
        self._move_pipeline_to_index(name, target_index)
        return "break"

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
        if visible:
            if self._move_pipeline_to_top(name):
                return
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

    def _update_pipeline_library_summary(self) -> None:
        selected_count = sum(
            1
            for pipeline in self.pipeline_rows
            if pipeline.available and self.pipeline_visibility.get(pipeline.name, False)
        )
        available_count = sum(
            1 for pipeline in self.pipeline_rows if pipeline.available
        )
        self.pipeline_library_summary_var.set(
            f"Selected: {selected_count}/{available_count}"
        )

    def select_all_pipelines(self) -> None:
        self._set_all_pipeline_visibility(True)

    def deselect_all_pipelines(self) -> None:
        self._set_all_pipeline_visibility(False)
