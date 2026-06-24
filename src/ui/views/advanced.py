"""Advanced-mode widget construction."""

from __future__ import annotations

import tkinter as tk
from tkinter import ttk

from ..widgets import Tooltip


def build_advanced_view(app, parent: ttk.Frame) -> None:
    parent.columnconfigure(0, weight=1)
    parent.rowconfigure(0, weight=1)
    app.notebook = ttk.Notebook(parent)
    app.notebook.grid(row=0, column=0, sticky="nsew")

    app.run_tab = ttk.Frame(app.notebook, padding=10)
    app.pipeline_library_tab = ttk.Frame(app.notebook, padding=10)
    app.notebook.add(app.run_tab, text="Run")
    app.notebook.add(app.pipeline_library_tab, text="Pipeline Library")
    _build_run_tab(app, app.run_tab)
    app.pipeline_library_controller.build_tab(app.pipeline_library_tab)


def _build_run_tab(app, parent: ttk.Frame) -> None:
    parent.columnconfigure(1, weight=1)
    parent.columnconfigure(2, weight=0)
    parent.rowconfigure(4, weight=1)
    _build_input_row(app, parent)
    _build_run_controls(app, parent)
    _build_run_log(app, parent)


def _build_input_row(app, parent: ttk.Frame) -> None:
    input_label = ttk.Label(parent, text="Input (.holo or .txt)")
    input_label.grid(row=0, column=0, sticky="w", padx=(0, 8))
    Tooltip(
        input_label,
        app.input_controller.reference_holo_tooltip_text,
        bg=app._surface_color,
        fg=app._text_fg,
    )
    app.holo_input_entry = ttk.Entry(parent, textvariable=app.holo_input_var)
    app.holo_input_entry.grid(row=0, column=1, sticky="ew", padx=(0, 4))
    app.holo_browse_button = ttk.Button(
        parent,
        text="Select files",
        command=app.run_controller.choose_holo_file,
    )
    app.holo_browse_button.grid(row=0, column=2, sticky="w")
    _build_status_row(app, parent)


def _build_status_row(app, parent: ttk.Frame) -> None:
    app.holo_status_frame = ttk.Frame(parent)
    app.holo_status_frame.grid(
        row=1, column=1, columnspan=2, sticky="w", pady=(2, 0)
    )
    app.holo_hd_status_label = tk.Label(
        app.holo_status_frame,
        textvariable=app.holo_hd_status_var,
        bg=app._bg_color,
        fg=app._muted_fg,
        justify="left",
        anchor="w",
    )
    app.holo_hd_status_label.pack(side="left")
    tk.Label(
        app.holo_status_frame,
        text=" | ",
        bg=app._bg_color,
        fg=app._muted_fg,
        justify="left",
        anchor="w",
    ).pack(side="left")
    app.holo_dv_status_label = tk.Label(
        app.holo_status_frame,
        textvariable=app.holo_dv_status_var,
        bg=app._bg_color,
        fg=app._muted_fg,
        justify="left",
        anchor="w",
    )
    app.holo_dv_status_label.pack(side="left")
    if not app.holo_hd_status_var.get() and not app.holo_dv_status_var.get():
        app.holo_status_frame.grid_remove()


def _build_run_controls(app, parent: ttk.Frame) -> None:
    controls = ttk.Frame(parent)
    controls.grid(row=2, column=0, columnspan=3, sticky="ew", pady=(12, 4))
    app.advanced_run_button = ttk.Button(
        controls,
        text="Run",
        command=app.run_controller.run_process,
    )
    app.advanced_run_button.grid(row=0, column=0, sticky="w")


def _build_run_log(app, parent: ttk.Frame) -> None:
    ttk.Label(parent, text="Run log").grid(
        row=3, column=0, sticky="nw", pady=(8, 2)
    )
    run_log_frame = ttk.Frame(parent)
    run_log_frame.grid(row=4, column=0, columnspan=3, sticky="nsew")
    run_log_frame.columnconfigure(0, weight=1)
    run_log_frame.rowconfigure(0, weight=1)
    app.run_log = tk.Text(
        run_log_frame,
        height=14,
        state="disabled",
        bg=app._text_bg,
        fg=app._text_fg,
        insertbackground=app._text_fg,
    )
    run_log_scroll = ttk.Scrollbar(
        run_log_frame,
        orient="vertical",
        command=app.run_log.yview,
    )
    app.run_log.configure(yscrollcommand=run_log_scroll.set)
    app.run_log.grid(row=0, column=0, sticky="nsew")
    run_log_scroll.grid(row=0, column=1, sticky="ns")
