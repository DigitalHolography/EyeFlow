"""Minimal-mode widget construction."""

from __future__ import annotations

import tkinter as tk
import tkinter.font as tkfont
from tkinter import ttk


def build_minimal_view(app, parent: ttk.Frame) -> None:
    parent.columnconfigure(0, weight=1)
    parent.rowconfigure(0, weight=1)
    parent.grid_anchor("n")

    minimal_wraplength = 480
    content = ttk.Frame(parent, padding=(24, 24, 24, 24))
    content.grid(row=0, column=0, pady=(8, 12))
    content.columnconfigure(0, minsize=minimal_wraplength)
    app.minimal_content = content

    app.minimal_title_label = ttk.Label(
        content,
        text="EyeFlow",
        font=_get_minimal_title_font(app),
    )
    app.minimal_title_label.grid(row=0, column=0, pady=(0, 10))
    _build_logo(app, content)
    _build_input_controls(app, content, minimal_wraplength)
    _build_progress(app, content, minimal_wraplength)


def _build_logo(app, content: ttk.Frame) -> None:
    minimal_logo = app.resource_controller.load_scaled_logo_image(
        max_width=360, max_height=144
    )
    if minimal_logo is None:
        return
    app._minimal_logo_image = minimal_logo
    app.minimal_logo_label = ttk.Label(content, image=minimal_logo)
    app.minimal_logo_label.grid(row=1, column=0, pady=(0, 18))


def _build_input_controls(
    app,
    content: ttk.Frame,
    minimal_wraplength: int,
) -> None:
    app.minimal_holo_browse_button = ttk.Button(
        content,
        text="Select .holo or .txt",
        command=app.run_controller.choose_holo_file,
    )
    app.minimal_holo_browse_button.grid(row=2, column=0, pady=(0, 10))
    app.minimal_holo_input_path_label = tk.Label(
        content,
        textvariable=app.minimal_holo_input_path_var,
        bg=app._bg_color,
        fg=app._muted_fg,
        justify="center",
        wraplength=minimal_wraplength,
    )
    app.minimal_holo_input_path_label.grid(
        row=3, column=0, pady=(0, 4), sticky="ew"
    )
    _build_status_row(app, content)
    app.minimal_run_button = ttk.Button(
        content,
        text="Run",
        command=app.run_controller.run_process,
    )
    app.minimal_run_button.grid(row=5, column=0, pady=(0, 18))


def _build_status_row(app, content: ttk.Frame) -> None:
    app.minimal_holo_status_frame = ttk.Frame(content)
    app.minimal_holo_status_frame.grid(row=4, column=0, pady=(0, 18))
    app.minimal_holo_hd_status_label = tk.Label(
        app.minimal_holo_status_frame,
        textvariable=app.holo_hd_status_var,
        bg=app._bg_color,
        fg=app._muted_fg,
        justify="center",
    )
    app.minimal_holo_hd_status_label.pack(side="left")
    tk.Label(
        app.minimal_holo_status_frame,
        text=" | ",
        bg=app._bg_color,
        fg=app._muted_fg,
        justify="center",
    ).pack(side="left")
    app.minimal_holo_dv_status_label = tk.Label(
        app.minimal_holo_status_frame,
        textvariable=app.holo_dv_status_var,
        bg=app._bg_color,
        fg=app._muted_fg,
        justify="center",
    )
    app.minimal_holo_dv_status_label.pack(side="left")
    if not app.holo_hd_status_var.get() and not app.holo_dv_status_var.get():
        app.minimal_holo_status_frame.grid_remove()


def _build_progress(
    app,
    content: ttk.Frame,
    minimal_wraplength: int,
) -> None:
    app.minimal_progress = ttk.Progressbar(
        content,
        orient="horizontal",
        mode="determinate",
        maximum=100,
        variable=app.run_progress_var,
        length=340,
        style=app._progress_primary_style,
    )
    app.minimal_progress.grid(row=6, column=0, sticky="ew")
    app.minimal_status_label = tk.Label(
        content,
        textvariable=app.minimal_status_var,
        bg=app._bg_color,
        fg=app._text_fg,
        justify="center",
        wraplength=minimal_wraplength,
    )
    app.minimal_status_label.grid(row=7, column=0, pady=(8, 0), sticky="ew")


def _get_minimal_title_font(app) -> tkfont.Font:
    if app._minimal_title_font is None:
        title_font = tkfont.nametofont("TkDefaultFont").copy()
        base_size = int(title_font.cget("size")) or 10
        title_font.configure(size=base_size * 2)
        app._minimal_title_font = title_font
    return app._minimal_title_font
