from __future__ import annotations

import tkinter as tk
import tkinter.font as tkfont
from tkinter import ttk

from .widgets import _Tooltip


class ViewBuilderMixin:
    def _build_ui(self) -> None:
        self._build_menu()

        container = ttk.Frame(self, padding=0)
        container.pack(fill="both", expand=True)
        self.main_container = container

        self.minimal_view = ttk.Frame(container, padding=0)
        self.advanced_view = ttk.Frame(container, padding=0)

        self._build_minimal_view(self.minimal_view)
        self._build_advanced_view(self.advanced_view)

    def _build_menu(self) -> None:
        self.ui_mode_var = tk.StringVar(value=self.ui_mode)
        menu_bar = tk.Menu(self)
        view_menu = tk.Menu(menu_bar, tearoff=False)
        view_menu.add_radiobutton(
            label="Minimal UI",
            value="minimal",
            variable=self.ui_mode_var,
            command=lambda: self._apply_ui_mode(self.ui_mode_var.get()),
        )
        view_menu.add_radiobutton(
            label="Advanced UI",
            value="advanced",
            variable=self.ui_mode_var,
            command=lambda: self._apply_ui_mode(self.ui_mode_var.get()),
        )
        menu_bar.add_cascade(label="View", menu=view_menu)
        self.configure(menu=menu_bar)

    def _build_minimal_view(self, parent: ttk.Frame) -> None:
        parent.columnconfigure(0, weight=1)
        parent.rowconfigure(0, weight=1)
        parent.grid_anchor("n")

        minimal_wraplength = 480

        content = ttk.Frame(parent, padding=(24, 24, 24, 24))
        content.grid(row=0, column=0, pady=(8, 12))
        content.columnconfigure(0, minsize=minimal_wraplength)
        self.minimal_content = content

        self.minimal_title_label = ttk.Label(
            content,
            text="EyeFlow",
            font=self._get_minimal_title_font(),
        )
        self.minimal_title_label.grid(row=0, column=0, pady=(0, 10))

        minimal_logo = self._load_scaled_logo_image(max_width=360, max_height=144)
        if minimal_logo is not None:
            self._minimal_logo_image = minimal_logo
            self.minimal_logo_label = ttk.Label(content, image=self._minimal_logo_image)
            self.minimal_logo_label.grid(row=1, column=0, pady=(0, 18))

        self.minimal_holo_browse_button = ttk.Button(
            content,
            text="Select .holo file",
            command=self.choose_holo_file,
        )
        self.minimal_holo_browse_button.grid(row=2, column=0, pady=(0, 10))
        self.minimal_holo_input_path_label = tk.Label(
            content,
            textvariable=self.minimal_holo_input_path_var,
            bg=self._bg_color,
            fg=self._muted_fg,
            justify="center",
            wraplength=minimal_wraplength,
        )
        self.minimal_holo_input_path_label.grid(
            row=3,
            column=0,
            pady=(0, 4),
            sticky="ew",
        )
        self.minimal_holo_status_frame = ttk.Frame(content)
        self.minimal_holo_status_frame.grid(row=4, column=0, pady=(0, 18))
        self.minimal_holo_hd_status_label = tk.Label(
            self.minimal_holo_status_frame,
            textvariable=self.holo_hd_status_var,
            bg=self._bg_color,
            fg=self._muted_fg,
            justify="center",
        )
        self.minimal_holo_hd_status_label.pack(side="left")
        self.minimal_holo_status_separator_label = tk.Label(
            self.minimal_holo_status_frame,
            text=" | ",
            bg=self._bg_color,
            fg=self._muted_fg,
            justify="center",
        )
        self.minimal_holo_status_separator_label.pack(side="left")
        self.minimal_holo_dv_status_label = tk.Label(
            self.minimal_holo_status_frame,
            textvariable=self.holo_dv_status_var,
            bg=self._bg_color,
            fg=self._muted_fg,
            justify="center",
        )
        self.minimal_holo_dv_status_label.pack(side="left")

        self.minimal_run_button = ttk.Button(
            content, text="Run", command=self.run_process
        )
        self.minimal_run_button.grid(row=5, column=0, pady=(0, 18))

        self.minimal_progress = ttk.Progressbar(
            content,
            orient="horizontal",
            mode="determinate",
            maximum=100,
            variable=self.run_progress_var,
            length=340,
            style=self._progress_primary_style,
        )
        self.minimal_progress.grid(row=6, column=0, sticky="ew")
        self.minimal_status_label = tk.Label(
            content,
            textvariable=self.minimal_status_var,
            bg=self._bg_color,
            fg=self._text_fg,
            justify="center",
            wraplength=minimal_wraplength,
        )
        self.minimal_status_label.grid(row=7, column=0, pady=(8, 0), sticky="ew")

    def _get_minimal_title_font(self) -> tkfont.Font:
        if self._minimal_title_font is None:
            title_font = tkfont.nametofont("TkDefaultFont").copy()
            base_size = int(title_font.cget("size")) or 10
            title_font.configure(size=base_size * 2)
            self._minimal_title_font = title_font
        return self._minimal_title_font

    def _build_advanced_view(self, parent: ttk.Frame) -> None:
        parent.columnconfigure(0, weight=1)
        parent.rowconfigure(0, weight=1)

        self.notebook = ttk.Notebook(parent)
        self.notebook.grid(row=0, column=0, sticky="nsew")

        self.run_tab = ttk.Frame(self.notebook, padding=10)
        self.pipeline_library_tab = ttk.Frame(self.notebook, padding=10)
        self.notebook.add(self.run_tab, text="Run")
        self.notebook.add(self.pipeline_library_tab, text="Pipeline Library")

        self._build_run_tab(self.run_tab)
        self._build_pipeline_library_tab(self.pipeline_library_tab)

    def _build_run_tab(self, parent: ttk.Frame) -> None:
        parent.columnconfigure(1, weight=1)
        parent.columnconfigure(2, weight=0)
        parent.rowconfigure(4, weight=1)

        input_label = ttk.Label(parent, text="Input (.holo)")
        input_label.grid(
            row=0,
            column=0,
            sticky="w",
            padx=(0, 8),
        )
        _Tooltip(
            input_label,
            self._reference_holo_tooltip_text,
            bg=self._surface_color,
            fg=self._text_fg,
        )
        self.holo_input_entry = ttk.Entry(parent, textvariable=self.holo_input_var)
        self.holo_input_entry.grid(row=0, column=1, sticky="ew", padx=(0, 4))
        self.holo_browse_button = ttk.Button(
            parent, text="Select file", command=self.choose_holo_file
        )
        self.holo_browse_button.grid(
            row=0,
            column=2,
            sticky="w",
        )
        self.holo_status_frame = ttk.Frame(parent)
        self.holo_status_frame.grid(
            row=1,
            column=1,
            columnspan=2,
            sticky="w",
            pady=(2, 0),
        )
        self.holo_hd_status_label = tk.Label(
            self.holo_status_frame,
            textvariable=self.holo_hd_status_var,
            bg=self._bg_color,
            fg=self._muted_fg,
            justify="left",
            anchor="w",
        )
        self.holo_hd_status_label.pack(side="left")
        self.holo_status_separator_label = tk.Label(
            self.holo_status_frame,
            text=" | ",
            bg=self._bg_color,
            fg=self._muted_fg,
            justify="left",
            anchor="w",
        )
        self.holo_status_separator_label.pack(side="left")
        self.holo_dv_status_label = tk.Label(
            self.holo_status_frame,
            textvariable=self.holo_dv_status_var,
            bg=self._bg_color,
            fg=self._muted_fg,
            justify="left",
            anchor="w",
        )
        self.holo_dv_status_label.pack(side="left")

        controls = ttk.Frame(parent)
        controls.grid(row=2, column=0, columnspan=3, sticky="ew", pady=(12, 4))

        run_btn = ttk.Button(controls, text="Run", command=self.run_process)
        run_btn.grid(row=0, column=0, sticky="w")

        ttk.Label(parent, text="Run log").grid(
            row=3, column=0, sticky="nw", pady=(8, 2)
        )
        run_log_frame = ttk.Frame(parent)
        run_log_frame.grid(row=4, column=0, columnspan=3, sticky="nsew")
        run_log_frame.columnconfigure(0, weight=1)
        run_log_frame.rowconfigure(0, weight=1)
        self.run_log = tk.Text(
            run_log_frame,
            height=14,
            state="disabled",
            bg=self._text_bg,
            fg=self._text_fg,
            insertbackground=self._text_fg,
        )
        run_log_scroll = ttk.Scrollbar(
            run_log_frame, orient="vertical", command=self.run_log.yview
        )
        self.run_log.configure(yscrollcommand=run_log_scroll.set)
        self.run_log.grid(row=0, column=0, sticky="nsew")
        run_log_scroll.grid(row=0, column=1, sticky="ns")
