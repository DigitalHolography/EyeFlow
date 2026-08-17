#!/usr/bin/env python3
"""Double-click Windows interface for the EyeFlow waveform fixture extractor."""

from __future__ import annotations

import argparse
import os
import queue
import threading
import tkinter as tk
from dataclasses import dataclass
from pathlib import Path
from tkinter import filedialog, messagebox, ttk
from typing import Callable, Sequence

from extract_waveform_fixture import (
    DEFAULT_MAX_OUTPUT_MB,
    FixtureExtractionError,
    discover_h5_files,
    extract_fixture,
    prepare_fixture,
)


APP_NAME = "EyeFlow Waveform Fixture Extractor"
APP_VERSION = "0.1.0"


@dataclass(frozen=True)
class JobSummary:
    discovered: int
    succeeded: int
    failed: int
    created: int
    reused: int
    dry_run: bool


def run_extraction_job(
    input_path: Path,
    output_dir: Path,
    *,
    recursive: bool,
    dry_run: bool,
    log: Callable[[str], None],
) -> JobSummary:
    """Run one validation/extraction job without depending on Tk state."""
    files = discover_h5_files(
        [input_path],
        recursive=recursive,
        output_dir=output_dir,
    )
    if not files:
        raise FixtureExtractionError("No HDF5 files were found in the selected input.")

    succeeded = 0
    failed = 0
    created = 0
    reused = 0
    for index, source_path in enumerate(files, start=1):
        try:
            if dry_run:
                prepared = prepare_fixture(
                    source_path,
                    max_output_mb=DEFAULT_MAX_OUTPUT_MB,
                )
                log(
                    f"Input {index}/{len(files)} is compatible "
                    f"({len(prepared.arrays)} datasets, "
                    f"{prepared.estimated_bytes / 1024:.1f} KiB)."
                )
            else:
                result = extract_fixture(
                    source_path,
                    output_dir,
                    max_output_mb=DEFAULT_MAX_OUTPUT_MB,
                )
                status = "Created" if result.created else "Reused"
                log(
                    f"{status} {result.output_path.name} "
                    f"({result.dataset_count} datasets, "
                    f"{result.estimated_bytes / 1024:.1f} KiB)."
                )
                created += int(result.created)
                reused += int(not result.created)
            succeeded += 1
        except FixtureExtractionError as exc:
            log(f"Input {index}/{len(files)} failed: {exc}")
            failed += 1

    return JobSummary(
        discovered=len(files),
        succeeded=succeeded,
        failed=failed,
        created=created,
        reused=reused,
        dry_run=dry_run,
    )


class ExtractorWindow:
    def __init__(self, root: tk.Tk) -> None:
        self.root = root
        self.events: queue.Queue[tuple[str, object]] = queue.Queue()
        self.worker: threading.Thread | None = None

        self.input_var = tk.StringVar()
        self.output_var = tk.StringVar()
        self.recursive_var = tk.BooleanVar(value=True)
        self.status_var = tk.StringVar(value="Choose an EyeFlow H5 file or folder.")

        self._configure_window()
        self._build_layout()
        self.root.after(100, self._drain_events)

    def _configure_window(self) -> None:
        self.root.title(f"{APP_NAME} {APP_VERSION}")
        self.root.geometry("780x560")
        self.root.minsize(680, 480)
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)

    def _build_layout(self) -> None:
        outer = ttk.Frame(self.root, padding=18)
        outer.grid(row=0, column=0, sticky="nsew")
        outer.columnconfigure(0, weight=1)
        outer.rowconfigure(6, weight=1)

        title = ttk.Label(
            outer,
            text=APP_NAME,
            font=("Segoe UI", 16, "bold"),
        )
        title.grid(row=0, column=0, sticky="w")
        ttk.Label(
            outer,
            text=(
                "Create compact waveform-only H5 fixtures without copying source "
                "filenames, images, masks, maps, or acquisition metadata."
            ),
            wraplength=720,
        ).grid(row=1, column=0, sticky="ew", pady=(4, 16))

        input_frame = ttk.LabelFrame(outer, text="Input", padding=10)
        input_frame.grid(row=2, column=0, sticky="ew")
        input_frame.columnconfigure(0, weight=1)
        ttk.Entry(input_frame, textvariable=self.input_var).grid(
            row=0,
            column=0,
            columnspan=2,
            sticky="ew",
        )
        ttk.Button(
            input_frame,
            text="Choose H5 file",
            command=self._choose_input_file,
        ).grid(row=1, column=0, sticky="w", pady=(8, 0))
        ttk.Button(
            input_frame,
            text="Choose folder",
            command=self._choose_input_folder,
        ).grid(row=1, column=1, sticky="w", padx=(8, 0), pady=(8, 0))
        ttk.Checkbutton(
            input_frame,
            text="Search subfolders",
            variable=self.recursive_var,
        ).grid(row=1, column=2, sticky="e", padx=(16, 0), pady=(8, 0))

        output_frame = ttk.LabelFrame(outer, text="Output folder", padding=10)
        output_frame.grid(row=3, column=0, sticky="ew", pady=(12, 0))
        output_frame.columnconfigure(0, weight=1)
        ttk.Entry(output_frame, textvariable=self.output_var).grid(
            row=0,
            column=0,
            sticky="ew",
        )
        ttk.Button(
            output_frame,
            text="Choose folder",
            command=self._choose_output_folder,
        ).grid(row=0, column=1, padx=(8, 0))
        self.open_output_button = ttk.Button(
            output_frame,
            text="Open output",
            command=self._open_output_folder,
        )
        self.open_output_button.grid(row=0, column=2, padx=(8, 0))

        actions = ttk.Frame(outer)
        actions.grid(row=4, column=0, sticky="ew", pady=(14, 0))
        self.validate_button = ttk.Button(
            actions,
            text="Validate inputs",
            command=lambda: self._start_job(dry_run=True),
        )
        self.validate_button.grid(row=0, column=0)
        self.extract_button = ttk.Button(
            actions,
            text="Extract fixtures",
            command=lambda: self._start_job(dry_run=False),
        )
        self.extract_button.grid(row=0, column=1, padx=(8, 0))
        self.progress = ttk.Progressbar(actions, mode="indeterminate", length=190)
        self.progress.grid(row=0, column=2, padx=(18, 0))

        ttk.Label(
            outer,
            textvariable=self.status_var,
        ).grid(row=5, column=0, sticky="w", pady=(10, 6))

        log_frame = ttk.LabelFrame(outer, text="Results", padding=8)
        log_frame.grid(row=6, column=0, sticky="nsew")
        log_frame.columnconfigure(0, weight=1)
        log_frame.rowconfigure(0, weight=1)
        self.log_text = tk.Text(
            log_frame,
            height=12,
            wrap="word",
            state="disabled",
            font=("Consolas", 9),
        )
        scrollbar = ttk.Scrollbar(
            log_frame,
            orient="vertical",
            command=self.log_text.yview,
        )
        self.log_text.configure(yscrollcommand=scrollbar.set)
        self.log_text.grid(row=0, column=0, sticky="nsew")
        scrollbar.grid(row=0, column=1, sticky="ns")

        ttk.Label(
            outer,
            text=(
                "Privacy notice: extracted waveforms remain patient-derived and "
                "potentially sensitive. Follow your institution's data-governance rules."
            ),
            wraplength=720,
            foreground="#8a3b12",
        ).grid(row=7, column=0, sticky="ew", pady=(12, 0))

    def _choose_input_file(self) -> None:
        selected = filedialog.askopenfilename(
            title="Choose an EyeFlow H5 file",
            filetypes=(("HDF5 files", "*.h5 *.hdf5"), ("All files", "*.*")),
        )
        if selected:
            self.input_var.set(selected)
            if not self.output_var.get():
                self.output_var.set(str(Path(selected).parent / "WaveformFixtures"))

    def _choose_input_folder(self) -> None:
        selected = filedialog.askdirectory(title="Choose the EyeFlow output folder")
        if selected:
            self.input_var.set(selected)
            if not self.output_var.get():
                self.output_var.set(str(Path(selected) / "WaveformFixtures"))

    def _choose_output_folder(self) -> None:
        selected = filedialog.askdirectory(title="Choose the fixture output folder")
        if selected:
            self.output_var.set(selected)

    def _open_output_folder(self) -> None:
        output = Path(self.output_var.get().strip()).expanduser()
        if not output.exists():
            messagebox.showinfo(APP_NAME, "The output folder does not exist yet.")
            return
        os.startfile(output)  # type: ignore[attr-defined]

    def _start_job(self, *, dry_run: bool) -> None:
        if self.worker is not None and self.worker.is_alive():
            return
        try:
            input_path, output_dir = self._validated_paths()
        except FixtureExtractionError as exc:
            messagebox.showerror(APP_NAME, str(exc))
            return

        if not dry_run:
            approved = messagebox.askokcancel(
                "Create patient-derived fixtures?",
                (
                    "The output removes source filenames and acquisition metadata, "
                    "but the waveforms remain patient-derived and potentially sensitive.\n\n"
                    "Create the fixtures now?"
                ),
            )
            if not approved:
                return

        self._set_running(True)
        action = "validation" if dry_run else "extraction"
        self.status_var.set(f"Running {action}…")
        self._append_log(f"Started {action}.")
        recursive = bool(self.recursive_var.get())
        self.worker = threading.Thread(
            target=self._worker_job,
            args=(input_path, output_dir, recursive, dry_run),
            daemon=True,
        )
        self.worker.start()

    def _validated_paths(self) -> tuple[Path, Path]:
        input_text = self.input_var.get().strip()
        output_text = self.output_var.get().strip()
        if not input_text:
            raise FixtureExtractionError("Choose an input H5 file or folder.")
        if not output_text:
            raise FixtureExtractionError("Choose an output folder.")
        input_path = Path(input_text).expanduser()
        output_dir = Path(output_text).expanduser()
        if not input_path.exists():
            raise FixtureExtractionError("The selected input does not exist.")
        return input_path, output_dir

    def _worker_job(
        self,
        input_path: Path,
        output_dir: Path,
        recursive: bool,
        dry_run: bool,
    ) -> None:
        try:
            summary = run_extraction_job(
                input_path,
                output_dir,
                recursive=recursive,
                dry_run=dry_run,
                log=lambda message: self.events.put(("log", message)),
            )
            self.events.put(("done", summary))
        except Exception as exc:
            self.events.put(("error", str(exc)))

    def _drain_events(self) -> None:
        while True:
            try:
                event, payload = self.events.get_nowait()
            except queue.Empty:
                break
            if event == "log":
                self._append_log(str(payload))
            elif event == "done":
                self._finish_job(payload)
            elif event == "error":
                self._set_running(False)
                self.status_var.set("Operation failed.")
                self._append_log(f"Failed: {payload}")
                messagebox.showerror(APP_NAME, str(payload))
        self.root.after(100, self._drain_events)

    def _finish_job(self, payload: object) -> None:
        if not isinstance(payload, JobSummary):
            self._set_running(False)
            messagebox.showerror(APP_NAME, "Unexpected extraction result.")
            return
        self._set_running(False)
        if payload.dry_run:
            summary_text = (
                f"Validated {payload.succeeded} of {payload.discovered} input(s); "
                f"{payload.failed} failed."
            )
        else:
            summary_text = (
                f"Created {payload.created}, reused {payload.reused}, "
                f"failed {payload.failed}."
            )
        self.status_var.set(summary_text)
        self._append_log(summary_text)
        if payload.failed:
            messagebox.showwarning(APP_NAME, summary_text)
        else:
            messagebox.showinfo(APP_NAME, summary_text)

    def _set_running(self, running: bool) -> None:
        state = "disabled" if running else "normal"
        self.validate_button.configure(state=state)
        self.extract_button.configure(state=state)
        if running:
            self.progress.start(10)
        else:
            self.progress.stop()

    def _append_log(self, message: str) -> None:
        self.log_text.configure(state="normal")
        self.log_text.insert("end", message + "\n")
        self.log_text.see("end")
        self.log_text.configure(state="disabled")


def run_self_test(input_path: Path, output_dir: Path) -> int:
    """Exercise the packaged extraction path without opening a window."""
    try:
        summary = run_extraction_job(
            input_path,
            output_dir,
            recursive=False,
            dry_run=False,
            log=lambda _message: None,
        )
    except Exception:
        return 1
    return 0 if summary.succeeded == 1 and summary.failed == 0 else 1


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--self-test", type=Path)
    parser.add_argument("--output-dir", type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.self_test is not None:
        if args.output_dir is None:
            return 2
        return run_self_test(args.self_test, args.output_dir)

    root = tk.Tk()
    ExtractorWindow(root)
    root.mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
