"""Injectable UI services for dialogs, file pickers, and OS file opening."""

from __future__ import annotations

import os
import subprocess
import sys
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path
from tkinter import filedialog, messagebox
from typing import Protocol


class DialogService(Protocol):
    def showwarning(self, title: str, message: str) -> object: ...

    def showerror(self, title: str, message: str) -> object: ...

    def showinfo(self, title: str, message: str) -> object: ...


class FileDialogService(Protocol):
    def askopenfilenames(self, **options) -> Sequence[str]: ...


class PathOpener(Protocol):
    def open_path(self, path: Path) -> None: ...


class TkDialogService:
    def showwarning(self, title: str, message: str) -> object:
        return messagebox.showwarning(title, message)

    def showerror(self, title: str, message: str) -> object:
        return messagebox.showerror(title, message)

    def showinfo(self, title: str, message: str) -> object:
        return messagebox.showinfo(title, message)


class TkFileDialogService:
    def askopenfilenames(self, **options) -> Sequence[str]:
        return filedialog.askopenfilenames(**options)


class SystemPathOpener:
    def open_path(self, path: Path) -> None:
        if sys.platform == "win32":
            os.startfile(str(path))  # type: ignore[attr-defined]
            return

        opener = "open" if sys.platform == "darwin" else "xdg-open"
        subprocess.Popen(
            [opener, str(path)],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )


@dataclass(frozen=True)
class UiServices:
    dialogs: DialogService = field(default_factory=TkDialogService)
    file_dialogs: FileDialogService = field(default_factory=TkFileDialogService)
    path_opener: PathOpener = field(default_factory=SystemPathOpener)


def services_for(app) -> UiServices:
    return getattr(app, "ui_services", UiServices())
