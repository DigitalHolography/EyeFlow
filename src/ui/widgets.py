from __future__ import annotations

import tkinter as tk
from collections.abc import Callable


class Tooltip:
    """Lightweight tooltip that shows on hover."""

    def __init__(
        self,
        widget: tk.Widget,
        text: str | Callable[[], str],
        bg: str = "#333333",
        fg: str = "#f7f7f7",
    ) -> None:
        self.widget = widget
        self.text = text
        self.bg = bg
        self.fg = fg
        self.tipwindow: tk.Toplevel | None = None
        widget.bind("<Enter>", self._show)
        widget.bind("<Leave>", self._hide)

    def _resolved_text(self) -> str:
        if callable(self.text):
            try:
                return str(self.text())
            except Exception:
                return ""
        return str(self.text)

    def _show(self, _event=None) -> None:
        text = self._resolved_text().strip()
        if self.tipwindow or not text:
            return
        x = self.widget.winfo_rootx() + 24
        y = self.widget.winfo_rooty() + self.widget.winfo_height() + 6
        self.tipwindow = tw = tk.Toplevel(self.widget)
        tw.wm_overrideredirect(True)
        tw.wm_geometry(f"+{x}+{y}")
        label = tk.Label(
            tw,
            text=text,
            justify="left",
            background=self.bg,
            foreground=self.fg,
            relief="solid",
            borderwidth=1,
            wraplength=360,
            padx=8,
            pady=6,
        )
        label.pack()

    def _hide(self, _event=None) -> None:
        if self.tipwindow:
            self.tipwindow.destroy()
            self.tipwindow = None
