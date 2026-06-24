"""Controller for UI resource lookup and image loading."""

from __future__ import annotations

import sys
import tkinter as tk
from pathlib import Path


class ResourceController:
    def __init__(self, app) -> None:
        self.app = app

    def resource_roots(self) -> list[Path]:
        roots: list[Path] = []
        frozen_root = getattr(sys, "_MEIPASS", None)
        if frozen_root:
            roots.append(Path(frozen_root))
        roots.append(Path(__file__).resolve().parents[2])
        roots.append(Path.cwd())
        return roots

    def resolve_logo_path(self) -> Path | None:
        for root in self.resource_roots():
            candidate = root / "EyeFlow_logo.png"
            if candidate.is_file():
                return candidate
        return None

    def load_logo_image(self) -> tk.PhotoImage | None:
        logo_path = self.resolve_logo_path()
        if logo_path is None:
            return None
        try:
            return tk.PhotoImage(file=str(logo_path))
        except tk.TclError:
            return None

    def load_scaled_logo_image(
        self,
        *,
        max_width: int,
        max_height: int,
    ) -> tk.PhotoImage | None:
        image = self.load_logo_image()
        if image is None:
            return None

        scale_x = max(1, (image.width() + max_width - 1) // max_width)
        scale_y = max(1, (image.height() + max_height - 1) // max_height)
        scale = max(scale_x, scale_y)
        if scale > 1:
            image = image.subsample(scale, scale)
        return image

    def set_window_icon(self) -> None:
        image = self.load_logo_image()
        if image is None:
            return
        self.app._window_icon_image = image
        try:
            self.app.iconphoto(True, self.app._window_icon_image)
        except tk.TclError:
            pass
