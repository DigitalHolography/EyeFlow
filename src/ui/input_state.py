from __future__ import annotations

import re
from collections.abc import Sequence
from pathlib import Path

from input_output import HOLO_SUFFIX, holo_input_status, reset_output_dir


class InputStateMixin:
    def _on_holo_input_changed(self, *_args) -> None:
        if not self._synchronizing_holo_input_var:
            self._selected_holo_input_paths = []
        self._update_minimal_path_labels()

    def _path_from_var(self, raw_value: str) -> Path | None:
        value = raw_value.strip()
        if not value:
            return None
        path = Path(value).expanduser()
        if not path.is_absolute():
            path = Path.cwd() / path
        return path

    def _set_holo_input_var(self, value: str) -> None:
        self._synchronizing_holo_input_var = True
        try:
            self.holo_input_var.set(value)
        finally:
            self._synchronizing_holo_input_var = False

    def _selected_holo_paths(self) -> list[Path]:
        selected_path = self._selected_holo_path()
        return [selected_path] if selected_path is not None else []

    def _selected_holo_path(self) -> Path | None:
        if self._selected_holo_input_paths:
            return self._selected_holo_input_paths[0]
        return self._path_from_var(self.holo_input_var.get() or "")

    def _assign_holo_input_path(self, input_path: Path) -> None:
        normalized_path = input_path.expanduser()
        if not normalized_path.is_absolute():
            normalized_path = Path.cwd() / normalized_path

        self._selected_holo_input_paths = [normalized_path]
        self._set_holo_input_var(str(normalized_path))
        self._apply_input_defaults(normalized_path)

    def _assign_holo_input_paths(self, input_paths: Sequence[Path]) -> None:
        if not input_paths:
            self._selected_holo_input_paths = []
            self._set_holo_input_var("")
            self._apply_input_defaults(None)
            return
        self._assign_holo_input_path(input_paths[0])

    def _normalized_input_token(self, input_path: Path) -> str:
        token = re.sub(r"[^A-Za-z0-9]+", "_", input_path.stem).strip("_")
        return token or input_path.stem or "output"

    def _default_work_h5_name_for_input(self, input_path: Path | None) -> str:
        base_name = (
            self._normalized_input_token(input_path)
            if input_path is not None
            else "output"
        )
        return f"{base_name}_eyeflow.h5"

    def _default_work_h5_name(self) -> str:
        return self._default_work_h5_name_for_input(self._selected_holo_path())

    def _reference_holo_tooltip_text(self) -> str:
        return "Pick one reference .holo file."

    def _set_holo_status_parts(
        self,
        *,
        hd_text: str,
        hd_color: str,
        dv_text: str,
        dv_color: str,
    ) -> None:
        self.holo_hd_status_var.set(hd_text)
        self.holo_dv_status_var.set(dv_text)
        for label_name, color in (
            ("minimal_holo_hd_status_label", hd_color),
            ("minimal_holo_dv_status_label", dv_color),
            ("holo_hd_status_label", hd_color),
            ("holo_dv_status_label", dv_color),
        ):
            label = getattr(self, label_name, None)
            if label is not None:
                label.configure(fg=color)

    def _holo_data_status(
        self,
        holo_path: Path,
        *,
        require_holo_file: bool,
    ):
        return holo_input_status(
            holo_path,
            require_holo_file=require_holo_file,
        )

    def _update_minimal_found_statuses(self, holo_paths: Sequence[Path]) -> None:
        if not holo_paths:
            self._set_holo_status_parts(
                hd_text="HD waiting",
                hd_color=self._muted_fg,
                dv_text="DV waiting",
                dv_color=self._muted_fg,
            )
            return

        normalized_holo = holo_paths[0].expanduser()
        if not normalized_holo.is_absolute():
            normalized_holo = Path.cwd() / normalized_holo

        if (
            not normalized_holo.exists()
            or not normalized_holo.is_file()
            or normalized_holo.suffix.lower() != HOLO_SUFFIX
        ):
            self._set_holo_status_parts(
                hd_text="HD unavailable",
                hd_color=self._error_color,
                dv_text="DV unavailable",
                dv_color=self._error_color,
            )
            return

        status = self._holo_data_status(
            normalized_holo,
            require_holo_file=True,
        )

        self._set_holo_status_parts(
            hd_text="HD found" if status.hd else "HD not found",
            hd_color=self._success_color if status.hd else self._error_color,
            dv_text="DV found" if status.dv else "DV not found",
            dv_color=self._success_color if status.dv else self._error_color,
        )

    def _update_minimal_path_labels(self) -> None:
        holo_paths = self._selected_holo_paths()
        if not holo_paths:
            self.minimal_holo_input_path_var.set("No input selected")
        else:
            self.minimal_holo_input_path_var.set(str(holo_paths[0]))
        self._update_minimal_found_statuses(holo_paths)

    def _default_output_dir_for_input(self, input_path: Path) -> Path:
        output_dir = input_path.parent if input_path.is_file() else input_path
        if input_path.is_file() and input_path.suffix.lower() == HOLO_SUFFIX:
            output_dir = input_path.parent / input_path.stem / f"{input_path.stem}_EF"
        return output_dir

    def _prepare_default_output_dir(self, input_path: Path) -> Path:
        output_dir = self._default_output_dir_for_input(input_path)
        if output_dir.exists():
            if output_dir.is_dir():
                reset_output_dir(output_dir)
            else:
                output_dir.unlink()
        output_dir.mkdir(parents=True, exist_ok=True)
        return output_dir

    def _apply_input_defaults(self, input_path: Path | None) -> None:
        del input_path
        self._reset_progress()
        self._set_minimal_status("Ready.")
