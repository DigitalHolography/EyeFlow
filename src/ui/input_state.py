from collections.abc import Sequence
from pathlib import Path

from input_output import (
    INPUT_LIST_SUFFIX,
    default_output_dir_for_input,
    holo_input_status,
    read_holo_input_list,
    reset_output_dir,
    stem_input_status,
)
from input_output.output_manager import OutputManager


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
        if self._selected_holo_input_paths:
            return list(self._selected_holo_input_paths)
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
        normalized_paths: list[Path] = []
        for input_path in input_paths:
            path = input_path.expanduser()
            if not path.is_absolute():
                path = Path.cwd() / path
            normalized_paths.append(path)

        if not normalized_paths:
            self._selected_holo_input_paths = []
            self._set_holo_input_var("")
            self._apply_input_defaults(None)
            return

        self._selected_holo_input_paths = normalized_paths
        display_value = str(normalized_paths[0])
        if len(normalized_paths) > 1:
            display_value += f" (+{len(normalized_paths) - 1} more)"
        self._set_holo_input_var(display_value)
        self._apply_input_defaults(normalized_paths[0])

    def _reference_holo_tooltip_text(self) -> str:
        return "Pick one or more reference .holo files, or a .txt path list."

    def _set_shared_holo_status(
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

    def _update_holo_input_found_statuses(self, holo_paths: Sequence[Path]) -> None:
        if not holo_paths:
            self._set_shared_holo_status(
                hd_text="Awaiting HD",
                hd_color=self._muted_fg,
                dv_text="Awaiting DV",
                dv_color=self._muted_fg,
            )
            return
        if len(holo_paths) == 1 and holo_paths[0].suffix.lower() == INPUT_LIST_SUFFIX:
            self._update_input_list_found_statuses(holo_paths[0])
            return

        total_files = len(holo_paths)
        hd_found_count = 0
        dv_found_count = 0
        missing_hd: list[str] = []
        missing_dv: list[str] = []

        for holo_path in holo_paths:
            normalized_holo = holo_path.expanduser()
            if not normalized_holo.is_absolute():
                normalized_holo = Path.cwd() / normalized_holo

            status = self._holo_data_status(normalized_holo, require_holo_file=True)
            if status.hd:
                hd_found_count += 1
            else:
                missing_hd.append(normalized_holo.stem)
            if status.dv:
                dv_found_count += 1
            else:
                missing_dv.append(normalized_holo.stem)

        hd_text = _found_status_text("HD", hd_found_count, total_files, missing_hd)
        if missing_hd:
            hd_color = self._error_color
        else:
            hd_color = self._success_color

        dv_text = _found_status_text("DV", dv_found_count, total_files, missing_dv)
        if missing_dv:
            dv_color = self._error_color
        else:
            dv_color = self._success_color

        self._set_shared_holo_status(
            hd_text=hd_text,
            hd_color=hd_color,
            dv_text=dv_text,
            dv_color=dv_color,
        )

    def _update_input_list_found_statuses(self, input_list_path: Path) -> None:
        try:
            input_list = read_holo_input_list(input_list_path)
        except (OSError, ValueError) as exc:
            self._set_shared_holo_status(
                hd_text=str(exc),
                hd_color=self._error_color,
                dv_text="DV unavailable",
                dv_color=self._error_color,
            )
            return

        inputs = list(input_list.path_stem_pairs)
        stems = [stem for _root_dir, stem in inputs]
        statuses = [stem_input_status(stem, root_dir) for root_dir, stem in inputs]
        missing_hd = [stem for stem, status in zip(stems, statuses) if not status.hd]
        missing_dv = [stem for stem, status in zip(stems, statuses) if not status.dv]
        self._set_shared_holo_status(
            hd_text=_found_status_text(
                "HD", len(stems) - len(missing_hd), len(stems), missing_hd
            ),
            hd_color=self._error_color if missing_hd else self._success_color,
            dv_text=_found_status_text(
                "DV", len(stems) - len(missing_dv), len(stems), missing_dv
            ),
            dv_color=self._error_color if missing_dv else self._success_color,
        )

    def _update_minimal_path_labels(self) -> None:
        holo_paths = self._selected_holo_paths()
        if not holo_paths:
            self.minimal_holo_input_path_var.set("No input selected")
        elif len(holo_paths) == 1:
            self.minimal_holo_input_path_var.set(str(holo_paths[0]))
        else:
            self.minimal_holo_input_path_var.set(
                f"{holo_paths[0]} (+{len(holo_paths) - 1} more)"
            )
        self._update_holo_input_found_statuses(holo_paths)

    def _default_output_dir_for_input(self, input_path: Path) -> Path:
        return default_output_dir_for_input(input_path)

    def _prepare_default_output_dir(self, input_path: Path) -> Path:
        output_dir = self._default_output_dir_for_input(input_path)
        if output_dir.exists():
            if output_dir.is_dir():
                reset_output_dir(output_dir)
            else:
                output_dir.unlink()
        output_dir.mkdir(parents=True, exist_ok=True)
        return output_dir

    def _output_manager_for_input(self, input_path: Path) -> OutputManager:
        return OutputManager.from_holo(input_path)

    def _prepare_output_manager_for_input(self, input_path: Path) -> OutputManager:
        output_manager = self._output_manager_for_input(input_path)
        output_manager.prepare(replace=True)
        return output_manager

    def _apply_input_defaults(self, input_path: Path | None) -> None:
        del input_path
        self._reset_progress()
        self._set_minimal_status("Ready.")


def _found_status_text(
    label: str,
    found_count: int,
    total_count: int,
    missing_stems: Sequence[str],
) -> str:
    if total_count == 1:
        return f"{label} found" if found_count else f"{label} not found"
    text = f"Found {found_count}/{total_count} {label}"
    if missing_stems:
        text += ": missing " + ", ".join(missing_stems)
    return text
