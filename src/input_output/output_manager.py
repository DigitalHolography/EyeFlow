"""App-scoped output manager for one EyeFlow run."""

from dataclasses import dataclass
from enum import Enum
from pathlib import Path
import shutil

from .inputs import HoloRunLayout


class OutputType(Enum):
    H5 = "h5"
    PNG = "png"
    AVI = "avi"
    PDF = "pdf"
    EPS = "eps"


@dataclass(frozen=True)
class OutputManager:
    layout: HoloRunLayout

    def prepare(self, *, replace: bool = False) -> None:
        output_dir = self.layout.ef_dir
        if replace and output_dir.exists():
            _reset_output_dir(output_dir)
        else:
            output_dir.mkdir(parents=True, exist_ok=True)

    def dir_for(self, output_type: OutputType) -> Path:
        return self.layout.ef_dir / output_type.value

    def path_for(
        self,
        output_type: OutputType,
        filename: str | None = None,
    ) -> Path:
        return self.dir_for(output_type) / self._filename_for(output_type, filename)

    def open_h5(self, filename: str | None = None, mode: str = "w"):
        from .writers.h5 import open_h5

        path = self.path_for(OutputType.H5, filename)
        path.parent.mkdir(parents=True, exist_ok=True)
        return open_h5(path, mode)

    def write_png(self, output, filename: str | None = None) -> Path:
        from .writers.png import write_png_file

        return write_png_file(self.path_for(OutputType.PNG, filename), output)

    def report_path(self) -> Path:
        """Return the canonical PDF report path for this run."""
        return self.path_for(OutputType.PDF, f"{self.layout.stem}_report.pdf")

    def _filename_for(self, output_type: OutputType, filename: str | None) -> str:
        if filename:
            return filename
        if output_type is OutputType.H5:
            return f"{self.layout.stem}_EF.h5"
        return self.layout.stem


def output_stem(output) -> str:
    """Return the run stem exposed by an output namespace or manager."""
    manager = getattr(output, "manager", None)
    layout = getattr(manager, "layout", None)
    stem = getattr(layout, "stem", None)
    if stem is None:
        stem = getattr(getattr(output, "layout", None), "stem", None)
    return str(stem or "eyeflow")


def artifact_path(
    output,
    output_type: OutputType,
    suffix: str,
    *,
    stem: str | None = None,
    subfolder: str | None = None,
) -> Path:
    """Resolve and prepare a stem-prefixed artifact path consistently."""
    filename = f"{stem or output_stem(output)}_{suffix}"
    if subfolder:
        filename = f"{subfolder}/{filename}"
    path = output.path_for(output_type, filename)
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def _reset_output_dir(path: str | Path) -> None:
    """Replace one run output directory, preserving the existing error message."""
    path_obj = Path(path)
    try:
        if path_obj.exists():
            if path_obj.is_dir():
                shutil.rmtree(path_obj)
            else:
                path_obj.unlink()
            if path_obj.exists():
                raise OSError(f"Output path still exists after removal: {path_obj}")
        path_obj.mkdir(parents=True, exist_ok=False)
    except OSError as exc:
        raise RuntimeError(
            "Could not replace the existing output directory. Close any File Explorer "
            f"window, terminal, or application using this folder, then retry:\n{path_obj}"
        ) from exc
