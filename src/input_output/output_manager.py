"""App-scoped output manager for one EyeFlow run."""

from dataclasses import dataclass
from enum import Enum
from pathlib import Path

from .archives import reset_output_dir
from .holo_run_layout import HoloRunLayout
from .writers import open_h5, write_json_file


class OutputType(Enum):
    H5 = "h5"
    PNG = "png"
    JSON = "json"
    MP4 = "mp4"
    AVI = "avi"


@dataclass(frozen=True)
class OutputManager:
    layout: HoloRunLayout

    @classmethod
    def from_holo(
        cls,
        holo_path: str | Path,
        *,
        output_root: str | Path | None = None,
    ) -> "OutputManager":
        return cls(HoloRunLayout.from_holo(holo_path, output_root=output_root))

    def prepare(self, *, replace: bool = False) -> None:
        output_dir = self.layout.ef_dir
        if replace and output_dir.exists():
            reset_output_dir(output_dir)
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
        path = self.path_for(OutputType.H5, filename)
        path.parent.mkdir(parents=True, exist_ok=True)
        return open_h5(path, mode)

    def write(
        self,
        output,
        output_type: OutputType,
        filename: str | None = None,
    ) -> Path:
        if output_type is OutputType.JSON:
            return write_json_file(self.path_for(output_type, filename), output)
        if output_type is OutputType.H5:
            raise NotImplementedError("Use open_h5() for session-based H5 output.")
        raise NotImplementedError(
            f"Output writer for {output_type.value!r} is not implemented yet."
        )

    def _filename_for(self, output_type: OutputType, filename: str | None) -> str:
        if filename:
            return filename
        if output_type is OutputType.H5:
            return f"{self.layout.stem}_EF.h5"
        if output_type is OutputType.JSON:
            return f"{self.layout.stem}.json"
        return self.layout.stem
