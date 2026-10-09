"""App-scoped output manager for one EyeFlow run."""

from dataclasses import dataclass
from enum import Enum
from pathlib import Path

from .archives import reset_output_dir
from .holo_run_layout import HoloRunLayout
from .writers import open_h5, write_png_file
from .writers.artifact_names import prefixed_artifact_path


class OutputType(Enum):
    H5 = "h5"
    PNG = "png"
    MP4 = "mp4"
    AVI = "avi"
    PDF = "pdf"
    EPS = "eps"


@dataclass(frozen=True)
class OutputManager:
    layout: HoloRunLayout
    artifact_folder: str | None = None

    def for_artifact_namespace(self, namespace: str) -> "OutputManager":
        """Return a manager whose sidecars live below one variant namespace."""

        return OutputManager(self.layout, str(namespace))

    def for_workflow(self, method: str) -> "OutputManager":
        from .schema.eyeflow_output import VELOCITY_WORKFLOW_FOLDERS

        return self.for_artifact_namespace(VELOCITY_WORKFLOW_FOLDERS[method])

    def prepare(self, *, replace: bool = False) -> None:
        output_dir = self.layout.ef_dir
        if replace and output_dir.exists():
            reset_output_dir(output_dir)
        else:
            output_dir.mkdir(parents=True, exist_ok=True)

    def dir_for(self, output_type: OutputType) -> Path:
        directory = self.layout.ef_dir / output_type.value
        if self.artifact_folder is not None:
            directory /= self.artifact_folder
        return directory

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

    def write_png(self, output, filename: str | None = None) -> Path:
        return write_png_file(self.path_for(OutputType.PNG, filename), output)

    def _filename_for(self, output_type: OutputType, filename: str | None) -> str:
        if filename:
            return str(prefixed_artifact_path(filename, self.layout.stem))
        if output_type is OutputType.H5:
            return f"{self.layout.stem}_EF.h5"
        return f"{self.layout.stem}_output.{output_type.value}"
