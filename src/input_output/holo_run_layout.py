"""Path layout for one selected HOLO run."""

from dataclasses import dataclass
from pathlib import Path

from .schema import DOPPLER_VIEW_SCHEMA, H5SourceSchema, HOLODOPPLER_SCHEMA

HDF5_SUFFIXES = (".h5", ".hdf5")
INPUT_SCHEMAS = (HOLODOPPLER_SCHEMA, DOPPLER_VIEW_SCHEMA)


@dataclass(frozen=True)
class HoloRunLayout:
    _holo_path: Path
    _stem: str
    _root_dir: Path

    @classmethod
    def from_holo(
        cls,
        holo_path: str | Path,
        *,
        output_root: str | Path | None = None,
    ) -> "HoloRunLayout":
        path = _absolute(Path(holo_path).expanduser())
        root = _absolute(Path(output_root).expanduser()) if output_root else path.parent
        return cls(_holo_path=path, _stem=path.stem, _root_dir=root / path.stem)

    @property
    def holo_path(self) -> Path:
        return self._holo_path

    @property
    def stem(self) -> str:
        return self._stem

    @property
    def root_dir(self) -> Path:
        return self._root_dir

    @property
    def ef_dir(self) -> Path:
        return self._run_dir("EF")

    @property
    def hd_h5(self) -> Path:
        return self._input_h5(HOLODOPPLER_SCHEMA)

    @property
    def dv_h5(self) -> Path:
        return self._input_h5(DOPPLER_VIEW_SCHEMA)

    @property
    def has_hd_h5(self) -> bool:
        return self._has_h5(HOLODOPPLER_SCHEMA)

    @property
    def has_dv_h5(self) -> bool:
        return self._has_h5(DOPPLER_VIEW_SCHEMA)

    def with_root_dir(self, root_dir: str | Path) -> "HoloRunLayout":
        return HoloRunLayout(
            _holo_path=self._holo_path,
            _stem=self._stem,
            _root_dir=_absolute(Path(root_dir).expanduser()),
        )

    def require_inputs(self) -> None:
        if not self.root_dir.is_dir():
            raise FileNotFoundError(f"Could not find data folder:\n{self.root_dir}")
        errors: list[str] = []
        for schema in INPUT_SCHEMAS:
            try:
                self._input_h5(schema)
            except FileNotFoundError as exc:
                errors.append(str(exc))
        if errors:
            raise FileNotFoundError(
                "Missing required input data for the selected .holo file:\n\n"
                + "\n\n".join(errors)
            )

    def _run_dir(self, suffix: str) -> Path:
        return self._root_dir / f"{self._stem}_{suffix}"

    def _h5_dir(self, schema: H5SourceSchema) -> Path:
        return self._run_dir(schema.companion_suffix) / schema.h5_folder_name

    def _preferred_h5(self, schema: H5SourceSchema) -> Path:
        folder_name = f"{self._stem}_{schema.companion_suffix}"
        filename = schema.h5_filename_template.format(
            stem=self._stem,
            folder=folder_name,
            companion=schema.companion_suffix,
        )
        return self._h5_dir(schema) / filename

    def _input_h5(self, schema: H5SourceSchema) -> Path:
        files = self._h5_files(schema)
        if not files:
            raise FileNotFoundError(
                f"{schema.label} HDF5 file missing in expected folder:\n"
                f"{self._h5_dir(schema)}"
            )
        preferred = self._preferred_h5(schema)
        if preferred in files:
            return preferred
        if len(files) == 1:
            return files[0]
        candidates = "\n".join(str(path) for path in files)
        raise FileNotFoundError(
            f"Multiple {schema.label} HDF5 files found in:\n{self._h5_dir(schema)}\n\n"
            f"Expected one file, preferably named:\n{preferred.name}\n\n"
            f"Candidates:\n{candidates}"
        )

    def _has_h5(self, schema: H5SourceSchema) -> bool:
        return bool(self._h5_files(schema))

    def _h5_files(self, schema: H5SourceSchema) -> list[Path]:
        folder = self._h5_dir(schema)
        if not folder.is_dir():
            return []
        return sorted(
            path
            for path in folder.iterdir()
            if path.is_file() and path.suffix.lower() in HDF5_SUFFIXES
        )


def _absolute(path: Path) -> Path:
    return path if path.is_absolute() else Path.cwd() / path
