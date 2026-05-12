"""Resolved output layout for one selected HOLO run."""

from dataclasses import dataclass
from pathlib import Path

from .schema import HOLO_SUFFIX


@dataclass(frozen=True)
class OutputLayout:
    holo_path: Path
    stem: str
    root_dir: Path

    @classmethod
    def from_holo(
        cls,
        holo_path: str | Path,
        *,
        output_root: str | Path | None = None,
    ) -> "OutputLayout":
        path = _absolute(Path(holo_path).expanduser())
        stem = path.stem
        root_dir = _root_dir_for_holo(path, output_root=output_root)
        return cls(holo_path=path, stem=stem, root_dir=root_dir)


def _root_dir_for_holo(
    holo_path: Path,
    *,
    output_root: str | Path | None,
) -> Path:
    if output_root is not None:
        return _absolute(Path(output_root).expanduser()) / f"{holo_path.stem}_EF"
    if holo_path.suffix.lower() == HOLO_SUFFIX:
        return holo_path.parent / holo_path.stem / f"{holo_path.stem}_EF"
    return holo_path.parent


def _absolute(path: Path) -> Path:
    return path if path.is_absolute() else Path.cwd() / path
