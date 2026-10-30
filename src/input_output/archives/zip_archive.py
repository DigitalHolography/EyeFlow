"""Create and extract zip archives used by EyeFlow."""

from __future__ import annotations

import shutil
import zipfile
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from pathlib import Path
from uuid import uuid4

from ..writers.artifact_names import prefixed_artifact_path


@contextmanager
def extracted_zip_tree(zip_path: str | Path) -> Iterator[Path]:
    from pipeline_engine.storage import temporary_directory

    with temporary_directory() as tmp_dir:
        target_root = Path(tmp_dir).resolve()
        with zipfile.ZipFile(zip_path, "r") as archive:
            _safe_extract_all(archive, target_root)
        yield target_root


def _safe_extract_all(archive: zipfile.ZipFile, target_root: Path) -> None:
    """Extract an archive without allowing members to escape the target root."""

    for member in archive.infolist():
        member_path = Path(member.filename.replace("\\", "/"))
        if member_path.is_absolute() or ".." in member_path.parts:
            raise ValueError(f"Unsafe ZIP member path: {member.filename}")
        target = (target_root / member_path).resolve()
        try:
            target.relative_to(target_root)
        except ValueError as exc:
            raise ValueError(f"Unsafe ZIP member path: {member.filename}") from exc
        if member.is_dir():
            target.mkdir(parents=True, exist_ok=True)
            continue
        target.parent.mkdir(parents=True, exist_ok=True)
        with archive.open(member) as source, target.open("wb") as destination:
            shutil.copyfileobj(source, destination)


def create_zip_from_tree(
    tree_root: str | Path,
    zip_path: str | Path,
    *,
    progress_callback: Callable[[int, int, Path], None] | None = None,
    stem: str | None = None,
) -> Path:
    tree_root_path = Path(tree_root).expanduser().resolve()
    zip_path_obj = prefixed_artifact_path(zip_path, stem) if stem is not None else Path(zip_path)
    zip_path_obj.parent.mkdir(parents=True, exist_ok=True)

    files = sorted(
        (path for path in tree_root_path.rglob("*") if path.is_file()),
        key=lambda path: path.relative_to(tree_root_path).as_posix(),
    )

    with zipfile.ZipFile(
        zip_path_obj,
        "w",
        compression=zipfile.ZIP_DEFLATED,
        compresslevel=1,
    ) as archive:
        total_files = len(files)
        if progress_callback is not None:
            progress_callback(0, total_files, Path("."))
        for idx, file_path in enumerate(files, start=1):
            archive.write(file_path, file_path.relative_to(tree_root_path))
            if progress_callback is not None:
                progress_callback(
                    idx,
                    total_files,
                    file_path.relative_to(tree_root_path),
                )
    return zip_path_obj


def write_zip_artifact(
    tree_root: str | Path,
    path: str | Path,
    *,
    stem: str | None = None,
    progress_callback: Callable[[int, int, Path], None] | None = None,
) -> Path:
    """Name and publish an archive, retaining the previous file on failure."""
    target = prefixed_artifact_path(path, stem) if stem is not None else Path(path)
    staging = target.with_name(f".{target.name}.eyeflow-staging-{uuid4().hex}")
    try:
        create_zip_from_tree(tree_root, staging, progress_callback=progress_callback)
        staging.replace(target)
    finally:
        if staging.exists():
            staging.unlink()
    return target


__all__ = [
    "create_zip_from_tree",
    "extracted_zip_tree",
    "write_zip_artifact",
]
