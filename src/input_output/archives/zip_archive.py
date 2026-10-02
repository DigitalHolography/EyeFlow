"""Create and extract zip archives used by EyeFlow."""

from __future__ import annotations

from collections.abc import Callable, Iterator
from contextlib import contextmanager
from pathlib import Path
from tempfile import TemporaryDirectory
import shutil
import zipfile


@contextmanager
def extracted_zip_tree(zip_path: str | Path) -> Iterator[Path]:
    with TemporaryDirectory() as tmp_dir:
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
) -> None:
    tree_root_path = Path(tree_root).expanduser().resolve()
    zip_path_obj = Path(zip_path)
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


def reset_output_dir(path: str | Path) -> None:
    path_obj = Path(path)
    try:
        _remove_existing_output_path(path_obj)
        path_obj.mkdir(parents=True, exist_ok=False)
    except OSError as exc:
        raise RuntimeError(_locked_output_dir_message(path_obj)) from exc


def _remove_existing_output_path(path_obj: Path) -> None:
    if not path_obj.exists():
        return
    if path_obj.is_dir():
        shutil.rmtree(path_obj)
    else:
        path_obj.unlink()
    if path_obj.exists():
        raise OSError(f"Output path still exists after removal: {path_obj}")


def _locked_output_dir_message(path_obj: Path) -> str:
    return (
        "Could not replace the existing output directory. Close any File Explorer "
        f"window, terminal, or application using this folder, then retry:\n{path_obj}"
    )


__all__ = [
    "create_zip_from_tree",
    "extracted_zip_tree",
    "reset_output_dir",
]
