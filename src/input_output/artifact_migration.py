"""Rename existing EyeFlow artifacts using the shared writer naming contract.

Usage: python -m input_output.artifact_migration RESULTS_PARENT [--dry-run]
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

from .writers.artifact_names import prefixed_artifact_path


def plan_artifact_renames(output_dir: str | Path) -> tuple[tuple[Path, Path], ...]:
    """Check a complete output folder before moving any files."""
    root = Path(output_dir).expanduser().resolve()
    if not root.is_dir() or not root.name.endswith("_EF"):
        raise ValueError(f"Expected an existing <stem>_EF output directory: {root}")
    stem = root.name[:-3]
    planned = []
    destinations = set()
    for source in sorted(root.rglob("*")):
        if not source.is_file():
            continue
        destination = prefixed_artifact_path(source, stem)
        if destination == source:
            continue
        if not source.resolve().is_relative_to(root) or not destination.resolve().is_relative_to(
            root
        ):
            raise ValueError(f"Artifact path escapes its output directory: {source}")
        key = os.path.normcase(str(destination))
        if destination.exists() or key in destinations:
            raise FileExistsError(
                f"Artifact rename would overwrite an existing file: {destination}"
            )
        destinations.add(key)
        planned.append((source, destination))
    return tuple(planned)


def migrate_artifacts(
    output_dir: str | Path, *, dry_run: bool = False
) -> tuple[tuple[Path, Path], ...]:
    """Rename only files; retain main HDF5, directory layout, and file contents."""
    planned = plan_artifact_renames(output_dir)
    if dry_run:
        return planned
    moved = []
    try:
        for source, destination in planned:
            source.rename(destination)
            moved.append((source, destination))
    except OSError:
        for source, destination in reversed(moved):
            destination.rename(source)
        raise
    return planned


def find_output_dirs(parent: str | Path) -> tuple[Path, ...]:
    """Find acquisition results, pruning input companions and unrelated caches."""
    root = Path(parent).expanduser().resolve()
    if not root.is_dir():
        raise NotADirectoryError(root)
    if root.name.endswith("_EF"):
        return (root,)
    found = []
    for directory, children, _files in os.walk(root):
        current = Path(directory)
        children[:] = [
            name
            for name in children
            if name not in {".git", ".venv", "__pycache__"} and not name.endswith(("_HD", "_DV"))
        ]
        if current.name.endswith("_EF"):
            stem = current.name[:-3]
            if (current / "h5" / f"{stem}_EF.h5").is_file():
                found.append(current)
            children[:] = []
    return tuple(sorted(found))


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Prefix existing EyeFlow artifact filenames with acquisition stems."
    )
    parser.add_argument("results", type=Path, help="An EyeFlow output folder or its parent.")
    parser.add_argument(
        "--dry-run", action="store_true", help="List renames without changing files."
    )
    args = parser.parse_args()
    output_dirs = find_output_dirs(args.results)
    # Validate every run before starting the batch.
    for output_dir in output_dirs:
        plan_artifact_renames(output_dir)
    total = 0
    for output_dir in output_dirs:
        changes = migrate_artifacts(output_dir, dry_run=args.dry_run)
        for source, destination in changes:
            print(f"{source} -> {destination.name}")
        total += len(changes)
    print(
        f"{'Planned' if args.dry_run else 'Renamed'} {total} files in {len(output_dirs)} output folders."
    )


if __name__ == "__main__":
    main()
