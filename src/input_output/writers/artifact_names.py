"""Acquisition naming shared by all artifact formats."""

from pathlib import Path


def prefixed_artifact_path(path: str | Path, stem: str) -> Path:
    """Prefix the basename exactly once, preserving directories and extensions."""
    if not stem or stem in {".", ".."} or "/" in stem or "\\" in stem:
        raise ValueError("An artifact acquisition stem must be a nonempty filename component.")
    target = Path(str(path).replace("\\", "/"))
    if not target.name:
        raise ValueError("An artifact filename must not be empty.")
    prefix = f"{stem}_"
    if target.name.startswith(prefix):
        return target
    return target.with_name(prefix + target.name)


def acquisition_stem(output, fallback: str = "eyeflow") -> str:
    """Resolve the actual run stem from an output namespace or manager."""
    manager = getattr(output, "manager", None) or output
    stem = getattr(getattr(manager, "layout", None), "stem", None)
    return str(stem if stem is not None else fallback)


def labeled_artifact_path(path: str | Path, stem: str, label: str) -> Path:
    """Preserve legacy writer labels/subfolders, then add the acquisition prefix."""
    candidate = Path(str(path).replace("\\", "/"))
    if label != stem:
        prefix = f"{stem}_"
        if candidate.name.startswith(prefix):
            candidate = candidate.with_name(candidate.name[len(prefix) :])
        candidate = f"{label}_{candidate}"
    return prefixed_artifact_path(candidate, stem)
