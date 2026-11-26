"""Local temporary directories shared by I/O and execution services."""

from __future__ import annotations

import os
import shutil
import tempfile
from dataclasses import dataclass
from pathlib import Path


def is_local_path(path: str | Path) -> bool:
    """Require a recognized local volume without searching other drives."""
    path = Path(path).resolve()
    if os.name == "nt":
        import ctypes

        # DRIVE_FIXED excludes UNC and mapped network drives.
        return ctypes.windll.kernel32.GetDriveTypeW(path.anchor) == 3

    import psutil

    try:
        mounts = sorted(
            psutil.disk_partitions(all=True), key=lambda mount: len(mount.mountpoint), reverse=True
        )
        for mount in mounts:
            if path.is_relative_to(Path(mount.mountpoint)):
                return mount.fstype.lower() in {
                    "ext2", "ext3", "ext4", "xfs", "btrfs", "zfs", "tmpfs", "overlay",
                    "apfs", "hfs", "hfs+", "hfsplus", "msdos", "vfat", "exfat",
                    "ntfs", "ntfs3", "f2fs", "ufs", "ufs2", "jfs", "reiserfs",
                }
    except (OSError, ValueError, psutil.Error):
        pass
    return False


def local_scratch_root() -> Path:
    """Use the configured local location or the OS temporary directory."""
    path = Path(os.environ.get("EYEFLOW_SCRATCH_DIR") or tempfile.gettempdir()).resolve()
    if not is_local_path(path):
        raise OSError(f"EyeFlow scratch must be local: {path}")
    try:
        path.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        raise OSError(f"Could not create EyeFlow scratch directory {path}: {exc}") from exc
    return path


def temporary_directory(prefix="eyeflow-", **kwargs):
    """Keep worker directories alive until the coordinator releases their files."""
    worker = os.environ.get("EYEFLOW_WORKER_DIRECTORY")
    if worker:
        root = Path(worker).resolve()
        return RetainedDirectory(tempfile.mkdtemp(prefix=prefix, dir=root), str(root))
    return tempfile.TemporaryDirectory(prefix=prefix, dir=local_scratch_root(), **kwargs)


@dataclass
class RetainedDirectory:
    """Explicitly owned worker directory, with no process-exit finalizer."""

    name: str
    root: str

    def cleanup(self):
        target, root = Path(self.name).resolve(), Path(self.root).resolve()
        if target == root or not target.is_relative_to(root):
            raise ValueError("Worker scratch cleanup escaped its owning directory.")
        try:
            shutil.rmtree(target)
        except FileNotFoundError:
            pass

    def __enter__(self):
        return self.name

    def __exit__(self, *exc):
        self.cleanup()
