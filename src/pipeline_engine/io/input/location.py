"""Storage classification without probing throughput or reading source data."""

import os
import ntpath
from pathlib import Path
import sys


def _windows_drive_type(root):
    import ctypes
    function = ctypes.windll.kernel32.GetDriveTypeW
    function.argtypes = [ctypes.c_wchar_p]
    function.restype = ctypes.c_uint
    return function(root)


def storage_location(filename):
    """Return local/network/unknown; an unrecognized filesystem stays unknown."""
    name = os.fspath(filename)
    windows = name.replace("/", "\\")
    if windows.lower().startswith("\\\\?\\unc\\"):
        return "network"
    if windows.startswith("\\\\") and not windows.startswith(("\\\\?\\", "\\\\.\\")):
        return "network"
    if sys.platform == "win32":
        try:
            resolved = str(Path(name).expanduser().resolve())
            if resolved.lower().startswith("\\\\?\\unc\\") or (
                    resolved.startswith("\\\\") and not resolved.startswith(("\\\\?\\", "\\\\.\\"))):
                return "network"
            drive, _ = ntpath.splitdrive(resolved)
            kind = _windows_drive_type(drive + "\\")
            return "network" if kind == 4 else "local" if kind in (2, 3, 5, 6) else "unknown"
        except (OSError, ValueError):
            return "unknown"
    if sys.platform.startswith("linux"):
        try:
            resolved = str(Path(name).expanduser().resolve())
            mounts = []
            for line in Path("/proc/self/mountinfo").read_text().splitlines():
                left, right = line.split(" - ", 1)
                mount = left.split()[4].replace("\\040", " ").replace("\\134", "\\")
                if resolved == mount or resolved.startswith(mount.rstrip("/") + "/"):
                    mounts.append((len(mount), right.split()[0]))
            filesystem = max(mounts)[1] if mounts else ""
            if filesystem in {"nfs", "nfs4", "cifs", "smb3", "ceph", "afs", "9p", "fuse.sshfs"}:
                return "network"
            if filesystem in {"ext2", "ext3", "ext4", "xfs", "btrfs", "zfs", "tmpfs", "overlay"}:
                return "local"
        except (OSError, ValueError, IndexError):
            pass
    if sys.platform == "darwin":
        import psutil
        try:
            resolved = Path(name).expanduser().resolve()
            for mount in sorted(psutil.disk_partitions(all=True), key=lambda item: len(item.mountpoint), reverse=True):
                if resolved.is_relative_to(mount.mountpoint):
                    filesystem = mount.fstype.lower()
                    if filesystem in {"nfs", "smbfs", "cifs", "afpfs", "webdav", "sshfs", "fuse.sshfs"}:
                        return "network"
                    return "local" if filesystem in {"apfs", "hfs", "hfs+", "msdos", "exfat"} else "unknown"
        except (OSError, ValueError):
            pass
    return "unknown"
