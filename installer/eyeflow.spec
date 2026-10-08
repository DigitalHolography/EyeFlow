# -*- mode: python ; coding: utf-8 -*-
import os
import sys
from pathlib import Path

from PyInstaller.utils.hooks import collect_data_files, collect_submodules

repo = Path(SPEC).resolve().parent.parent
sys.path.insert(0, str(repo))
from installer.bundle_policy import (
    DEVELOPMENT_EXCLUDES, GPU_EXCLUDES, cuda_payload, normalize_binaries,
)

app_name = os.environ["EYEFLOW_BUILD_NAME"]
wheel_binaries, cuda_headers = cuda_payload()
datas = [
    (str(repo / name), ".")
    for name in ("EyeFlow_logo.png", "default_settings.json", "pyproject.toml")
]
datas += cuda_headers
datas += collect_data_files("sv_ttk")
# TkinterDnD's other-platform binaries and demos are not Windows x64 inputs.
datas += collect_data_files("tkinterdnd2", includes=["tkdnd/win-x64/**"])

a = Analysis(
    [str(repo / "installer/entry_point.py")],
    pathex=[str(repo / "src"), str(repo / "installer")],
    binaries=wheel_binaries,
    datas=datas,
    hiddenimports=["eye_flow", "cli", "graphlib", "matplotlib.backends.backend_ps", "cupyx.scipy.ndimage"] + collect_submodules("pipelines"),
    hookspath=[str(repo / "installer/hooks")],
    hooksconfig={"matplotlib": {"backends": ["TkAgg", "Agg", "pdf", "ps"]}},
    runtime_hooks=[str(repo / "installer/hooks/runtime_cuda.py")],
    excludes=list(DEVELOPMENT_EXCLUDES + GPU_EXCLUDES),
    noarchive=False,
    optimize=0,
)
a.binaries = normalize_binaries(a.binaries, wheel_binaries)
pyz = PYZ(a.pure)
exe = EXE(
    pyz, a.scripts, [], exclude_binaries=True,
    name=app_name, debug=False, bootloader_ignore_signals=False,
    strip=False, upx=False, console=os.environ.get("EYEFLOW_BUILD_CONSOLE") == "1",
    disable_windowed_traceback=False, icon=[str(repo / "EyeFlow.ico")],
)
coll = COLLECT(exe, a.binaries, a.datas, strip=False, upx=False, name=app_name)
