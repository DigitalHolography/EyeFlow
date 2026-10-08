"""Make the single wheel-layout CUDA payload visible to Windows' native loader."""
import os
import sys
from pathlib import Path

if sys.platform == "win32":
    import ctypes

    # Keep handles alive: closing one removes its DLL search directory.
    _eyeflow_dll_directories = [
        os.add_dll_directory(str(Path(sys._MEIPASS) / directory))
        for directory in ("nvidia/cu13/bin/x86_64", "numpy.libs", "scipy.libs")
        if (Path(sys._MEIPASS) / directory).is_dir()
    ]
    # NVRTC loads its builtins internally. AddDllDirectory alone does not
    # cover native calls using the standard LoadLibrary search order.
    _eyeflow_cuda_bin = Path(sys._MEIPASS) / "nvidia/cu13/bin/x86_64"
    os.environ["PATH"] = str(_eyeflow_cuda_bin) + os.pathsep + os.environ.get("PATH", "")
    # Pin the compiler and its matching builtins to this bundle, before any
    # worker compiles a kernel or another CUDA installation can be selected.
    _eyeflow_nvrtc_libraries = [
        ctypes.WinDLL(str(_eyeflow_cuda_bin / name))
        for name in ("nvrtc-builtins64_130.dll", "nvrtc64_130_0.dll")
    ]
