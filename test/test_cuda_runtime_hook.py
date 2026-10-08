"""Exercise native DLL lookup without a developer CUDA search path."""
import importlib.metadata
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest


@pytest.mark.skipif(sys.platform != "win32", reason="Windows DLL search behavior")
def test_native_nvrtc_lookup_uses_bundled_builtins_with_clean_path(tmp_path):
    try:
        distribution = importlib.metadata.distribution("nvidia-cuda-nvrtc")
    except importlib.metadata.PackageNotFoundError:
        pytest.skip("Install the gpu extra to exercise real NVRTC DLL loading")
    internal = tmp_path / "_internal"
    cuda_bin = internal / "nvidia/cu13/bin/x86_64"
    cuda_bin.mkdir(parents=True)
    for name in ("nvrtc64_130_0.dll", "nvrtc-builtins64_130.dll"):
        source = Path(distribution.locate_file(f"nvidia/cu13/bin/x86_64/{name}"))
        shutil.copyfile(source, cuda_bin / name)
    hook = Path(__file__).resolve().parents[1] / "installer/hooks/runtime_cuda.py"
    script = """
import ctypes, json, runpy, sys
sys._MEIPASS = sys.argv[1]
hook_globals = runpy.run_path(sys.argv[2])
# winmode=0 uses standard native loading rather than Python's user-dir flags.
library = ctypes.WinDLL('nvrtc-builtins64_130.dll', winmode=0)
kernel32 = ctypes.WinDLL('kernel32', use_last_error=True)
get_filename = kernel32.GetModuleFileNameW
get_filename.argtypes = [ctypes.c_void_p, ctypes.c_wchar_p, ctypes.c_uint]
get_filename.restype = ctypes.c_uint
buffer = ctypes.create_unicode_buffer(32768)
assert get_filename(library._handle, buffer, len(buffer))
print(json.dumps(buffer.value))
"""
    environment = dict(os.environ)
    windows = Path(os.environ["SYSTEMROOT"])
    environment["PATH"] = os.pathsep.join([str(windows / "System32"), str(windows)])
    environment.pop("CUDA_PATH", None)
    environment.pop("CUDA_HOME", None)
    result = subprocess.run([sys.executable, "-c", script, str(internal), str(hook)],
                            cwd=tmp_path, env=environment, capture_output=True,
                            text=True, timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr
    assert Path(json.loads(result.stdout)).resolve() == (cuda_bin / "nvrtc-builtins64_130.dll").resolve()
