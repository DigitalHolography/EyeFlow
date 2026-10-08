"""Cython imports cannot be inferred from Python bytecode."""
from PyInstaller.utils.hooks import collect_submodules

unused = ("cupy_backends.cuda.libs.cublas", "cupy_backends.cuda.libs.curand", "cupy_backends.cuda.libs.cusolver", "cupy_backends.cuda.libs.cusparse", "cupy_backends.cuda.libs.cutensor")
hiddenimports = collect_submodules("cupy_backends", filter=lambda name: name not in unused)
