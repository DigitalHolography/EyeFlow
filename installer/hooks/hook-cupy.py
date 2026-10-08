"""Include Cython imports and JIT resources without collecting development APIs."""
from PyInstaller.utils.hooks import collect_data_files, collect_submodules, copy_metadata

# CuPy constructs a LazyLoader for cupy.testing during every import. Its small
# helper package must exist even though EyeFlow never runs its tests.
hiddenimports = collect_submodules("cupy")
datas = collect_data_files("cupy", includes=["_core/include/**", "cuda/*.cu", "cuda/*.h", "random/*.cu", "random/*.cuh", ".data/*.json"])
datas += copy_metadata("cupy-cuda13x", recursive=True)
