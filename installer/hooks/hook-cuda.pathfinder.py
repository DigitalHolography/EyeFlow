from PyInstaller.utils.hooks import collect_submodules, copy_metadata

hiddenimports = collect_submodules("cuda.pathfinder")
datas = copy_metadata("cuda-pathfinder")
