# Windows installer packaging

Build from the repository root with `powershell ./build_installer.ps1`. The
script prepares a locked installer environment with `uv`, then calls
`installer/build_installer.py` to freeze EyeFlow, stage pipeline sources,
run frozen-runtime checks, and compile the Inno Setup installer. Inno Setup 6
is required; `-InnoSetupCompiler` can point to it. `-Python` selects an
already prepared interpreter instead of syncing the installer environment.
The release workflow uses the same entry point; see
`.github/workflows/release.yml` for tag and publication rules.

## Running the installed CLI

The installed executable opens the GUI without arguments and accepts CLI
arguments in an existing terminal:

```powershell
& 'C:\path\to\EyeFlow <version>.exe' --help
& 'C:\path\to\EyeFlow <version>.exe' --data 'D:\Data\sample.holo' --output 'D:\Results'
$LASTEXITCODE
```

Use the actual installed version and path. The installer does not add a
separate `eyeflow-cli` command to `PATH`. Exit codes are 0 for success, 1 for
runtime failure, and 2 for invalid arguments. Saved settings govern pipeline
selection unless `--pipelines` supplies target names. `-Console` enables a
visible console for diagnostic desktop builds.

## Bundle policy and checks

`installer/eyeflow.spec` and build hooks include GUI/CLI entry points, built-in
pipelines, Tk resources, HDF5/media support, and the CuPy/CUDA parts EyeFlow
uses. The bundle omits unused CuPy APIs and duplicate CUDA libraries to control
size. Custom pipelines that need more GPU APIs must extend the spec and smoke
checks. The exact native library selection is implementation detail owned by
`installer/`, not a data contract.

Before compiling the installer, the builder checks frozen imports, pipeline
discovery, GUI resources, image/PDF/EPS support, HDF5, AVI, and SimpleITK. A
machine with CUDA also checks representative GPU calculations against CPU
references; a machine without it records GPU checks as skipped. CLI subprocess
checks cover help and error exit codes. Reports in `build/installer/` include
`size-report.json`, `bundle-manifest.json`, and `smoke-report.json`.

Default size budgets are 1 GiB for the bundle and 600 MiB for the installer.
`-IncludePipelineExtras` includes optional pipeline dependencies and disables
those budgets. Keep the frozen GPU check when changing CuPy, CUDA, or bundle
exclusions. The build script and installer tests own detailed thresholds and
paths.
