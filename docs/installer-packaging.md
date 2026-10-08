# Windows installer packaging

Build from the repository root with no options:

```powershell
powershell ./build_installer.ps1
```

The PowerShell entry point uses `uv sync --locked` to prepare the `gpu` and
pinned `installer` dependencies in `build/installer-env`, without development
tools. It then invokes that environment's Python directly; uv is not involved
in freezing, staging or running the installed application. The environment
is reused between builds and stays outside the cleaned staging directory.

`installer/build_installer.py` owns the complete build: clean generated output,
invoke PyInstaller, stage pipeline sources, verify the frozen runtime, generate
Inno Setup configuration from `installer/EyeFlow.iss.in`, compile the installer,
and report its size. No virtual environment is included in the installer.
An explicit `-Python` uses an already prepared interpreter instead of uv;
that interpreter must have the `gpu` and `installer` extras installed.
Inno Setup 6 must be installed or provided through `-InnoSetupCompiler`.

The release workflow refreshes `uv.lock` after setting the release version,
then uses the same option-free PowerShell entry point and `dist/installer`
artifact path. It retains the runtime and size reports as workflow artifacts,
including when a build fails. GPU execution is checked when a CUDA device
is available; hosted runners without one record that check as skipped.

## Collection policy

`installer/eyeflow.spec` defines the application bundle. Its hooks preserve
CuPy's Cython imports (including `_softlink`), `graphlib`, JIT source/header
resources, the small `cupy.testing` helper required by CuPy's lazy loader, and
the `cuda.pathfinder` library/header locator. GUI and CLI entry
points and all built-in pipelines are included. The externally discoverable
pipeline source directory is copied once into the bundle without bytecode
caches, and Inno Setup installs that bundle.

EyeFlow's GPU work uses elementwise kernels, reductions, image transforms,
and FFT. It does not call CuPy's BLAS, random, sparse or solver APIs. Their
backend extensions, including the lazily imported `cupyx._cusolver`, are
excluded. Those extensions would otherwise collect very large native
libraries even when no EyeFlow calculation uses them. This is an EyeFlow
payload, not an installation of every public CuPy API; custom pipelines using
additional GPU APIs must extend this policy and its smoke checks.

The native payload comes only from the locked NVIDIA wheels:

| Runtime | Purpose |
| --- | --- |
| CUDA runtime | Device allocation and runtime calls |
| cuFFT | Velocity FFT and Fourier resampling |
| NVRTC and its builtins | Runtime compilation of CuPy kernels |
| nvJitLink | Runtime linking used by CUDA libraries |

The CUDA headers needed for compilation remain under
`_internal/nvidia/cu13/include`. Static libraries, alternative NVRTC binaries,
unused CUDA libraries, test runners, developer tools and other-platform TkDnD
resources are not selected. OpenCV and SimpleITK remain because EyeFlow uses
them for video and displacement products. Matplotlib retains Tk/Agg, PNG,
PDF and PostScript/EPS support.

`installer/build_installer.py` selects wheel binaries during analysis, then
normalizes CUDA entries to one wheel-provided copy under
`_internal/nvidia/cu13/bin/x86_64`. A runtime hook adds that directory to the
Windows DLL search path and keeps its handle alive for the process lifetime.
The hook also prepends that bundled directory to the process PATH and loads
the matching NVRTC compiler and `nvrtc-builtins64_130.dll` by absolute path.
This covers NVRTC's internal native DLL loading as well as Python imports;
no additional copy of either DLL is installed.
Unexpected nested NVIDIA binaries fail the build rather than being silently
discarded. The build's CUDA paths point at the wheel payload, preventing
accidental collection from the developer's system CUDA Toolkit.
NumPy/SciPy BLAS libraries similarly remain once in their original `.libs`
directory, with Windows DLL search handles retained by the runtime hook.

## Verification and size reports

Before Inno Setup runs, the frozen executable checks GUI and CLI imports,
pipeline discovery, Tk drag-and-drop/theme resources, PNG/PDF/EPS, HDF5,
AVI round-trips and SimpleITK. On a GPU machine it also checks fused and
staged transforms and EyeFlow velocity FFT, prohibits CPU fallback in the
GPU calculations, and compares results with CPU references. It uses a fresh
kernel cache and only Windows system directories on the inherited PATH.
GPU kernels first compile on a worker thread, matching GUI analysis.
A machine without an
NVIDIA driver/device records the GPU check as skipped; import and CPU checks
still must pass.

To explicitly require GPU verification on a target machine:

```powershell
& '.\build\installer\pyinstaller-dist\EyeFlow 1.17.0\EyeFlow 1.17.0.exe' `
    --packaging-check "$PWD\gpu-check.json" --require-gpu
```

Use the built application's actual version in its path. This diagnostic writes
only the requested report and temporary test artifacts; it does not process
or replace an EyeFlow result directory.

Each build writes:

- `build/installer/size-report.json`: a compact summary of bundle bytes,
  groups and file count; installer bytes are added after Inno Setup completes.
- `build/installer/bundle-manifest.json`: every bundled file's size, hash,
  source and PyInstaller collection kind, for investigating size regressions.
- `build/installer/smoke-report.json`: the frozen runtime check results.
- `build/installer/pyinstaller-work/eyeflow/xref-eyeflow.html`: the import
  graph for investigating why a Python dependency was included.
- `build/installer/pyinstaller-work/eyeflow/Analysis-00.toc` and
  `COLLECT-00.toc`: analysis and native/data collection records.

The default build rejects duplicate CUDA runtime filenames, an incomplete
CUDA payload, a bundle above 1 GiB, or an installer above 600 MiB. Budget
failures retain their reports so a deliberate dependency change can be
reviewed before adjusting thresholds. Keep the frozen GPU smoke check when
changing CuPy, CUDA or the packaging exclusions.

`-IncludePipelineExtras` prepares the optional pipeline dependencies and
disables the default size budgets; the collection and runtime checks still run.
