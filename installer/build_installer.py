"""Freeze, validate and package EyeFlow using the current Python environment."""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

if __package__:
    from .bundle_policy import cuda_payload, write_manifest
else:
    from bundle_policy import cuda_payload, write_manifest

REPO = Path(__file__).resolve().parent.parent


def project_metadata(repo: Path):
    text = (repo / "pyproject.toml").read_text(encoding="utf-8-sig")
    project = re.search(r"(?ms)^\[project\]\s*\n(.*?)(?=^\[|\Z)", text)
    if project is None:
        raise ValueError("Missing [project] in pyproject.toml")
    values = {}
    for key in ("name", "version"):
        match = re.search(rf'^\s*{key}\s*=\s*"([^"]+)"', project.group(1), re.MULTILINE)
        if match is None:
            raise ValueError(f"Missing [project].{key}")
        values[key] = match.group(1)
    return values["name"], values["version"]


def clean_directory(repo: Path, target: Path):
    """Resolve the target before removing any generated directory."""
    resolved = target.resolve()
    if resolved == repo.resolve() or not resolved.is_relative_to(repo.resolve()):
        raise ValueError(f"Refusing to clean outside the repository: {resolved}")
    if target.exists():
        shutil.rmtree(target)


def resolve_inno(explicit: str | None = None):
    if explicit:
        candidate = Path(explicit).resolve()
        if not candidate.is_file():
            raise FileNotFoundError(f"Inno Setup compiler not found: {candidate}")
        return candidate
    command = shutil.which("ISCC.exe")
    if command:
        return Path(command)
    for variable in ("ProgramFiles(x86)", "ProgramFiles"):
        directory = os.environ.get(variable)
        if directory:
            candidate = Path(directory) / "Inno Setup 6/ISCC.exe"
            if candidate.is_file():
                return candidate
    raise FileNotFoundError("Install Inno Setup 6 or pass -InnoSetupCompiler to build_installer.ps1")


def write_inno_script(repo: Path, build: Path, bundle: Path, output: Path, name: str, version: str):
    version_dir = re.sub(r'[<>:"/\\|?*]+', "-", version).rstrip(" .")
    if not version_dir:
        raise ValueError("Version cannot be converted to an installer directory name")
    if not version_dir.lower().startswith("v"):
        version_dir = "v" + version_dir
    version_info = ""
    if re.fullmatch(r"\d+(\.\d+){0,3}", version):
        parts = version.split(".")
        version_info = "VersionInfoVersion=" + ".".join(parts + ["0"] * (4 - len(parts)))
    values = {
        "APP_NAME": name, "APP_VERSION": version, "DISPLAY_NAME": f"{name} {version}",
        "VERSION_DIR": version_dir, "BUNDLE_DIR": str(bundle.resolve()),
        "OUTPUT_DIR": str(output.resolve()), "SETUP_NAME": f"{name}-Setup-{version}",
        "LICENSE": str(repo / "LICENSE"), "README": str(repo / "README.md"),
        "NOTICES": str(repo / "THIRD_PARTY_NOTICES"), "ICON": str(repo / "EyeFlow.ico"),
        "VERSION_INFO": version_info,
    }
    template = (repo / "installer/EyeFlow.iss.in").read_text(encoding="utf-8")
    # A single replacement pass avoids treating values as template syntax.
    script = re.sub(r"@([A-Z_]+)@", lambda match: values[match.group(1)].replace('"', '""'), template)
    path = build / "EyeFlow.iss"
    path.write_text(script, encoding="utf-8")
    return path


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inno-setup-compiler")
    parser.add_argument("--console", action="store_true")
    parser.add_argument("--skip-clean", action="store_true")
    parser.add_argument("--include-pipeline-extras", action="store_true")
    args = parser.parse_args(argv)
    if sys.platform != "win32":
        parser.error("The installer must be built on Windows")
    name, version = project_metadata(REPO)
    iscc = resolve_inno(args.inno_setup_compiler)
    binaries, _ = cuda_payload()
    wheel_bin = Path(binaries[0][0]).parent
    build = REPO / "build/installer"
    output = REPO / "dist/installer"
    if not args.skip_clean:
        for directory in (build, output):
            clean_directory(REPO, directory)
    build.mkdir(parents=True, exist_ok=True)
    output.mkdir(parents=True, exist_ok=True)
    display_name = f"{name} {version}"
    environment = dict(os.environ)
    environment["CUDA_PATH"] = str(wheel_bin.parent.parent)
    environment.pop("CUDA_HOME", None)
    environment["PATH"] = os.pathsep.join([
        str(wheel_bin),
        *(entry for entry in environment.get("PATH", "").split(os.pathsep)
          if "NVIDIA GPU Computing Toolkit" not in entry),
    ])
    environment["EYEFLOW_BUILD_NAME"] = display_name
    environment["EYEFLOW_BUILD_CONSOLE"] = "1" if args.console else "0"
    print(f"Building {display_name} installer...", flush=True)
    subprocess.run([
        sys.executable, "-m", "PyInstaller", "--noconfirm", "--clean",
        "--distpath", str(build / "pyinstaller-dist"),
        "--workpath", str(build / "pyinstaller-work"), str(REPO / "installer/eyeflow.spec"),
    ], cwd=REPO, env=environment, check=True)
    bundle = build / "pyinstaller-dist" / display_name
    executable = bundle / f"{display_name}.exe"
    if not executable.is_file():
        raise FileNotFoundError(f"PyInstaller did not produce {executable}")
    shutil.copytree(REPO / "src/pipelines", bundle / "pipelines", dirs_exist_ok=True,
                    ignore=shutil.ignore_patterns("__pycache__", "*.pyc", "*.pyo", "AGENTS.md"))
    # Optional pipeline packages can legitimately exceed the default size budget.
    budget = sys.maxsize if args.include_pipeline_extras else 1024**3
    report_path = build / "size-report.json"
    manifest_path = build / "bundle-manifest.json"
    try:
        write_manifest(bundle, build / "pyinstaller-work/eyeflow/COLLECT-00.toc",
                       manifest_path, max_bundle_bytes=budget)
    finally:
        if manifest_path.is_file():
            report = json.loads(manifest_path.read_text(encoding="utf-8"))
            report["file_count"] = len(report.pop("files"))
            report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"Bundle size: {report['bundle_bytes'] / 1024**2:.1f} MiB", flush=True)
    # Verify the frozen app with no build-only CUDA paths or cached kernels.
    smoke_environment = dict(os.environ)
    smoke_environment.pop("CUDA_PATH", None)
    smoke_environment.pop("CUDA_HOME", None)
    windows = Path(os.environ["SYSTEMROOT"])
    smoke_environment["PATH"] = os.pathsep.join([str(windows / "System32"), str(windows)])
    cache = build / "smoke-kernel-cache"
    clean_directory(REPO, cache)
    smoke_environment["CUPY_CACHE_DIR"] = str(cache)
    smoke_environment["EYEFLOW_PIPELINES_DIR"] = str(bundle / "pipelines")
    smoke = build / "smoke-report.json"
    smoke.unlink(missing_ok=True)
    result = subprocess.run([str(executable), "--packaging-check", str(smoke)],
                            cwd=bundle, env=smoke_environment, timeout=120, check=False)
    if result.returncode or not smoke.is_file():
        raise RuntimeError(f"Frozen application smoke check failed; see {smoke}")
    smoke_report = json.loads(smoke.read_text(encoding="utf-8"))
    if smoke_report.get("result") != "passed" or not smoke_report.get("frozen"):
        raise RuntimeError(f"Frozen application smoke check failed; see {smoke}")
    print(f"Frozen smoke check passed ({smoke_report['gpu']}); report: {smoke}", flush=True)
    script = write_inno_script(REPO, build, bundle, output, name, version)
    subprocess.run([str(iscc), str(script)], cwd=REPO, check=True)
    installer = output / f"{name}-Setup-{version}.exe"
    report["installer_bytes"] = installer.stat().st_size
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    if not args.include_pipeline_extras and report["installer_bytes"] > 600 * 1024**2:
        raise RuntimeError(f"Installer exceeds the 600 MiB size budget; see {report_path}")
    print(f"Installer created: {installer}\nInstaller size: {report['installer_bytes'] / 1024**2:.1f} MiB")


if __name__ == "__main__":
    main()
