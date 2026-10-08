"""Explicit Windows GPU payload and auditable PyInstaller collection policy."""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
from collections import defaultdict
from pathlib import Path


# EyeFlow uses CUDA kernels and FFT, not CuPy's BLAS, random or sparse solvers.
CUDA_DLLS = (
    "cudart64_13.dll",
    "cufft64_12.dll",
    "nvJitLink_130_0.dll",
    "nvrtc64_130_0.dll",
    "nvrtc-builtins64_130.dll",
)
GPU_EXCLUDES = (
    "cupy_backends.cuda.libs.cublas",
    "cupy_backends.cuda.libs.curand",
    "cupy_backends.cuda.libs.cusolver",
    "cupy_backends.cuda.libs.cusparse",
    "cupy_backends.cuda.libs.cutensor",
    "cupyx._cusolver",
    "cupyx.distributed",
    "cupyx.tools",
)
DEVELOPMENT_EXCLUDES = ("ipdb", "IPython", "jedi", "pytest", "_pytest")


def cuda_payload(root: Path | None = None):
    """Use only the locked wheels, never a system CUDA Toolkit."""
    if root is None:
        root = Path(importlib.metadata.distribution("nvidia-cuda-runtime").locate_file("nvidia/cu13"))
    bin_dir = root / "bin" / "x86_64"
    binaries = []
    for name in CUDA_DLLS:
        path = bin_dir / name
        if not path.is_file():
            raise FileNotFoundError(f"Required CUDA wheel runtime is missing: {path}")
        binaries.append((str(path), "nvidia/cu13/bin/x86_64"))
    include_dir = root / "include"
    if not (include_dir / "cuda_runtime.h").is_file():
        raise FileNotFoundError(f"CUDA kernel compilation headers are missing: {include_dir}")
    datas = [
        (str(path), str(Path("nvidia/cu13/include") / path.relative_to(include_dir).parent))
        for path in sorted(include_dir.rglob("*")) if path.is_file()
    ]
    return binaries, datas


def normalize_binaries(entries, wheel_binaries):
    """Keep one wheel-provided copy of each CUDA DLL before COLLECT."""
    wheels = {Path(source).name.lower(): (source, destination) for source, destination in wheel_binaries}
    result = []
    seen = set()
    for destination, source, kind in entries:
        name = Path(destination.replace("\\", "/")).name.lower()
        if name in wheels:
            source, directory = wheels[name]
            destination = (Path(directory) / Path(source).name).as_posix()
        elif "nvidia" in Path(destination.replace("\\", "/")).parts or name.startswith(("cublas", "curand", "cusolver", "cusparse")):
            raise ValueError(f"Unexpected NVIDIA binary collected: {destination}")
        elif Path(source).parent.name in {"numpy.libs", "scipy.libs"}:
            destination = (Path(Path(source).parent.name) / Path(source).name).as_posix()
        if destination.lower() in seen:
            continue
        seen.add(destination.lower())
        result.append((destination, source, kind))
    return result


def write_manifest(bundle: Path, toc_path: Path, output: Path, max_bundle_bytes: int = 1024**3):
    """Report sources and hashes, and reject CUDA duplication or size regressions."""
    import ast

    toc = ast.literal_eval(toc_path.read_text(encoding="utf-8"))[0]
    sources = {destination.replace("\\", "/"): (source, kind) for destination, source, kind in toc}
    files = []
    groups = defaultdict(int)
    cuda_names = {name.lower() for name in CUDA_DLLS}
    cuda_seen = set()
    for path in sorted(bundle.rglob("*")):
        if not path.is_file():
            continue
        relative = path.relative_to(bundle).as_posix()
        key = relative.removeprefix("_internal/")
        source, kind = sources.get(key, (None, "PIPELINE_SOURCE" if relative.startswith("pipelines/") else "UNKNOWN"))
        size = path.stat().st_size
        groups[key.split("/")[0]] += size
        name = path.name.lower()
        if name in cuda_names:
            if name in cuda_seen:
                raise ValueError(f"Duplicate CUDA runtime in bundle: {name}")
            cuda_seen.add(name)
        with path.open("rb") as stream:
            digest = hashlib.file_digest(stream, "sha256").hexdigest() if hasattr(hashlib, "file_digest") else _sha256(stream)
        files.append(dict(path=relative, bytes=size, sha256=digest, source=source, kind=kind))
    total = sum(item["bytes"] for item in files)
    report = dict(bundle_bytes=total, groups=dict(sorted(groups.items(), key=lambda item: -item[1])), files=files)
    output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    if cuda_seen != cuda_names:
        raise ValueError(f"CUDA runtime payload incomplete: {sorted(cuda_names - cuda_seen)}")
    if total > max_bundle_bytes:
        raise ValueError(f"Bundle exceeds {max_bundle_bytes / 1024**2:.0f} MiB budget: {total / 1024**2:.1f} MiB; see {output}")
    return report


def _sha256(stream):
    digest = hashlib.sha256()
    for block in iter(lambda: stream.read(1024 * 1024), b""):
        digest.update(block)
    return digest.hexdigest()
