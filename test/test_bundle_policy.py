"""Native-library source selection and bundle regression guards."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from installer.bundle_policy import CUDA_DLLS, cuda_payload, normalize_binaries, write_manifest


def test_system_cuda_and_nested_copy_resolve_to_one_wheel_library():
    selected = [("/wheels/cufft64_12.dll", "nvidia/cu13/bin/x86_64")]
    entries = [
        ("cufft64_12.dll", "/system/cufft64_12.dll", "BINARY"),
        ("nvidia/cu13/bin/x86_64/cufft64_12.dll", "/wheels/cufft64_12.dll", "BINARY"),
        ("python310.dll", "/python/python310.dll", "BINARY"),
    ]
    assert normalize_binaries(entries, selected) == [
        ("nvidia/cu13/bin/x86_64/cufft64_12.dll", "/wheels/cufft64_12.dll", "BINARY"),
        ("python310.dll", "/python/python310.dll", "BINARY"),
    ]


def test_unexpected_nested_nvidia_binary_is_rejected():
    with pytest.raises(ValueError, match="Unexpected NVIDIA"):
        normalize_binaries([("nvidia/cu13/bin/unused.dll", "/wheel/unused.dll", "BINARY")], [])


def test_unused_cuda_binary_cannot_reenter_through_native_dependencies():
    with pytest.raises(ValueError, match="Unexpected NVIDIA"):
        normalize_binaries([("cusolver64_12.dll", "/system/cusolver64_12.dll", "BINARY")], [])


def test_blas_dependency_is_kept_once_in_its_package_directory():
    entries = [
        ("blas.dll", "/wheel/numpy.libs/blas.dll", "BINARY"),
        ("numpy.libs/blas.dll", "/wheel/numpy.libs/blas.dll", "BINARY"),
    ]
    assert normalize_binaries(entries, []) == [("numpy.libs/blas.dll", "/wheel/numpy.libs/blas.dll", "BINARY")]


def make_payload(tmp_path):
    root = tmp_path / "wheel"
    binaries = root / "bin/x86_64"
    binaries.mkdir(parents=True)
    for name in CUDA_DLLS:
        (binaries / name).write_bytes(b"runtime")
    headers = root / "include"
    headers.mkdir()
    (headers / "cuda_runtime.h").write_text("header")
    (root / "lib").mkdir()
    (root / "lib/unused.lib").write_bytes(b"static library")
    (binaries / "unused.alt.dll").write_bytes(b"unused alternative")
    return root


def test_payload_omits_static_libraries_and_alternative_binaries(tmp_path):
    binaries, datas = cuda_payload(make_payload(tmp_path))
    assert {Path(source).name for source, _ in binaries} == set(CUDA_DLLS)
    assert {Path(source).name for source, _ in datas} == {"cuda_runtime.h"}


def test_missing_required_wheel_runtime_fails_before_freezing(tmp_path):
    root = make_payload(tmp_path)
    (root / "bin/x86_64" / CUDA_DLLS[0]).unlink()
    with pytest.raises(FileNotFoundError, match="Required CUDA wheel runtime"):
        cuda_payload(root)


def make_bundle(tmp_path):
    bundle = tmp_path / "bundle"
    internal = bundle / "_internal"
    internal.mkdir(parents=True)
    entries = []
    for name in CUDA_DLLS:
        path = internal / name
        path.write_bytes(b"runtime")
        entries.append((name, f"/wheel/{name}", "BINARY"))
    toc = tmp_path / "COLLECT-00.toc"
    toc.write_text(repr((entries,)))
    return bundle, toc, tmp_path / "report.json"


def test_manifest_reports_hashes_and_binary_provenance(tmp_path):
    bundle, toc, output = make_bundle(tmp_path)
    report = write_manifest(bundle, toc, output)
    assert report["bundle_bytes"] == len(CUDA_DLLS) * len(b"runtime")
    assert all(len(item["sha256"]) == 64 for item in report["files"])
    assert all(item["source"].startswith("/wheel/") for item in report["files"])
    assert json.loads(output.read_text()) == report


def test_manifest_rejects_duplicate_cuda_runtime(tmp_path):
    bundle, toc, output = make_bundle(tmp_path)
    other = bundle / "_internal/nvidia/cu13/bin/x86_64"
    other.mkdir(parents=True)
    (other / CUDA_DLLS[0]).write_bytes(b"runtime")
    with pytest.raises(ValueError, match="Duplicate CUDA"):
        write_manifest(bundle, toc, output)


def test_size_regression_writes_report_and_fails(tmp_path):
    bundle, toc, output = make_bundle(tmp_path)
    with pytest.raises(ValueError, match="size|budget"):
        write_manifest(bundle, toc, output, max_bundle_bytes=1)
    assert output.is_file()
