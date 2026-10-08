"""Exercise the packaged runtime without processing or replacing user data."""
from __future__ import annotations

import argparse
import json
import os
import sys
import tempfile
import traceback
from pathlib import Path


def gpu_check(report, require_gpu=False):
    import cupy
    from cupyx.scipy import ndimage

    try:
        count = cupy.cuda.runtime.getDeviceCount()
    except cupy.cuda.runtime.CUDARuntimeError as exc:
        if not require_gpu and exc.status in (35, 100):
            report["gpu"] = f"skipped: NVIDIA driver/device unavailable ({exc.status})"
            return
        raise
    if not count:
        if require_gpu:
            raise RuntimeError("GPU execution required but no CUDA device found")
        report["gpu"] = "skipped: no CUDA device"
        return
    import numpy as np
    from calculations.compute_backend import optional_cupy_backend
    from calculations.topology import transforms
    from pipelines.velocity_analysis.analysis.segments import _VelocityProfileFftAccumulator

    backend = optional_cupy_backend()
    if backend is None:
        raise RuntimeError("EyeFlow selected CPU despite an available CUDA device")
    values = np.random.default_rng(91).normal(size=(8, 29, 29)).astype(np.float32)
    values[:, :4, :] = np.nan

    def reject_cpu(*args, **kwargs):
        raise AssertionError("unexpected CPU fallback")

    affine = transforms.ndi.affine_transform
    transforms.ndi.affine_transform = reject_cpu
    try:
        gpu = transforms.resample_rotate_segment(values, 31.7)
    finally:
        transforms.ndi.affine_transform = affine
    select_backend = transforms.optional_cupy_backend
    transforms.optional_cupy_backend = lambda: None
    try:
        cpu = transforms.resample_rotate_segment(values, 31.7)
    finally:
        transforms.optional_cupy_backend = select_backend
    np.testing.assert_allclose(gpu, cpu, rtol=5e-5, atol=2e-5, equal_nan=True)
    # Exercise the staged ndimage kernels as well as the fused transform.
    device = cupy.asarray(values)
    ndimage.zoom(device, (1, 1.2, 1.2), order=1).get()
    ndimage.rotate(device, 17, axes=(1, 2), order=1, reshape=False).get()
    stack = np.random.default_rng(3402).normal(size=(7, 51, 5)).astype(np.float32)
    stack[:, :4, 0] = np.nan
    mask = np.zeros((51, 5), dtype=bool)
    mask[20:31, 1:4] = True
    options = dict(frame_count=7, ring_count=1, branch_count=1, canvas_side=5,
                   cycle_boundary_indexes=np.asarray([0, 3, 6], dtype=np.int32), index_base=0)
    expected = _VelocityProfileFftAccumulator(**options)
    expected.observe(0, 0, stack, mask)
    actual = _VelocityProfileFftAccumulator(**options)
    actual._write_beat_cpu = reject_cpu
    actual.observe(0, 0, cupy.asarray(stack), mask)
    np.testing.assert_allclose(actual.unmasked, expected.unmasked, rtol=1e-6, atol=1e-6, equal_nan=True)
    np.testing.assert_allclose(actual.masked, expected.masked, rtol=1e-6, atol=1e-6, equal_nan=True)
    cupy.cuda.get_current_stream().synchronize()
    report["gpu"] = "passed: fused/staged transforms, FFT, CPU parity; CPU fallback prohibited"
    report["cupy_version"] = cupy.__version__


def check(report, require_gpu=False):
    import numpy as np
    import cli
    import eye_flow
    import h5py
    import cv2
    import SimpleITK as sitk
    from pipelines import load_pipeline_catalog
    from utils.logger import Logger

    Logger.configure()
    available, missing = load_pipeline_catalog()
    if missing:
        raise RuntimeError(f"Unavailable packaged pipelines: {[(p.name, p.error_msg) for p in missing]}")
    report["pipelines"] = [p.name for p in available]
    report["entry_points"] = "passed: GUI and CLI imports"
    np.testing.assert_allclose(np.linalg.solve(np.diag([2., 4.]), np.asarray([2., 8.])), [1., 2.])
    # Instantiate the actual DnD Tcl extension and theme, then close immediately.
    from tkinterdnd2 import TkinterDnD
    import sv_ttk
    window = TkinterDnD.Tk()
    window.withdraw()
    try:
        sv_ttk.set_theme("dark")
    finally:
        window.destroy()
    report["ui"] = "passed: Tk, drag-and-drop extension and theme"
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    with tempfile.TemporaryDirectory(prefix="eyeflow-packaging-") as folder:
        root = Path(folder)
        figure, axis = plt.subplots()
        axis.plot([0, 1, 2], [1, 0, 2])
        try:
            for suffix in ("png", "pdf", "eps"):
                path = root / f"plot.{suffix}"
                figure.savefig(path)
                if not path.stat().st_size:
                    raise RuntimeError(f"Empty {suffix} export")
        finally:
            plt.close(figure)
        with h5py.File(root / "sample.h5", "w") as h5:
            h5["sample"] = np.arange(8, dtype=np.float32)
        with h5py.File(root / "sample.h5") as h5:
            np.testing.assert_array_equal(h5["sample"][:], np.arange(8))
        path = root / "sample.avi"
        writer = cv2.VideoWriter(str(path), cv2.VideoWriter_fourcc(*"MJPG"), 5, (32, 32))
        if not writer.isOpened():
            raise RuntimeError("AVI writer unavailable")
        writer.write(np.zeros((32, 32, 3), dtype=np.uint8))
        writer.release()
        reader = cv2.VideoCapture(str(path))
        try:
            if not reader.read()[0]:
                raise RuntimeError("AVI round-trip failed")
        finally:
            reader.release()
        sitk.GetArrayFromImage(sitk.GetImageFromArray(np.zeros((8, 8), dtype=np.float32)))
    report["exports"] = "passed: PNG, PDF, EPS, HDF5, AVI and SimpleITK"
    # The GUI runs analysis off the main thread. Compile the first kernels
    # there too, with the builder providing a fresh cache and clean PATH.
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=1) as pool:
        pool.submit(gpu_check, report, require_gpu).result()


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("report", type=Path)
    parser.add_argument("--require-gpu", action="store_true")
    args = parser.parse_args(argv)
    report = dict(executable=sys.executable, frozen=bool(getattr(sys, "frozen", False)))
    os.environ["EYEFLOW_COMPUTE_BACKEND"] = "auto"
    try:
        check(report, args.require_gpu)
        report["result"] = "passed"
        code = 0
    except Exception:
        report["result"] = "failed"
        report["traceback"] = traceback.format_exc()
        code = 1
    args.report.write_text(json.dumps(report, indent=2), encoding="utf-8")
    return code


if __name__ == "__main__":
    raise SystemExit(main())
