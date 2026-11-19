"""Optional numeric acceleration backends used by large batched operations."""

from __future__ import annotations

import os
from dataclasses import dataclass
from functools import lru_cache


@dataclass(frozen=True)
class CuPyBackend:
    cupy: object
    ndimage: object


def optional_cupy_backend() -> CuPyBackend | None:
    """Return a usable CuPy backend unless CPU execution was requested.

    ``EYEFLOW_COMPUTE_BACKEND`` accepts ``auto`` (the default), ``cpu``, or
    ``cupy``. Auto mode falls back silently when CuPy or a CUDA device is not
    available; explicit CuPy mode raises a clear runtime error instead.
    """
    preference = os.environ.get("EYEFLOW_COMPUTE_BACKEND", "auto").strip().lower()
    if preference not in {"auto", "cpu", "cupy"}:
        raise ValueError(
            "EYEFLOW_COMPUTE_BACKEND must be 'auto', 'cpu', or 'cupy'."
        )
    if preference == "cpu":
        return None
    backend, error = _probe_cupy_backend()
    if backend is None and preference == "cupy":
        raise RuntimeError(
            "CuPy acceleration was requested but no usable CuPy/CUDA backend "
            "is available."
        ) from error
    return backend


@lru_cache(maxsize=1)
def _probe_cupy_backend() -> tuple[CuPyBackend | None, Exception | None]:
    """Check one real image kernel and FFT, once in each probing process."""
    try:
        import cupy
        from cupyx.scipy import ndimage

        if int(cupy.cuda.runtime.getDeviceCount()) < 1:
            raise RuntimeError("no CUDA device is available")
        sample = cupy.asarray([[1.0, 2.0], [3.0, 4.0]], dtype=cupy.float32)
        ndimage.affine_transform(
            sample,
            cupy.eye(2),
            output_shape=(2, 2),
            order=1,
            mode="grid-constant",
            cval=0.0,
            prefilter=False,
        )
        ndimage.zoom(
            sample,
            (2.0, 2.0),
            order=1,
            mode="grid-constant",
            cval=0.0,
            prefilter=False,
            grid_mode=True,
        )
        cupy.fft.fft(sample, axis=0)
        cupy.cuda.get_current_stream().synchronize()
    except Exception as exc:
        return None, exc
    return CuPyBackend(cupy=cupy, ndimage=ndimage), None
