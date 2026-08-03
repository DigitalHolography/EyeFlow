"""Temporary HDF5 workspace for large waveform intermediates."""

from __future__ import annotations

import os
import tempfile
from contextlib import contextmanager
from pathlib import Path

import h5py


@contextmanager
def waveform_scratch_h5(ctx):
    """Yield a scratch H5 and always remove its exact temporary path."""

    output_filename = getattr(ctx.runtime.work_h5, "filename", None)
    preferred_dir = None
    if output_filename:
        candidate = Path(output_filename).resolve().parent
        if candidate.is_dir():
            preferred_dir = candidate
    descriptor, filename = tempfile.mkstemp(
        prefix=".eyeflow-waveform-",
        suffix=".scratch.h5",
        dir=preferred_dir,
    )
    os.close(descriptor)
    path = Path(filename)
    try:
        with h5py.File(path, "w") as scratch:
            scratch.attrs["temporary"] = True
            scratch.attrs["purpose"] = "EyeFlow waveform intermediates"
            yield scratch
    finally:
        path.unlink(missing_ok=True)
