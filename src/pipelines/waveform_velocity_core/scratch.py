"""Temporary HDF5 workspace for large waveform intermediates."""

from __future__ import annotations

import os
import tempfile
from contextlib import contextmanager
from pathlib import Path

import h5py

SCRATCH_CHUNK_CACHE_BYTES = 128 * 1024 * 1024


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
        with h5py.File(
            path,
            "w",
            rdcc_nbytes=SCRATCH_CHUNK_CACHE_BYTES,
            rdcc_nslots=1_000_003,
            rdcc_w0=0.75,
        ) as scratch:
            scratch.attrs["temporary"] = True
            scratch.attrs["purpose"] = "EyeFlow waveform intermediates"
            yield scratch
    finally:
        path.unlink(missing_ok=True)
