"""RAM-backed HDF5 workspace for heartbeat intermediates."""

from __future__ import annotations

from contextlib import contextmanager
from uuid import uuid4

import h5py

SCRATCH_CHUNK_CACHE_BYTES = 128 * 1024 * 1024
SCRATCH_BLOCK_BYTES = 64 * 1024 * 1024


@contextmanager
def heartbeat_scratch_h5(_ctx):
    """Yield a non-persistent HDF5 file allocated entirely in RAM."""

    filename = f"eyeflow-heartbeat-{uuid4().hex}.h5"
    with h5py.File(
        filename,
        "w",
        driver="core",
        backing_store=False,
        block_size=SCRATCH_BLOCK_BYTES,
        rdcc_nbytes=SCRATCH_CHUNK_CACHE_BYTES,
        rdcc_nslots=1_000_003,
        rdcc_w0=0.75,
    ) as scratch:
        scratch.attrs["temporary"] = True
        scratch.attrs["storage"] = "memory"
        scratch.attrs["purpose"] = "EyeFlow heartbeat intermediates"
        yield scratch
