"""Temporal image reductions on CPU or CUDA from bounded sequential reads."""

from contextlib import closing
from math import prod
from typing import Literal

import numpy as np

from calculations.compute_backend import optional_cupy_backend


_BLOCK_BYTES = 32 * 1024 * 1024
_LAST_HISTOGRAM_BIN = 65535
_LOW_TAIL = 0.0001
_HIGH_TAIL = 0.9999
DisplayScale = Literal["raw", "global", "per_frame"]


def frame_blocks(video):
    """Yield native-dtype blocks valid until the next read, aligned when possible."""
    if not video.shape[0]:
        return
    frame_samples = prod(video.shape[1:])
    frames = min(video.shape[0], max(1, _BLOCK_BYTES // (frame_samples * 8)))
    chunks = getattr(video, "chunks", None)
    if chunks and chunks[0] <= frames:
        frames -= frames % chunks[0]

    stream = getattr(video, "iter_blocks", None)
    if callable(stream) and frame_samples * video.dtype.itemsize <= _BLOCK_BYTES:
        with stream(frames=frames) as blocks:
            for lease in blocks:
                with lease as block:
                    yield block
        return

    direct_read = getattr(video, "read_direct", None)
    if callable(direct_read):
        buffer = np.empty((frames, *video.shape[1:]), dtype=video.dtype)
    for start in range(0, video.shape[0], frames):
        stop = min(start + frames, video.shape[0])
        if callable(direct_read):
            block = buffer[:stop-start]
            direct_read(block, source_sel=np.s_[start:stop])
            yield block
        else:
            yield video[start:stop]


def temporal_image_statistics(video, *, display: DisplayScale = "raw") -> dict[str, np.ndarray]:
    """Return frame-axis mean, sample deviation, and source-intensity entropy.

    Display reductions follow HoloDoppler's flash subtraction and contrast
    scaling. Raw reductions omit nonfinite samples. Entropy always uses finite
    source magnitudes normalized across frames at each pixel and channel.
    """
    backend = optional_cupy_backend()
    xp = np if backend is None else backend.cupy
    shape = video.shape[1:]
    count = xp.zeros(shape, dtype=xp.int64)
    mean = xp.zeros(shape, dtype=xp.float64)
    squared_deviations = xp.zeros_like(mean)
    mass = xp.zeros_like(mean)
    mass_log_mass = xp.zeros_like(mean)
    source_seen = xp.zeros(shape, dtype=bool)
    low, high = -np.inf, np.inf
    if display == "global":
        low, high = np.inf, -np.inf

    with closing(frame_blocks(video)) as blocks:
        for host in blocks:
            # A leased HDF5 buffer can be reused as soon as we advance the iterator.
            source = host if backend is None else xp.asarray(host, blocking=True)
            source = source.astype(xp.float64, copy=False)
            finite = xp.isfinite(source)
            values = xp.where(finite, source, 0.0)
            source_seen |= finite.any(axis=0)
            magnitude = xp.abs(values)
            mass += magnitude.sum(axis=0)
            mass_log_mass += (magnitude * xp.log(magnitude + (magnitude == 0))).sum(axis=0)

            if display != "raw":
                spatial_axes = (1, 2)
                frame_count = finite.sum(axis=spatial_axes, keepdims=True)
                frame_mean = values.sum(axis=spatial_axes, keepdims=True) / xp.maximum(frame_count, 1)
                values = values - frame_mean
                if display == "per_frame":
                    values = _stretch_frames(values, finite, xp)
                else:
                    low = xp.minimum(low, xp.where(finite, values, xp.inf).min())
                    high = xp.maximum(high, xp.where(finite, values, -xp.inf).max())

            if display == "per_frame":
                block_count = values.shape[0]
                block_mean = values.mean(axis=0)
                centered = values - block_mean
            else:
                block_count = finite.sum(axis=0)
                block_mean = xp.where(finite, values, 0.0).sum(axis=0) / xp.maximum(block_count, 1)
                centered = xp.where(finite, values - block_mean, 0.0)

            new_count = count + block_count
            fraction = block_count / xp.maximum(new_count, 1)
            delta = block_mean - mean
            squared_deviations += xp.einsum("t...,t...->...", centered, centered)
            # Merge second moments without a difference of squared raw moments.
            squared_deviations += delta * delta * count * fraction
            mean += delta * fraction
            count = new_count

    if display == "global" and video.shape[0]:
        # Apply the common scale once, then merge the invalid samples as zeros.
        has_range = high > low
        span = xp.where(has_range, high - low, 1.0)
        scaled_mean = (mean - xp.where(has_range, low, 0.0)) / span
        fraction = count / video.shape[0]
        squared_deviations = squared_deviations / span**2 + scaled_mean**2 * count * (1 - fraction)
        mean = xp.where(has_range, scaled_mean * fraction, 0.0)
        deviation = xp.where(has_range, xp.sqrt(squared_deviations / max(video.shape[0] - 1, 1)), 0.0)
    else:
        mean = xp.where(count > 0, mean, xp.nan)
        deviation = xp.where(count > 0, xp.sqrt(squared_deviations / xp.maximum(count - 1, 1)), xp.nan)

    # H(w / sum(w)) = log(sum(w)) - sum(w * log(w)) / sum(w).
    safe_mass = xp.where(mass > 0, mass, 1.0)
    entropy = xp.maximum(xp.log(safe_mass) - mass_log_mass / safe_mass, 0.0)
    entropy = xp.where(source_seen, entropy, xp.nan)
    images = {"avg": mean, "std": deviation, "entropy": entropy}
    return images if backend is None else {name: xp.asnumpy(image) for name, image in images.items()}


def _stretch_frames(values, finite, xp):
    flat = values.reshape(values.shape[0], -1)
    valid = finite.reshape(flat.shape)
    low = xp.where(valid, flat, xp.inf).min(axis=1)
    high = xp.where(valid, flat, -xp.inf).max(axis=1)
    active = high > low
    origin = xp.where(active, low, 0.0)
    span = xp.where(active, high - low, 1.0)
    samples = xp.clip((flat - origin[:, None]) / span[:, None], 0, 1)

    # Below 10,000 samples the two crossings are the minimum and maximum.
    if flat.shape[1] >= 10_000:
        low, high = _stretch_limits(samples, valid, xp)
        samples = xp.clip((samples - low[:, None]) / (high - low)[:, None], 0, 1)
    return xp.where(valid & active[:, None], samples, 0.0).reshape(values.shape)


def _stretch_limits(samples, valid, xp):
    """Reproduce the two crossings of the normalized 65,536-bin histogram."""
    if xp is np:
        # Quantization is monotone, so only the two order statistics need bins.
        low = np.zeros(samples.shape[0])
        high = np.ones_like(low)
        for frame, (row, mask) in enumerate(zip(samples, valid)):
            values = row[mask]
            if values.size < 10_000:
                continue
            first = int(np.floor(_LOW_TAIL * values.size))
            last = int(np.ceil(_HIGH_TAIL * values.size)) - 1
            values.partition((first, last))
            low[frame], high[frame] = np.floor(values[[first, last]] * _LAST_HISTOGRAM_BIN + 0.5) / _LAST_HISTOGRAM_BIN
    else:
        # Put each frame into a disjoint range so one CUDA histogram covers the block.
        bins_per_frame = _LAST_HISTOGRAM_BIN + 2
        bins = xp.floor(samples * _LAST_HISTOGRAM_BIN + 0.5).astype(xp.int32)
        bins = xp.where(valid, bins, bins_per_frame - 1)
        offsets = xp.arange(samples.shape[0], dtype=xp.int32) * bins_per_frame
        hist = xp.bincount((bins + offsets[:, None]).ravel(), minlength=samples.shape[0] * bins_per_frame)
        counts = hist.reshape(-1, bins_per_frame)[:, :-1].cumsum(axis=1)
        total = counts[:, -1:]
        low = (counts > _LOW_TAIL * total).argmax(axis=1) / _LAST_HISTOGRAM_BIN
        high = (counts >= _HIGH_TAIL * total).argmax(axis=1) / _LAST_HISTOGRAM_BIN

    separated = high > low
    return xp.where(separated, low, 0.0), xp.where(separated, high, 1.0)
