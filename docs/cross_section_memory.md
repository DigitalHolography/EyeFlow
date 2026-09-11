# Cross-section memory and numerical behavior

Cross-section radius measurements and annulus masks now share pixel distances
normalized by the image half-diagonal. This corrects mask-based orientation;
it can intentionally change measured angles and resulting signals.

## Sparse payloads

Generated signals, transverse/longitudinal profiles, centered profiles, rotated
movies, masks, and mean images use `SegmentArray`. Their logical shape remains
`(ring, branch, ...)`, but only measured intersections own payload arrays.
Small coordinate, fit, and validity metadata arrays remain rectangular.

Use `array[ring, branch]` to read a segment without materializing other segments.
Missing entries return a read-only broadcast of NaN (False for masks).
`array.segment_indexes` lists stored entries; `array.nbytes` counts stored payload
bytes. Assignments require an integer ring/branch pair. This object is not a full
NumPy ndarray replacement: explicitly use `np.asarray(array)` when an external
consumer requires a dense ndarray. General slicing also materializes a dense array.

The existing rectangular HDF5 export schema is preserved. Segment-map temporal
interpolation consumes sparse movies directly and skips absent segments, but its
returned export array is dense. Profile and mask exporters may also materialize
dense arrays. Export memory is outside the calculation working-memory budget.

## Working-memory budget

```python
settings = CrossSectionSignalSettings(
    hydrodynamic_diameters=True,
    velocity_profile_threshold=0.5,
    rotate_from_mask=True,
    pixel_size_mm=0.0191,
    working_memory_mb=512.0,
)
```

`working_memory_mb` is a conservative estimate-based budget for concurrent
calculation scratch arrays, not a process RSS limit. It excludes the input cube,
cached geometry, retained results, library allocator caches, and export arrays.
The existing `EYEFLOW_MAX_PARALLEL_JOBS` cap also applies. GPU work is serial.

Input windows are extracted inside workers, with at most one submitted job per
worker. Completed workers store directly into separate segment slots rather than
queuing large movies behind a slower earlier result.

When a full temporal window would exceed the budget, processing uses two passes:
one accumulates the global resized mean, then another projects temporal batches
using the same global angle and integration limits. The additional resize pass
trades computation for lower peak scratch memory. A budget too small for even one
frame raises `MemoryError`. Global batch means accumulate in float64; very small
floating-point differences from the full-window float32 reduction are expected.

## GPU execution

The segment GPU path keeps resized, masked, padded, and rotated arrays on the
device through profile reduction. Small mean images cross to the CPU for
orientation and quadratic fitting. Only the unmasked movie is downloaded;
the masked movie is reduced and discarded on-device. Temporal batching also
applies to this path; the global-mean pass uses the shared resize helper.

In automatic backend mode, a device-operation failure emits a warning and
disables further cross-section GPU attempts for the process. Restart the process
to retry automatic GPU execution. Explicit `EYEFLOW_COMPUTE_BACKEND=cupy` raises
on device-operation failure instead of silently switching backends.

## Missing data and fits

- All-missing frames remain NaN rather than becoming zero flow.
- Empty results consistently have zero branches.
- Disconnected masks use the largest component to choose the window center;
  ties retain component-label order. The mask itself retains all intersection pixels.
- Pixel size must be positive and finite; the profile threshold must lie strictly
  between zero and one; memory and percentile settings are validated.
- Quadratic limits are fit in pixel units and reject rank-deficient, flat,
  upward-opening, or nonfinite fits. Failed fits retain the existing full-width
  fallback. The orientation search strategy is unchanged.

## Validation

Run `python -m pytest test -q` with the project dependencies and pytest installed.
The optimization regressions include sparse storage/export compatibility,
rectangular/off-center geometry, temporal batching, missing frames, fit edge cases,
worker caps, lazy extraction, and GPU fallback behavior. A NumPy/SciPy device stand-in
checks the GPU algorithm and download boundaries; a separate CuPy test runs only
when CuPy and a usable CUDA device are available. Local RTX 4090 hardware results are recorded in
[cross_section_gpu_benchmark.md](cross_section_gpu_benchmark.md). Production-data
performance still needs to be measured on representative recordings.
