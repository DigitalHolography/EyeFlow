# RTX 4090 cross-section benchmark

Measured locally on September 9, 2026 using an NVIDIA GeForce RTX 4090 (24 GB),
Windows 11/WDDM, driver 610.88, Python 3.12.14, CuPy 14.2.0, NumPy 2.5.3,
SciPy 1.18.1 and CUDA toolkit 13.1. The driver reports CUDA 13.3 support.

Five synchronized wall-clock measurements per case after a warm-up. Timings
include host/device copies and host results. No concurrent test suite ran during
the recorded benchmark. Inputs are synthetic pulsatile profiles with either
shared or varying validity masks. This is one-segment computation, excluding
branch labeling, final profile postprocessing, and HDF5 export.

## Memory-bounded processing (default 512 MiB scratch budget)

This path includes lazy window extraction, result storage and temporal batching.
The 256-frame cases use two passes and temporal batches. Small geometry metadata
and final outputs are outside the scratch budget.

| Frames / native window | Validity | CPU median | RTX 4090 median | Speedup |
|---|---|---:|---:|---:|
| 64 / 33x33 | shared | 258.1 ms | 25.6 ms | 10.1x |
| 256 / 65x65 | shared | 1209.6 ms | 116.4 ms | 10.4x |
| 256 / 65x65 | varying | 1554.8 ms | 116.1 ms | 13.4x |

## Direct segment measurement (without memory-budget orchestration)

| Frames / validity | Original main CPU | Updated CPU | Updated GPU |
|---|---:|---:|---:|
| 64 / shared | 264.1 ms | 256.1 ms | 21.5 ms |
| 256 / shared | 1043.3 ms | 993.3 ms | 47.0 ms |
| 256 / varying | 1207.0 ms | 1156.2 ms | 47.1 ms |

## Correctness and limitations

Hardware tests exposed unsupported `where=` arguments to `cupy.divide`.
The device path and legacy GPU resize/rotation helpers now use safe denominators
and explicit missing-value assignment. Explicit CuPy mode completed without fallback.

CPU/GPU signals and transverse/longitudinal profiles passed `rtol=1e-5, atol=1e-5`,
including NaN locations; integration limits agreed in all benchmark cases.
Maximum absolute profile discrepancy was approximately 1.96e-5; maximum masked
signal discrepancy was approximately 1.91e-6. Hardware regression coverage also
checks a two-frame batch budget with a fully missing frame.

The full test suite passed: 187 tests, no skips. Ruff passed on changed source/tests.

These are warmed synthetic timings, not a production-data or whole-pipeline speedup.
Individual GPU timings vary under Windows WDDM. CUDA pool reservations reported in
the JSON are allocator-reserved bytes, not a peak total-VRAM measurement.
The original-main comparison is CPU-only; no claim is made about old/new GPU speedup.

## Reproduce

```sh
git show 064f54705fc9428a43c4a9a80423a5bd146407b5:src/calculations/blood_flow_velocity/cross_section/generate_cross_section_signals.py > original_cross_section.py
python benchmarks/benchmark_cross_section_gpu.py --baseline-file original_cross_section.py --output benchmark.json --repeats 5
```

CuPy was installed only in the workspace virtual environment, following the
[CuPy installation documentation](https://docs.cupy.dev/en/stable/install.html).
No GPU driver or system CUDA installation was modified.
