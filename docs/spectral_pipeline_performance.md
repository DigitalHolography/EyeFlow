# Spectral pipeline performance audit

Audit date: 2026-10-28. The registered pipeline is `low_freq_spectra`.
This audit adds a profiling tool; it does not change production calculations,
output selection, or scientific definitions.

## Measurements and scope

The requested Justine pressure 0 input under `sandbox/Pressure` has no
`/spectrograms` dataset. The current spectral runner skips it in approximately
0.0003 seconds. Pressure 1–3 also lack that dataset. No cohort archive was read
or tested during this audit, and the pressure files were opened read-only.

Following the supplied representative dimensions, the benchmark generated a
fully allocated `(640,64,128,128)` float32 HDF5 `/spectrograms` dataset:
2,684,354,560 bytes, or **2.5 GiB**. It used contiguous, uncompressed storage.
The generated spatial power is reproducible and varies by pixel and frequency,
with a time-dependent overall gain. It is not a physiological simulation.

The benchmark calls the production spatial reducer, gap-analysis functions,
HDF5 writer, QC calculator and figure exporter. Upstream topology/velocity work
is excluded: it supplies a 75%-coverage mask and eight equally spaced synthetic
beat intervals. Therefore these are timings of the spectral stages, not a full
clinical run. Different masks, beat counts and source layouts can change them.
The execution is serial, as in the current runner.

Environment: Windows, Python 3.13.7, NumPy 2.4.4, h5py 3.16.0, 16 logical CPUs.
These are single-run wall times, without cProfile active during the main run.
Detailed cProfile passes run separately. Source generation took 3.25 seconds
and is excluded from the analysis time. The source had just been written, so
the OS cache may be warm. No cold-disk performance claim is made.

| Non-overlapping work | Seconds | Share of measured total |
|---|---:|---:|
| Spatial reduction and memory operations, excluding source reads | 32.63 | 67.3% |
| HDF5 spectral frame reads, included in the spatial stage below | 1.70 | 3.5% |
| Figure construction and PNG/EPS export | 12.30 | 25.4% |
| HDF5 result writes and final flush | 1.00 | 2.1% |
| Gap analysis, including metrics, Fourier filter and SVD | 0.106 | 0.2% |
| QC calculations | 0.104 | 0.2% |
| Remaining setup/import/dispatch work inside the timed interval | 0.63 | 1.3% |
| **Total** | **48.47** | **100%** |

The complete spatial preprocessing stage took **34.33 seconds**, including
the 1.70 seconds of source access. Do not add this inclusive stage time to
its two components above. HDF5 source timings include h5py allocation/copy
overhead; they are not measurements of physical device service time.

## Main bottleneck: spatial medians and array copies

[`spectrograms.py`](../src/pipelines/low_freq_spectra/spectrograms.py) reads
each HD frame once, then performs all three spatial reductions. For every frame:

1. `np.asarray(align_spectrogram_frame(...), dtype=np.float64)` converts the
   complete 4 MiB float32 frame into an 8 MiB float64 array. The alignment itself
   is a flip/transpose view, but the dtype conversion allocates storage.
2. `frame.reshape(freq,-1)` needs another spatial-order copy for this layout.
3. `full[valid]` creates separate large temporary arrays for the full-frame
   mean and full-frame median, even when every frequency is valid.
4. `np.median` uses its own working copy by default.
5. `frame[:,mask]` gathers gap pixels into another array, followed by
   `samples[gap_valid]` and the median's working storage.

A separate four-frame breakdown at the representative dimensions found the
largest individual operations to be the full-frame median (11.9 ms/frame),
gap median (9.3 ms/frame), and gap pixel gather (7.9 ms/frame). Other copies,
reshape, dtype conversion and finite masks also contribute. These are a small
separate sample, not an exact allocation of the 640-frame total.

The global mean and median are always calculated, including when gap metrics
are the main interest. A source search found no downstream scientific consumer
of `Global/SpatialAverage` or `Global/SpatialMedian`; they are persisted outputs.
Making those outputs optional would therefore merit consideration, but dropping
them would be an explicit product change, not a transparent optimization.

Optimization priority: reduce repeated copies, improve spatial/mask memory
access and reuse median work buffers. Check float64 mean and even-length median
behavior, overflow, included NaN/Inf handling, negative-value clamping, and
the canonical spatial transform against current numerical references.
Do not substitute a mean or approximate quantile for the spatial median.

A small comparison of gathering from the already-created flattened frame
selected exactly the same values, but improved the gather only from about
7.51 to 7.24 ms at 128x128. That isolated indexing substitution is insufficient
to resolve the overall bottleneck.

## Second bottleneck: automatic figure export

[`artifacts.py`](../src/pipelines/low_freq_spectra/artifacts.py) exports
**30 PNG/EPS pairs plus one QC PNG: 61 files** in this workload. Everything is
serial. Saving the 31 PNG files occupied 4.51 seconds, and the 30 EPS files
3.29 seconds; construction, layout and other export work account for the rest.
The largest single pairs were the color gap spectrogram (1.13 s), monochrome
spectrogram (0.77 s), and the three normalized spectral overlays (0.56–0.61 s).

[`writers/article.py`](../src/input_output/writers/article.py) explicitly draws
each figure, freezes its layout, then saves PNG and EPS. Those saves can render
the same artists again. The paired output policy is useful for publication,
but means figure cost is paid even in analysis runs that only need HDF5.

The recently added shape metrics account for **11 additional pairs**, taking
3.11 seconds here. Their numerical packers took only 0.043 seconds combined
(0.005 s vertical and 0.038 s distribution). Their plotting cost is materially
larger than their calculation cost.

Optimization priority: make figure generation and output formats selectable,
or defer export until after numerical analysis. Any render/layout changes
should preserve the shared article styles and PNG/EPS output contracts.

## HD reading and smaller costs

The main preprocessing pass made **640 source reads** and returned exactly
2.5 GiB, confirming one complete frame per time index. Mean, median and gap
reductions share that read. SVD, metrics, QC and exports operate on reduced
arrays, not on the original four-dimensional HD video. The legacy helper
`preprocess_gap_spectrograms` is not separately called by the current runner.

For this storage/cache condition, even eliminating source-read time entirely
would save at most about 3.5% of the measured total. Reading is therefore not
the leading optimization target in this benchmark. However, an actual HD
export could differ: compression, chunks spanning multiple time frames, an
insufficient chunk cache or slower storage can increase source access time.
Inspect and measure the real dataset before tuning batch sizes or chunk cache.
The profiler records chunk shape, compression, storage bytes and cache settings.

Other observed costs:

- The SVD took 0.018 seconds. It sees a `(1024,64)` reduced beat/phase matrix
  in this benchmark, rather than the 2.5 GiB spatial video.
- Shannon entropy, quantiles and Wasserstein packing took 0.038 seconds.
  Repeated frequency sorting in per-phase quantiles is an optimization
  opportunity, but small in absolute terms.
- Result HDF5 writes took roughly 1 second. Many small datasets and attributes
  contribute; raw input reads and numerical calculations should be distinguished
  from this metadata serialization cost.
- QC took 0.104 seconds for eight beats. Leave-one-beat-out SEI sensitivity
  repeatedly reduces almost the same beat array and also computes unused
  flooded metrics. Its work grows with beat count and deserves reevaluation
  for long recordings, but it is not a leading cost here.

## Reproduce the audit

From the repository root in an editable EyeFlow environment:

```powershell
.\.venv\Scripts\python.exe benchmarks/profile_low_freq_spectra.py `
  --generate-synthetic-hd 640 64 128 128 `
  --output .tmp/spectral_profile_new

.\.venv\Scripts\python.exe benchmarks/profile_low_freq_spectra.py `
  --hd 'C:\path\to\scan_HD_output.h5' `
  --output .tmp/spectral_profile_real
```

Use a fresh output directory for source generation; the generated input is
created exclusively and will not overwrite an existing HDF5 file. The first
command needs at least 2.5 GiB plus the reserved free space. The tool writes
`summary.json`, `timings.csv`, isolated result HDF5/figures, and separate
`gap_analysis`/`qc` cProfile binary and text reports. Nested detailed timers are
inclusive and must not be added to their parent stages.

`--frames N` limits a real-file pass to its first N frames. The profiler still
uses a synthetic mask and synthetic beat boundaries, so use it to isolate HD
read/reduction costs rather than to reproduce that patient's numerical metrics.
`--no-figures` measures the same calculations and HDF5 persistence without
image export. `--synthetic-shape` benchmarks generated in-memory frames only.

The measured full-size run and generated source are under
`.tmp/spectral_performance/representative_640/`. The actual pressure-0 skip
inspection is under `.tmp/spectral_performance/justine0/`. No production
optimization has been applied by this audit.
