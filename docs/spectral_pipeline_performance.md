# Spectral pipeline performance audit

Audit date: 2026-10-28. The registered pipeline is `low_freq_spectra`.
The original audit below records the baseline implementation. The optimization
follow-up at the end records the implemented native-coordinate reductions and
bounded reads. Scientific outputs and figure selection retain their definitions;
metadata now identifies native HD mask application.

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

## Baseline bottleneck: spatial medians and array copies

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
inspection is under `.tmp/spectral_performance/justine0/`.

## Optimization follow-up

The production reducer now uses native HD spatial order. It prepares native
gap support once, keeps the large float32 sample buffers at float32, flattens
spatial axes as a view, and partitions owned buffers in place for medians.
Float64 is used for mean accumulation and central-value midpoint arithmetic,
preventing overflow without converting all pixels. Global and gap products
share reads and validity checks; when the gap covers the entire field the
median is evaluated once and reused. An arbitrary subset median cannot replace
the whole-field median, so partial coverage still needs two exact selections.

A reusable 64 MiB input buffer receives sequential HDF5 time slabs via
`read_direct`. The standard shape fits 16 frames per slab. Arrays passed by a
caller are copied into this owned buffer before partitioning; source data is
never changed. Metadata records the native HD mask axes, sample dtype, read
budget and median method. Published segmentation remains in EF coordinates.

Replaying the same 2.5 GiB generated fixture on local storage gave:

| Measurement | Baseline | Optimized |
|---|---:|---:|
| Spatial preprocessing, including reads | 34.33 s | 14.91 s |
| Source-access time | 1.70 s | 1.07 s |
| HDF5 source calls | 640 | 40 |
| Total spectral stages, including figures | 48.47 s | 27.87 s |

This is approximately 2.3x faster preprocessing and a 42% shorter total on
these single local measurements. Figure export was not changed and remains
about 11–12 seconds. All **149 persisted datasets**, including all 127 metrics,
were exactly equal between the baseline and optimized fixture runs.

These figures do not establish NAS performance. A dataset slice invokes HDF5
and can cause remote reads; an open file is not a fully resident dataset.
HDF5/OS/SMB caching can combine or satisfy reads, so h5py call counts are not
one-to-one counts of NAS operations. All 2.5 GiB still need to be scanned for
the global outputs. The slab reader reduces API overhead and exposes larger
sequential reads without relying on the entire dataset being cached. Actual
NAS throughput and latency should be measured on the generated HD export.

Optimized results are under `.tmp/spectral_performance/native_optimized_640/`.
The current profiler meters direct slab reads and reports the new block
operations separately. Reuse the generated source with:

```powershell
.\.venv\Scripts\python.exe benchmarks/profile_low_freq_spectra.py `
  --hd .tmp/spectral_performance/representative_640/synthetic_HD.h5 `
  --synthetic-hd-input --output .tmp/spectral_profile_optimized_new
```

## CPU parallelization

Threaded NumPy reductions give a substantial additional improvement. Each
worker owns disjoint, contiguous `(time,frequency)` rows of the existing slab;
it computes the exact median over each entire spatial vector. Threads never
partition overlapping input regions. Global/gap calculations and the validity
rules within a row are unchanged. All futures complete before the next HDF5
read can reuse the buffer, and worker exceptions propagate after pool cleanup.

NumPy documents that many numerical operations release the GIL, permitting
CPU work in threads; this benchmark used ordinary Python with the GIL enabled.
See [NumPy thread safety](https://numpy.org/doc/stable/reference/thread_safety.html).
HDF5 reads stay sequential: h5py's global lock serializes its API calls, so the
workers only receive arrays already in RAM. See
[h5py threading](https://docs.h5py.org/en/stable/threads.html).

The same `(640,64,128,128)` generated HDF5 and 75%-coverage mask were scanned
twice per worker count, reversing the order on the second pass. The times
below are the median of two runs and include source reading and pool startup:

| CPU workers | Spatial pass | Speedup versus one worker |
|---:|---:|---:|
| 1 | 14.74 s | 1.0x |
| 2 | 8.16 s | 1.8x |
| 4 | 4.92 s | 3.0x |
| 8 | 3.33 s | 4.4x |
| 16 | 2.70 s | 5.5x |

The 16-worker result is from a separate two-pass comparison of 8 and 16
workers; the repeated 8-worker result was 3.35 seconds. Mean and both median
arrays were exactly equal across all worker counts. Reading remained roughly
1.05–1.09 seconds, with 40 slab reads and the same 2.5 GiB transfer volume.
These are local cached-storage measurements, not NAS measurements.

Regular execution now defaults to at most eight available logical CPUs, keeps
blocks below 8 MiB serial, and reuses one pool for the whole pass. Eight workers
capture most of the measured gain; sixteen saved another approximately
0.6 seconds while doubling the worker count. The explicit `workers` keyword
and profiler's `--spatial-workers` option allow further machine-specific tuning.
Median reduction introduces no new dependencies or child processes. The 64 MiB
input budget remains unchanged; workers use row views rather than copies of
the whole block. Their gap scratch arrays cover separate row subsets.

A production-path run with automatic eight-worker selection took **3.38 s**
for preprocessing and **16.50 s** including calculations, writes and the same
61 figure files, versus **27.87 s** for the preceding serial implementation.
All **149 persisted datasets** were exactly equal. Figure export now dominates
this local workload at **11.32 s**. Detailed `spatial_operation_timings` in the
profiler remain separate serial kernel measurements, not per-thread accounting.

Reproduce worker scaling with:

```powershell
.\.venv\Scripts\python.exe benchmarks/benchmark_spectral_median_threads.py `
  --hd .tmp/spectral_performance/representative_640/synthetic_HD.h5 `
  --workers 1 2 4 8 --repeats 2 --output .tmp/spectral_threads_new
```

Use `--workers 8 16` for the higher-count comparison. The profiler can compare
complete stages with `--spatial-workers 1` or `--spatial-workers 8`.
Measurements are saved under `.tmp/spectral_performance/thread_scaling/`,
`.tmp/spectral_performance/thread_scaling_16/`, and
`.tmp/spectral_performance/cpu_parallel_640/`.

## Background figure exports

Regular execution now snapshots the small, reduced spectral outputs and gap
mask, then submits figure export to the shared isolated rendering pool. Subsequent
numerical pipelines continue immediately. The snapshot does not include the
HD video or share open HDF5 handles. All original PNG/EPS outputs and article
formatting are retained.

A process isolates Matplotlib from the main pipeline's plotting operations.
Matplotlib documents that it is not thread-safe; using a background thread
would require serializing its access with every other plotting pipeline.
See [Matplotlib threading guidance](https://matplotlib.org/stable/users/faq.html#work-with-threads).
The worker uses Agg and conservative native numeric thread limits. The frozen
launcher calls `multiprocessing.freeze_support()` before importing launch code,
as required by [PyInstaller multiprocessing support](https://pyinstaller.org/en/stable/common-issues-and-pitfalls.html#multi-processing).
The Windows Python spawn path was tested; an installer rebuild was not run.

Background jobs start after the producing descriptor or graph handoff commits. Rendering
uses a private local staging directory outside the EyeFlow result directory.
Completed files are copied to the output filesystem and published with atomic
per-file replacement under explicit artifact ownership. Partial publication is
removed if a job fails. This keeps unrelated descriptor rollback from deleting
concurrently completed exports and supports output directories on another
drive or NAS.

`pdf_report` declares `wait_for_artifacts=True`. The engine waits and publishes
before artifact consumers and before run completion, so reports and archives
cannot race an unfinished renderer. A failed export retains numerical results,
lets independent calculations finish, skips artifact consumers and reports the
run as incomplete. Root `artifact_status`/`artifact_failures` JSON attributes
record those outcomes separately from numerical pipeline status. Worker log
messages and completion are replayed asynchronously through the main run logger.
Direct callers without the managed runtime retain synchronous export.

Shared writers now apply this behavior across pipelines, with at most two
rendering processes reused across files. Analysis of the next file can overlap
unfinished exports, bounded to two unfinished inputs. Final file completion and
the batch return still wait for exports. The spectral pipeline now declares its
source, calibration, spatial, beat, QC and figure stages in its own package;
handoffs commit the reduced products without copying the HD source. See
[pipeline stage graphs](pipeline_stage_graphs.md). The measurements below are
from the earlier spectral-only implementation, not a timing claim for the new
graph/checkpoint and batch execution paths.

An isolated export-only comparison reused the reduced HDF5 result from the
`(640,64,128,128)` fixture, without rescanning its HD input:

| Measurement | Local wall time |
|---|---:|
| Synchronous export, all 61 files | 13.13 s |
| Snapshot reduced products, 1.81 MiB | 0.087 s |
| Snapshot plus process submission, blocking caller time | 0.101 s |
| Simulated later pipeline work, idle interval | 5.00 s |
| Final wait and publication after that interval | 11.69 s |
| Background submission through publication, including the interval | 16.79 s |

All **31 PNG files were pixel-identical**. All **30 EPS files matched** after
excluding creation date/title metadata. This change overlaps rendering with
later work; it does not reduce the work needed to render each figure. The
separate process adds startup/renderer initialization overhead. With no later
work to overlap, that overhead can make completion slower. The five-second
interval above is a synthetic scheduling example, not an actual downstream CPU
workload or an end-to-end EyeFlow measurement. NAS throughput, simultaneous
CPU contention and frozen-installer startup remain unmeasured.

Reproduce with a new isolated output directory:

```powershell
.\.venv\Scripts\python.exe benchmarks/benchmark_spectral_exports.py `
  --h5 .tmp/spectral_performance/cpu_parallel_640/profile.h5 `
  --output .tmp/spectral_exports_new --overlap-seconds 5
```

The clean comparison is saved under
`.tmp/spectral_performance/background_exports_640_quiet/`.
