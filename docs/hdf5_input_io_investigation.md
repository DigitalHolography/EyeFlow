# HDF5 input I/O investigation

Investigation target: `refactor/pipeline-execution-ownership`, commit
`084673172710bc0906848060f79291822542d546`. Production defaults are unchanged.
Source files are opened read-only; benchmark results use private local folders.

## Evidence and scope

The representative acquisition is `260616_ABS1850789_R_1`, supplied as matching
SSD and NAS `.holo` paths. Its HD file is 11,009,977,196 bytes; the local DV file
is 2,052,154,188 bytes. The input byte count is **not** the combined file size:
EyeFlow requests only selected datasets. No vessel segmentation is performed.

Early exploratory measurements ran on `dev` (`e7220ba9`) before the branch
correction. They are kept separately under `.tmp/input_io/ssd` and
`.tmp/input_io/nas`; they must not be presented as complete-workflow results
for this refactor. Scientific pipeline code and the spectral reducer are
unchanged between those revisions (apart from displacement scratch imports),
but execution ownership, immutable-array reuse and resource control differ.
Correct-branch runs record their revision in `results.json`.

## Verified current path

`runtime._execute_input` opens HD/DV, initializes the output model, walks source
metadata and creates one `InputService`. Its serving thread handles requests
serially, independently of compute-slot admission. Workers receive `InputFile`
metadata and clients, never open source handles. `InputDataset._value` resolves
the canonical source identity, so hard links share replay and client mappings.

`InputDataset.__getitem__` applies the selection **after** `InputClient.read`
has obtained the entire dataset. `_stage` persists numeric/fixed-string inputs
of at least 4 KiB; object metadata goes through a value reference. For a large
dataset, `ProductStorage._write_array` reads first-axis slabs of about 16 MiB,
writes a local `.npy`, flushes and closes the writable mapping, and only then
publishes a read-only reference. The first small selection therefore encounters
a whole-dataset barrier. Subsequent selections and consumers reuse local data.

The thread is request driven: it does not currently prefetch an upcoming
dataset. Dataset dependencies are absent from pipeline descriptors; `input_slot`
does not describe actual dataset requirements. The compute DAG alone cannot
safely predict the exact input sequence.

`ProcessExecutor` owns workers and transactional acceptance; `scheduler` owns
readiness. Input mappings close after workers finish. `PendingOutput` retains
product storage through exports and final publication; publication is implemented
by `input_output/writers/h5_publication.py`. Do not mix those responsibilities
into an input prefetcher.

## Dataset/consumer audit

Sizes and chunks below are actual HD/DV metadata from R_1. All listed HD
datasets have float32 samples, no compression or other filters. HD cache access
properties observed with h5py 3.16 were 8 MiB / 8191 slots / w0=0.75. These are
environment observations, not hard-coded EyeFlow settings. Dataset-size formulas
are `T*Y*X*itemsize` and `T*F*Y*X*itemsize`, respectively.

| Dataset / input family | Consumers | Access pattern | Decoded size | Chunking | Reuse | Optimization opportunity |
|---|---|---|---|---|---|---|
| HD `moment0` / `M0`, `moment2` / `M2` | `velocity`; M0 also optional topology reference and displacement | Velocity: paired sequential 32-frame blocks, whole time axis eventually; topology orientation can rescan blocks for each unresolved segment | 640 MiB each, `(640,512,512)` | `(32,32,32)`, 128 KiB/chunk | Service deduplicates M0 across consumers; derived velocity is retained for downstream spatial work | Start paired blocks earlier; preserve replay where topology/displacement rescan. Direct-to-mapping reads can avoid temporary slabs without changing science |
| HD lowest/highest `band_<i>_<lo>_<hi>` | Band-ratio `velocity`; lowest band is topology fallback reference | Paired 32-frame sequential blocks; validate finite samples, clamp negatives independently, compute calibrated ratio | 640 MiB each; here `band_0_3000_6000`, `band_1_9000_18000` | `(32,32,32)` | Shared cache across source adapters/variants; moments and bands are distinct inputs | Resolve actual names from numeric bounds, not hard-coded 3000–9000 LF name; same paired-block opportunity |
| HD `moment0ff` / `M0FF` | `topology_core`; optional `spatial_gradient_moment0` for gradient BVR | Topology: spatial tiles over **all** time samples, float64 temporal nanmedian. Gradient: ordered temporal blocks with six-frame periodic halo | 640 MiB | `(32,32,32)` | Local replay is useful for incompatible spatial/temporal access; successful continuous gradient extraction shared between velocity variants | Keep replay provisionally. Frame-only streaming cannot replace temporal medians; direct spatial tiles may reread many chunks |
| HD `spectrograms` | Hidden, run-once `low_freq_spectra` after velocity/topology | Single sequential time pass; global spatial mean/median and gap median in one pass. Default 64 MiB working blocks (16 frames here) | 2560 MiB, `(640,64,128,128)` | `(64,4,16,16)`, 256 KiB/chunk | Source needed once for reduction; subsequent beats, filters, QC and metrics use `(time,freq)` reductions | Strongest bounded streaming candidate. Begin reducing first block without full replay; retain dtype, gap gathering and owned mutable partition buffer |
| DV retinal artery/vein/labeled masks; vessel/choroid/disc masks and retinal logits | Canonical source assembly; topology, velocity and displacement | Whole small 2-D arrays; Boolean conversions/alignment, independent copies where edited | Int64 masks 2 MiB each; Boolean masks 256 KiB; float32 logits 1 MiB | Contiguous, uncompressed in this acquisition | Multiple adapters/workflows; one source staging per identity | Low priority; cache already prevents repeated source transfers; preserve caller isolation |
| HD/DV JSON/scalar metadata, timing, pitch, disc geometry | Source assembly and provenance | Small whole values | Small; variable-length JSON `nbytes` measures descriptors, not text payload | Mostly contiguous | Reused source parameters and client cache | Measure metadata traversal/open latency separately before optimizing |
| HD `registration`, optional Zernike pass-through | Output initialization when `copy_hd_passthrough` is enabled | Selected objects copied through output initialization, outside input-client replay | Registration 5120 bytes here; actual configured pass-through availability varies | Registration `(1024,2)` | Separate output-owned persistence | Keep optional; do not copy entire HD just to obtain metadata |
| Derived velocity, per-beat signals, segment/profile products | `velocity_analysis`, `absolute_waveform_metrics`, `lowrank_waveform_decomposition`, `waveform_shape_metrics`, `blood_volume_rate` | Shared derived state and packed outputs, not repeated HD volume reads | Selection dependent; not source-HDF5 bytes | Local arrays/output storage | Reused within each isolated velocity workflow | Do not attribute downstream compute/exports to NAS input. Gradient BVR adds M0ff through its declared dependency |
| Other HD color/directional videos, M1, averaged spectrum | No consumers in the audited requested families | Unrequested | Color and directional videos 1920 MiB each; M1 640 MiB; averaged spectrum 256 MiB | Export dependent | None in these workloads | Whole-file staging transfers these unused bytes |

Code evidence: `pipelines/velocity/estimation.py` (`_estimator_chunk_slices`,
`_read_volume_chunk`); `pipelines/topology_core/gap_masks.py`
(`temporal_median_image`); `calculations/vessel_segments/sampling/preparation.py`
(`resolve_segment_rotations`); `pipelines/low_freq_spectra/spatial_reduction.py`
(`spectral_batches`, `reduce_spatial_block`); `pipelines/vessel_inputs.py`;
`input_output/schema/holodoppler.py`; the downstream pipelines' `runner.py` files.

### When requirements become known

Run resolution knows targets, option closure and velocity-method selection.
Source metadata reveals root aliases, optional M0ff/spectra and numeric band
names. Canonical assembly resolves source validity and mask alignment; topology
determines whether orientation fallback needs a reference scan. Gradient BVR's
option requires M0ff; displacement's selected root moment is configuration
dependent. Velocity methods can fail independently, so speculative reads must
not turn an unused input's error into a failure of another workflow.

Explicit source requirements could eventually be declared separately from DAG
keys. Initially restrict any lookahead to a validated consumer's known sequence.
Avoid eagerly loading all DAG-associated datasets, especially both velocity
methods when only one was selected.

### Copies and lifetimes

Necessary working copies include dtype conversion when required, temporal
nanmedian's finite-cleaned float64 tiles, in-place spatial-median partition
buffers, gathered gap samples, independently mutated masks, and writable
snapshots at process transport boundaries. Velocity's float32 `asarray` is
often a view for native float32 input; negative-band clamping allocates only
when needed. `dataset[()]` need not make an extra full-sized NumPy copy.

Potentially avoidable costs include HDF5 temporary slabs before replay writes,
full spectral replay before its one-pass reduction, and displacement's additional
float32 frame replay. The latter supports stable reference/range/registration
passes and must not be removed without checking their lifetime and range logic.
Displacement replay still scans all frames for its sampled-range behavior even
when retaining only a prefix.

Source proxy objects survive in `velocity_source` and velocity-analysis source
models across workflow descendants. `InputClient.cache` holds mappings until
client close; workers close their own attached mappings. Export jobs receive
persisted snapshots and can outlive their numerical invocation. No shared-memory
allocation may be recycled merely because a predicted DAG last use has passed.

## Chunk/cache implications

The current 16 MiB stage reads four spectral frames, versus 64 frames per
storage chunk; a complete temporal chunk plane is 256 MiB across 1024 spatial/
frequency chunks. The observed 8 MiB raw chunk cache cannot hold that plane.
The default spectral working slab is 16 frames and also splits these chunks.
Moment staging reads 16 frames against 32-frame chunks. Cache eviction and
repeated chunk access are plausible costs, distinct from hash collisions and
`phil` contention. These exports are **uncompressed**, so repeated decompression
cannot explain their timing. Requested logical bytes do not establish physical
chunk/SMB bytes: unfiltered HDF5 can optimize partial I/O and the OS can cache it.

A 256 MiB spectral block aligns 64-frame chunks but increases buffer/scratch
pressure; benchmark before choosing it. Cache tuning is relevant to repeated
partial chunk access, not inherently to one full read. Do not globally enlarge
rdcc based on file size.

The h5py interpreter-wide `phil` serializes HDF5 calls, including calls to
different files. It permits overlapping HDF5 input with non-HDF5 computation;
a reader thread is therefore a valid candidate. Final HDF5 publication is a
separate phase. See [h5py threading](https://docs.h5py.org/en/stable/threads.html),
[dataset reads/chunks](https://docs.h5py.org/en/stable/high/dataset.html), and
[raw chunk-cache controls](https://docs.h5py.org/en/stable/high/file.html).

## Benchmark method

`scripts/benchmark_hdf5_input_io.py` runs the actual `SpatialReducer` on identical
source values, a fixed nontrivial gap support, unchanged reductions and identical
gzip-9 reduced-output HDF5 publication. It compares existing staging, experimental
direct-to-mapping staging, full RAM (`[()]` and `read_direct` separately), threaded
full load, spawned shared-memory full load, synchronous slabs, one-block-ahead
threaded slabs and whole-HD-file local copying. Full-load async modes have no
artificial compute head start; they measure first-use cost, not the maximum
possible lookahead benefit in a different schedule. Each trial gets a fresh
interpreter and private local scratch; trial order is deterministically shuffled.

Every variant's reduced values are checked exactly, including NaNs, against
the first variant's stored reference. Nonfinite and compressed synthetic cases
are additionally checked against the production reducer in focused tests.
Threaded slabs publish independently owned arrays; advancing or closing the
iterator cannot overwrite a retained prior block. The spawned reader opens its
own HDF5, reads into a coordinator-owned dedicated shared allocation, and sends
only control messages. Failed/incomplete reads never publish readiness; the
reader is joined/terminated before unmapping. No buffer recycling is used.

Metrics include total through reduced publication, HDF5 call spans, first block,
input wait, compute, interval-intersection overlap, logical bytes/calls, logical
scratch-file size, sampled process-tree RSS, reader startup and publication.
RSS is sampled at about 20 ms and summing shared mappings across processes can
double count pages. It excludes the system filesystem cache; scratch-file size
is not physical SSD write volume. Read time includes cache/transport/HDF5 work;
it is not isolated decompression or CPU time. There are no measured SMB byte or
physical disk counters. No cache purge, server-cache control or per-core trace
was performed. Fresh processes do not imply cold OS caches.

`scripts/benchmark_hdf5_workflow.py` pins an exported settings JSON, resolves
through `run_service`, executes the real DAG, joins exports and completes final
publication. Its optional direct-staging prototype patches only HDF5 dataset
staging inside that benchmark process; ndarray/product storage remains unchanged.
It records actual staged dataset names/sizes/times and final HDF5 publication.
It never changes active user settings or writes source-adjacent results.

## Results and decision

Correct-branch measurements and the final architecture decision are recorded
below after the focused runs. Earlier observations are provisional: an initial
NAS full read took approximately 28 seconds, whereas a later NAS slab pass after
whole-file copying took approximately 3 seconds. That is direct evidence of a
cache-state confound; it is not a 9x architecture speedup or a cold-cache trial.

## Reproduction

Use the repository environment and a new local output folder for every command:

```powershell
$env:PYTHONPATH = 'src'
$py = '.\.venv\Scripts\python.exe'
$localHolo = 'D:\Data\Matlab vs Python\260616_ABS1850789_R_1.holo'
$nasHolo = 'Z:\260616\260616_ABS1850789_R_1.holo'
$localHD = 'D:\Data\Matlab vs Python\260616_ABS1850789_R_1\260616_ABS1850789_R_1_HD\h5\260616_ABS1850789_R_1_HD_output.h5'
$nasHD = 'Z:\260616\260616_ABS1850789_R_1\260616_ABS1850789_R_1_HD\h5\260616_ABS1850789_R_1_HD_output.h5'
& $py scripts/benchmark_hdf5_input_io.py --source $localHD --output .tmp/io-ssd-new --repeats 3 --workers 4
& $py scripts/benchmark_hdf5_input_io.py --source $nasHD --output .tmp/io-nas-new --repeats 3 --workers 4
& $py scripts/benchmark_hdf5_input_io.py --source $nasHD --output .tmp/io-aligned-new --modes slab-sync slab-thread --buffer-mib 256 --repeats 3 --workers 4
# Use the same exported settings.json for each full-workflow comparison:
& $py scripts/benchmark_hdf5_workflow.py --holo $localHolo --settings settings.json --output .tmp/workflow-ssd-baseline-new --strategy npy
& $py scripts/benchmark_hdf5_workflow.py --holo $localHolo --settings settings.json --output .tmp/workflow-ssd-direct-new --strategy npy-direct
& $py scripts/benchmark_hdf5_workflow.py --holo $nasHolo --settings settings.json --output .tmp/workflow-nas-baseline-new --strategy npy
& $py scripts/benchmark_hdf5_workflow.py --holo $nasHolo --settings settings.json --output .tmp/workflow-nas-direct-new --strategy npy-direct
& $py -m pytest test/test_hdf5_input_benchmark.py -q
```

Run sequentially to avoid benchmark-induced contention. Repeat workflows in
alternating order with separate destinations; retain `timing.json`, settings and
outputs. Identify cache labels explicitly. A reboot/controlled cache experiment
is separate work; do not describe these commands as evicting any OS/NAS cache.
Use `--rdcc-nbytes` and `--rdcc-nslots` only for a targeted cache experiment.
The component benchmark has a full-array admission limit (`--max-ram-mib`);
exceeding it fails explicitly instead of silently changing strategies.
