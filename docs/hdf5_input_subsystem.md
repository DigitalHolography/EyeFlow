# Declarative HDF5 inputs, backend policy and continuous prefetch

The coordinator owns one `InputService`, its authoritative source-reading thread,
and read-only HD/DV handles. Scientific scheduling, calculations, axes, product
acceptance and output publication are unchanged. `auto`, `direct`, and `replay`
remain public comparison modes. End-to-end time and consumer blocking take
precedence over minimizing RAM or scratch use.

## Access declarations and representation

Runtime reuses the selected pipeline dependency closure and stage graphs built
for product registration. `configure(graphs, run_settings, execution_order=...)`
resolves semantic fields through the typed HD/DV schemas, canonicalizes hard-link
aliases, filters `when` by the selected velocity methods, and records shape,
dtype, chunks, compression, required and potential consumers. It never evaluates
scientific callbacks. Invalid/missing candidates remain `unresolved`; stages own
existing branch validation, fallbacks and failure isolation.

`InputUse` describes metadata, small arrays, sequential slices, streaming blocks,
or random slices. Random includes full-time-axis spatial selections. A required
consumer count is **not** a physical pass count. Potential consumers also include
optional uses; random access indicates possible repeated regions even for a single
stage. Optionality affects speculative eligibility separately from representation.
No exact pass count is invented.

The deterministic `auto` policy is:

| Source/access | Representation | Logged reason |
|---|---|---|
| Local SSD | Worker-local direct HDF5 | `local_storage` |
| Unknown | Coordinator direct, or worker direct when explicitly configured | `unknown_fallback` |
| Network, fixed-size numeric input >=4 KiB with random access or multiple potential consumers | Progressive local replay | `potential_repeated_regions` |
| Network, single-use stream >=1 MiB | Bounded shared-memory buffers | `single_pass_bounded_prefetch` |
| Network, sequential input >=1 MiB, including one consumer | Progressive local replay | `sequential_prefetch` |
| Other small/undeclared inputs | Direct; normal bounded streaming remains available through `iter_blocks` | `small_or_undeclared` |

`bulk_bytes` exposes the 1 MiB threshold for experiments. Shape/dtype reject empty,
scalar or object inputs from bulk planning. Chunks and compression are logged and
influence physical read geometry. `expected_consumers` describes required graph
uses; `potential_consumers` and `reuse` describe risk, not measured accesses.
This is a conservative policy, not a calibrated throughput predictor.

`direct` forces worker-local reads on local storage and bounded coordinator
transfers on NAS storage. It creates no mandatory `.npy` input staging. Declared
stream consumers can still use bounded asynchronous buffers, including before
dispatch. `replay` forces reusable local backing; streaming consumers then copy
ready regions into exclusive working buffers. Metadata-only declarations never
initiate bulk prefetch. No dataset-name allowlist chooses a backend.

### Why moment0 and moment0ff differ from the old policy

Topology sampling calls `resolve_segment_rotations` separately for artery and
vein plans. For each unresolved segment it scans reference frames; the helper's
`_frame_slice` requests complete spatial frames before extracting that segment.
The number of these passes depends on actual unresolved geometry. Velocity then
reads its own temporal chunks. Topology's optional random moment0 declaration
therefore expresses real potential reuse that the retired required-count-only
policy ignored. The new policy stages those remote source regions once.

`temporal_median_image` intentionally computes the exact finite-sample temporal
median in approximately 8 MiB float64 spatial tiles. With roughly 640 frames,
that produces sides around 40 pixels. These are calculation working sets, not a
reason to issue the same remote selections. The gap input declaration now uses
random access. NAS replay reads temporal source slabs, then serves spatial tiles
from the ready local mapping. The median algorithm and its numerical policy are
unchanged. An optional-only input starts no speculative source transfer.

Replay and stream windows align down to temporal HDF5 chunk boundaries when the
chunk fits the requested window and byte bound. Oversized chunks retain bounded
partial reads. Slabs span the spatial axes; oversized frames use bounded tiles.
Native dataset handles and their HDF5 chunk caches survive between operations.
Chunk-intersection estimates are not physical cache-miss or NAS-byte counts;
actual decompression/cache effects still require measurements.

## Scheduling is independent of representation

At source-operation boundaries, queued demand and consumer lookahead win. Before
starting queued future work, the planner can admit higher-priority active prefetch;
it does not accumulate additional same-priority jobs behind that work. With an
empty queue, it selects the next useful replay **or streaming** operation. It uses
pipeline order plus local stage topological order, with current activation and
scheduler-supplied viability updates. Computational dependencies do not prevent
independent external input reads. Prefetch never authorizes early computation.

| Priority | Work |
|---|---|
| 0 | Immediate demand, including promotion of an existing job |
| 1 | Current consumer lookahead, paired windows and attached stream filling |
| 2 | Active-pipeline prefetch |
| 3 | Future-pipeline prefetch |

Within a priority class, pipeline/stage first use wins **before** unused-region
count. Sources in the same first-use group are balanced by unused ready regions,
with paired inputs adjacent and declaration order breaking ties. A distant empty
stream therefore cannot reserve the allowance before nearer replay sources advance
through additional windows. Current paired lookahead still precedes the next pair;
stale lookahead is discarded.
A running bounded HDF5 read is not interruptible. Queue draining and priority
promotion take effect before the next source operation.

For example, during topology gaps computation, moment0ff demand stays ahead of
moment0/moment2/band speculation. The latter are eligible through required selected
velocity consumers even though moment0 is optional in topology. A future spectra
stage can have its initial bounded stream buffers ready too. On velocity dispatch,
remaining reads become active priority; slice requests promote current regions.
The exact order follows declarations and configured capacity, not an acquisition ID.
The topology gaps stage declares its present flat-field video as required for
that work, allowing imminent active prefetch before a blocking request. An absent
video still follows `_prepare_gap_masks`' existing legacy-mask fallback; optional
sampling inputs remain optional. The selected spectral spatial stage declares its streaming input as required
for that work; its source/inventory stage still handles absent or legacy-group
videos with the existing skip path. Marking both declarations optional previously
excluded a valid present spectral video from all future planning.

### Future stream ownership

The planner creates an unopened pool through the existing `_open_stream` and
`_fill_free` lifecycle with no client owner. It fills only the existing bounded
number of buffers. Actual `stream_open` adopts a compatible canonical pool and
promotes queued fills, preserving ready data and its reservation. Smaller blocks
still satisfy caller upper bounds. Incompatible stricter requests discard the
unopened pool and open a correctly bounded one.

Once attached, the normal FREE/FILLING/READY/IN_USE protocol applies. A leased
buffer cannot be overwritten until release. Mutable median partitioning changes
only exclusive working storage. Demand can evict unopened speculative pools to
obtain stream capacity or input memory; attached pools and leases are never evicted.
A canonical input gets at most one speculative stream attempt per run, so closing
a completed/cancelled pool cannot repeatedly prefetch its beginning. Later explicit
consumers can always open a fresh stream. Eviction may discard bounded unused data;
it deliberately does not create an unbounded retry/read loop.

`CoordinatorScheduler` supplies viable consumers after outcomes. Completed pipelines
and failed/skipped dependency closures lose speculation; independent consumers
retain it. A failed variant does not remove data still needed by other invocations
of that pipeline. Cancellation removes unowned pools and queued speculative replay,
while preserving urgent demand and completed reusable regions.

The obsolete activation-only first-region scheduler and required-count-only
backend decision were removed. There is one source thread and one job scheduler.

## Readiness, limits and APIs

| InputPolicy setting | Default |
|---|---:|
| Replay window | 32 frames |
| Stream window / pool size | 16 frames / 2 buffers |
| Source operation byte bound | 64 MiB |
| Active plus unopened stream pools | 2 |
| Queued speculative replay/unopened-stream jobs | 16 |
| Unconsumed planned replay regions per source | 16 |
| Unconsumed replay bytes plus unopened pool reservations | 1000 MiB |
| Large sequential/stream threshold (`bulk_bytes`) | 1 MiB |
| Future prefetch | Enabled |

The planner generates work progressively, rather than building volume-sized queues.
Limits constrain additional planned work; urgent demand and normal consumer lookahead
retain priority. Full conversions queue one missing replay region at a time. Stream
pools shrink buffer counts/frame counts under memory pressure and fail explicitly
when even one frame cannot fit. Speculative memory refusal backs off for 20 ms;
control traffic and demand remain responsive. Existing input-account headroom and
memory reservations remain authoritative. Mapped pages/OS/HDF5 caches are not an
application RAM cap.

Replay regions are never exposed until fully read and flushed. Production still
flushes every I/O unit, including large-frame tiles. Workers map only completed
extents; full references require every region ready. Aliases share readiness and
ready regions are never rewritten. Read errors preserve prior completed regions
and propagate without falsely publishing success. `.npy` files retain their full
logical shape; the speculative-byte bound limits useful reads/writes, not logical
file allocation. Local placement and free space are checked before staging.

Use `dataset.read_slice(...)`, indexing, `read_direct`, or read-only `np.asarray`.
For exclusive mutable working blocks:

```python
with dataset.iter_blocks(max_bytes=64 * 1024**2, writable=True) as blocks:
    for lease in blocks:
        with lease as block:
            process_in_place(block)
```

Block views must not escape their lease. Workers close/reap clients before the
reader stops; unfinished replay files, owned buffers and reservations are cleaned
up before source handles close. Forced worker termination uses the same reclamation.
Local `read_direct` fills supported caller destinations natively. Safe numeric dtype
conversion occurs in HDF5; unsupported selectors, strided destinations and unsafe
narrowing casts preserve the previous NumPy fallback and nonfinite/overflow policy.

Set `$env:EYEFLOW_FUTURE_PREFETCH='0'` or `InputPolicy(future_prefetch=False)` for
comparisons. Active planning/lookahead remain. `lookahead=False` disables automatic
planning. `cancel_prefetch()` also suppresses immediate planner refill, globally or
for a canonical source.

## Startup attribution and diagnostics

The old batch path called `HoloRunLayout.hd_h5` and `.dv_h5` for logging and again
for dispatch. Each call rescanned directories and checked HDF5 signatures. The
resolved `RunRequest` now snapshots these validated paths before batch execution;
HDF5 itself is still opened/validated at runtime. There is no global discovery cache
that could silently stale UI selections. Directory enumeration now uses `scandir`
entry metadata and filters suffixes before file/signature checks.

Development-only `Logger.perf_stage` boundaries cover DAG resolution, input discovery
and validation, companion-path snapshot, budget creation/activation, output directory
reset, runtime output preparation, product/run-array scratch, input thread startup,
graph construction and scheduler setup. Existing HDF5-open, output initialization,
input metadata/manifest, pipeline dispatch and publication timings remain. Discovery
and plan validation occur before `run/batch`; output reset and runtime setup occur
between `run/batch` and `input_hdf5_open`.

`input_dataset` logs location, bytes, access, required/potential consumers, reuse
classification, chunks/compression and the decision reason. `input_source_hdf5_read`
logs source regions, reason/priority, enqueue/start/end times and logical bytes;
`input_local_replay_mapping` identifies local mappings separately. Enqueue, selection,
promotion and cancellation records explain scheduling. Aggregate metrics include
source bytes/durations/counts per dataset, flush/write time, stream adoption/reuse,
consumed/unused prefetch, queue waits and input memory peaks.

`input_idle_exhausted_seconds` measures queue-wait idle without more eligible work.
The existing `input_idle_resource_seconds` remains inclusive of policy/resource
limits; new `input_idle_ahead_limit_seconds` and `input_idle_constraint_seconds`
split that total into read-ahead/pool limits and genuine memory/storage constraints.
These counters exclude planning, mapping and flush work. Transition-only
`input_idle_begin`/`input_idle_end` records show the reason, held replay bytes,
unopened-stream reservations, blocked next-region indices and applicable limits.
Their wall-clock spans include intervening control handling; the counters measure
only actual queue waits. No extra timer or polling loop was added.

`input_region_ready` exposes the next missing region for reevaluation at the next
operation boundary. `input_plan_selection` records first use, balanced depth and
pair identity; `input_prefetch_consumed` names the speculative regions or stream
frames later delivered. No new consumer request is necessary to refill planned
work: the existing source loop replans after each completed queued operation. I/O/computation
overlap uses unions of measured intervals, excluding each consumer's input waits.
Histories are bounded to 8,192 intervals; dropped samples invalidate the overlap
result instead of fabricating one. Logical HDF5 bytes/time are neither actual NAS
wire bytes nor isolated network-transfer duration. Ordinary executable logging does
not acquire performance traces.
