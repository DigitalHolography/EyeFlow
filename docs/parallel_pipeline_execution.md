# Coordinator and process execution

`pipeline_engine.runtime` owns a run's resolved DAG, workflow variants, input
service, output structure, accepted stage products, resource budget and final
publication. GUI and CLI continue to enter through `run_service`.

## Scheduling

Declarations may opt into `@pipeline(..., execution_mode="process")`; the default
is `"coordinator"`. Process entry points must be importable module-level registered
pipelines. Dependency-ready work starts without a barrier between graph levels.
Displacement is the first production process pipeline and runs once, outside
velocity workflow fan-out. Small waveform metrics stay sequential. Workflow
variant descendants currently retain coordinator execution, even if marked as
process work; variant scheduling and state isolation keep their established
semantics during this first rollout.

The coordinator can run a calculation while independent processes execute. Its
input-service and result-receiver threads remain responsive during that work.
Only the coordinator merges results or changes the authoritative output model.
Input-service and result queues are bounded. Queue messages contain small events
and local result-file references, rather than videos or large result objects.

## Inputs and local storage

The coordinator opens the HD/DV sources once and reuses their typed adapters,
schema resolution, sidecar configuration and canonical alignment. Workers receive
dataset adapters, never open HDF5 handles. A request for an uncached dataset is
handled by the input service. Numeric datasets are staged in bounded sequential
reads into a local `.npy`; subsequent slices and other pipelines reuse that file.
This first implementation stages the complete requested dataset, not individual
arbitrary slices. Datasets never requested are not copied.

Workers and coordinator calculations see read-only input data. A calculation
requiring writable input must make its own bounded copy. Dataset proxy identities
remain stable within a run, so pure calculation caches can share results across
velocity workflows.

`EYEFLOW_SCRATCH_DIR` selects a local scratch directory. Network locations are
rejected. Otherwise the storage manager prefers a writable local SSD on Windows
(using the device seek-penalty query when available), then the normal local
application-data/temp location. A configured local directory takes precedence.
Media classification can be unavailable, in which case the default local path
is retained. Large arrays use `.npy`; snapshots use `.pickle` or `.h5`; artifacts
use their ordinary suffixes. No engine temporary file uses `.tmp`.

Scratch exhaustion is reported explicitly. It does not trigger an unbounded RAM
allocation or fallback to NAS scratch. Source data, retained products, and worker
directories remain local until their consumers finish. File mappings and worker
processes close before their owning scratch directory is removed.

## Simple commits and failures

`OutputModel` holds groups, attributes, payload references and HDF5 storage
options in memory. Numeric payloads at least 64 KiB are stored in local `.npy`
files. Existing scientific serialization policies (dtype conversion, compression,
NaN chunk storage, mirrored FFT virtual datasets and provenance) remain owned by
the shared I/O writers and are applied when serializing the final file.

Managed checkpoints journal only changed output records. Handoff commits advance
the relevant dataset paths and their parent attributes; later failures discard
uncommitted records. They do not visit/snapshot the whole output hierarchy.
Standalone callers using real HDF5 handles retain the existing HDF5 checkpoint
adapter. State checkpoints track bindings, not arbitrary in-place edits to nested
mutable objects; committed products must not be mutated in place.

Each process has its own output structure and private state. It publishes only
changed state bindings, avoiding replacement of newer coordinator state with an
old input snapshot. No Git repository, content hashing, conflict reconciliation or
large version-history metadata is created. Known workflow namespaces determine
placement; duplicate dataset paths raise a programming error.

The event protocol is:

```text
WorkerStarted
StageStarted
ProductCommit
StageFinished / StageFailed / StageSkipped
InvocationFinished
```

A product commit references a closed local snapshot containing stable products,
state changes, owned artifacts and the stage summary. The result receiver
registers and acknowledges it before the worker continues. This acknowledgement
also prevents loss of a queue feeder's buffered handoff on abrupt process exit.
The final invocation message returns all completed work, including partial
results after a handled failure. If a process dies, the coordinator retains the
last acknowledged stable snapshot and records its exit failure. Unpublished
work cannot be recovered.

The existing intra-pipeline DAG API continues to distinguish required products,
optional products and successful-stage dependencies. A failed stage does not
cancel consumers of products it already committed. Other nodes and exports
continue. Displacement prepares one replay and then has separate artery/vein
nodes: a missing or failed vessel does not erase the other vessel's AVI or field
state. Legacy pipelines are initially one transactional stage. Existing
inter-pipeline dependencies still require whole-pipeline success.

## Resources

One budget spans a batch, including overlapping file completion and renderers.
It supersedes the former independent fixed limits in `runtime_limits.py`; that
module now bootstraps native libraries and delegates worker allowances.

The default compute allowance is 75% of currently unused logical CPU capacity,
respecting CPU affinity. External load is estimated from utilization, excluding
EyeFlow's process tree, and smoothed for subsequent admissions. Software thread
counts are not subtracted from logical CPU counts. At least one compute slot is
allowed. Running native operations are not interrupted when the allowance falls.

The run RAM reservation limit is 50% of available memory at budget creation. It
is not repeatedly halved as EyeFlow allocates memory. Shared counters separate
transient calculations, retained outputs, input data and transfers. Retained allocations leave working-memory headroom for downstream consumers.
Worker-owned
reservations are reclaimed after a crashed process exits. Large retained arrays
spill locally. Internal Python pools respect the invocation allowance, while
native libraries are configured to avoid multiplying Python worker counts.
Process declarations may use `native_parallelism=True` (displacement does): their
native libraries receive the allowance and Python pools are capped to one. Other
workers give the allowance to Python pools and keep native libraries at one.
Export admission uses the same budget, with separate dispatch and publication
threads. Artifact barriers release compute reservations while waiting.

These are admission limits and allocation accounting, not an OS-enforced hard
limit on arbitrary native-library allocations or mapped file pages. A monitor
tracks process-tree RSS quietly; one observed peak RAM summary is logged when
the batch finishes. Scheduling delays are logged as events. The initial process
and renderer working-memory reservations are estimates
(up to 64 MiB each, bounded by one quarter of the run budget); representative data is needed to tune
them. Helper/IPC/monitor threads are not active compute slots.

## Reports and final publication

Reports run after numerical work and required artifact publication. They read the
assembled output through `output_source`, without reopening a future result
filename. Existing disk-HDF5 report callers remain supported.

The coordinator serializes the accepted structure to a local `.h5` after exports
finish, then copies it to a destination-side staging name ending in `.h5` and
renames it to the final filename. There is no extra read-back verification pass.
An incomplete run still publishes its accepted datasets and failure records.
Batch analysis may overlap the preceding input's exports/final serialization;
the final path does not exist until that input's publication finishes. Required
payload storage is retained until serialization completes.

## Verification

`test_parallel_coordinator.py` exercises real Windows spawn, simultaneous
independent work, one-read input replay, immutable local array transfers, partial
stage failure, abrupt worker exit, reservation reclamation and the production
displacement pipeline. Existing stage-graph, failure-isolation, writer, spectral,
dual-workflow, report and batch-export tests cover the adapters and consumers.

Concurrency correctness does not establish a NAS speedup. Measure representative
runs through final HDF5 publication, including input bytes read, local scratch
usage, peak RAM, native thread utilization and all exports.
