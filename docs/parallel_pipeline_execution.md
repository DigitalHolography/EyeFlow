# Pipeline execution design

## Purpose and scope

Commit `6aa4bfd8` introduced process execution, local input replay, resource
admission, partial-result recovery, and deferred HDF5 publication. Its other
changes (velocity workflow selection, gradient memoization, CUDA packaging and
release configuration) are independent features and retain their behavior.
The stage graph API predates this commit; its scientific contract remains in
[stage graphs](pipeline_stage_graphs.md).

The implementation must let independent displacement work overlap coordinator
calculations without reading source videos repeatedly or losing successful
products after a failure. GUI and CLI still use `run_service`; calculations,
HDF5 paths, axes, units, workflow namespaces and pipeline selection do not change.

## Ownership

| Component | Owns | Must not own |
|---|---|---|
| `runtime.py` | Per-input lifecycle and workflow fan-out | IPC, file-copy transactions |
| `invocation.py` | Context construction, one descriptor and its checkpoints in either mode | Batch scheduling |
| `scheduler.py` | Dependency readiness, artifact barriers, progress order | Worker serialization, output merging |
| `process_executor.py` | Spawned workers, event receipt/acknowledgement, result acceptance, shutdown | Workflow fan-out, scientific work |
| `worker.py` / `protocol.py` | Spawn entry point and explicit invocation/event/result messages | Authoritative run state |
| `input_service.py` | Source metadata, lazy staging and read-only input adapters | Output serialization |
| `storage.py` | Immutable payload references and mapping lifetime | Hardware selection, DAG policy |
| `resources.py` | Shared CPU/RAM admission and process memory accounts | Pipeline dependencies |
| `output_model.py` | Output tree, touched-record journals and collision checks | Pipeline scheduling |
| `input_output/writers/h5_publication.py` | HDF5 serialization and destination publication | Output ownership and scheduling |
| `finalization.py` | Export join, final publication and payload-store lifetime | Input or worker lifetime |

These are concrete runtime services, not interchangeable backend frameworks.
Workers import the invocation executor directly; they do not import the run
coordinator to execute one descriptor. The scheduler delegates process work
through a small executor interface and has no knowledge of HDF5 records.

## Execution and failure contracts

1. Resolve settings, options, dependencies and destinations before execution.
   Start ready process descriptors before selecting ready coordinator work.
   Discovery order breaks ties. Pipeline-name dependencies require successful
   whole invocations; failed dependencies skip only their consumers.
2. Process execution is opt-in (`execution_mode="process"`) with an importable
   registered entry point. Base displacement registration runs once independently;
   its profile/map process waits for registration, topology and shared timing.
   Missing vessel fields are skipped; one successful vessel keeps registration
   successful and permits profile analysis. Artifact consumers,
   variant producers and their descendants execute on the coordinator.
   `run_once` stays outside fan-out. Artifact consumers run after numerical jobs
   and required exports.
   Process pipelines consume sources and declared state, not earlier pipelines'
   output trees; pipelines reading assembled output stay on the coordinator.
3. Each invocation owns private state bindings, output edits and new artifacts.
   Checkpoints journal changed records only. Stage handoffs advance declared
   products and parent attributes. In-place mutation of existing state objects
   or published arrays is outside the checkpoint contract.
4. Workers send small typed events over a bounded queue. Product, terminal-stage
   and final events reference closed local snapshots containing records, changed
   state bindings and owned artifacts. LogMessage events are forwarded immediately
   by the receiver, including during long coordinator calculations, rather than
   buffered in snapshots. The receiver acknowledges
   receipt before the worker continues, even while the coordinator calculates.
   Stage events describe graph nodes; invocation events never impersonate stages.
   A process crash retains the last acknowledged snapshot and records failure.
   Collect worker completions before yielding coordinator work when CPU/RAM
   admission is full.
5. Acceptance validates output collisions and every artifact destination before
   exposing changes. A failed acceptance restores model/state changes and removes
   only files created by that acceptance. An earlier producer's data survives.
   Worker state is a binding delta, so an unchanged old binding cannot replace
   a newer coordinator value. Deleting shared bindings is unsupported.
6. Continue independent jobs after calculation, spawn, transfer or export failure.
   Publish accepted partial results and failure metadata before reporting an
   incomplete run. Always join/terminate workers and release their reservations.

## Data and resource lifetime

The coordinator opens each HD/DV file once. A requested numeric dataset is staged
once by source identity, in bounded sequential reads, then replayed read-only
from local storage. Hard links share that replay; unrequested datasets are not
copied. Metadata and small nonnumeric values use the existing typed adapters.

Large immutable arrays travel as local `.npy` references. Repeated references
within and across snapshots reuse persisted bytes and preserve aliases when
loaded. Writable arrays are copied at transport boundaries; identity alone must
never cache a mutable value. No full video travels on an IPC queue.

`EYEFLOW_SCRATCH_DIR` may select a writable local directory; otherwise use the
OS local temporary directory. Reject network scratch and report allocation
failure. Do not scan other drives or silently spill to the source/output NAS.
Input mappings close after workers finish. Product storage stays alive through
final publication and closes before its directory is removed.

One batch-wide budget covers calculations, workers, renderers and overlapping
finalization. The default CPU limit is 75% of affinity-available logical CPUs
(at least one); an explicit limit is deterministic. RAM reservations are capped
at 50% of available memory sampled at budget creation. Retained arrays leave
working headroom and spill locally when needed. Admission uses the larger of
explicit reservations and observed process-tree RSS above the run baseline.
Both displacement processes declare 512 MiB of working memory; when the budget
cannot grant that amount, reservations retain consumer headroom and profile
chunking uses the granted allowance. These are admission estimates,
not hard limits on native allocations or mapped pages.

The budget exposes admission operations; callers do not inspect private locks or
hardware samples. Coordinator native pools stay serial; invocation allowances
control Python pools. A worker configures native pools once at startup, assigning
its allowance either to native libraries or Python pools. Overlapping threads
must not temporarily change and restore process-wide native limits. Artifact
barriers release compute reservations while waiting. Failed workers' owned
memory reservations are reclaimed only after exit. Peak process-tree RSS is
reported once at batch completion.

## Completion

Reports read the accepted output model and published artifacts. Finalization
joins exports, writes their statuses, serializes through the existing HDF5 writer
policies, copies to a destination-side staging file, and renames it to the final
name. The final path is absent until publication succeeds. Finalization owns
its model/storage explicitly and may overlap the next input's analysis; batch
completion waits for every finalizer. Serialization failure must not return a
missing output as successful.

## Verification

Focused tests cover readiness/failure propagation, real spawn and abrupt exit,
atomic acceptance, immutable transfer reuse, one-read source replay, resource
reclamation, deferred publication, report reads and production displacement.
Existing workflow/stage/writer tests protect scientific and serialization
contracts. Performance claims require representative runs measured through final
publication; concurrency tests alone do not establish a NAS speedup.
