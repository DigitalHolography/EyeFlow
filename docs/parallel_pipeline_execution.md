# Pipeline execution and cross-process product exchange

`runtime.py` coordinates one input and its variants. `CoordinatorScheduler`
chooses invocations; `ProcessExecutor` creates/supervises workers; each
invocation uses the same checkpoint and `StageGraphExecutor` path. The local
stage planner and global `PipelineDAG` have distinct semantics.

## Product-level readiness

Pipeline descriptors can declare a `graph_factory(options=...)` and explicit
`ProductRequirement` inputs. The same factory defines exposed products during
preflight and builds the executable graph. Requirements select producer closure
in the global DAG but are separate from successful whole-invocation edges.

```python
ProductRequirement(
    alias="registration", producer_pipeline="displacement_map",
    product="registered_fields", variant="shared",
)
```

`PipelineInvocationId(run_id, pipeline, variant)` plus a product name is the
global identity. `shared` means no workflow qualifier; `same` uses the consuming
workflow's identity. Runtime replaces pending shared placeholders with workflow
identities when a variant producer publishes them.

The registry distinguishes `PENDING`, `AVAILABLE` and `FAILED`. Scheduler
admission waits for all product requirements to resolve, as well as whole
dependencies and CPU/RAM admission. A failed product does not skip an entire
graph: stages requiring it skip, and independent/optional branches can run.
Producer termination finalizes every unproduced exposed product. Ordinary
whole-pipeline dependencies still require successful invocation completion.

Displacement profiles consume `registered_fields` through a product requirement
and retain whole dependencies on topology and velocity workflows. Their
registration relationship no longer waits for the entire producer invocation.

## Receipt, acceptance, publication

1. The stage establishes stable ownership and locally accepts its product/output.
2. The worker closes an immutable **incremental** commit manifest in the
   coordinator-owned run directory and sends a typed event containing stable
   references and small metadata. No large array goes into an IPC queue.
3. The receiver forwards logs, queues events and returns `TransportReceipt`.
   This is receipt only; it does not make a product globally available.
4. On the coordinator thread, `ProductAcceptance.accept_event` validates the
   producer, commit, storage shape/dtype/owner, record merge and artifacts. The
   run already controls backing lifetime, so acquiring it requires no array copy.
5. The existing output checkpoint transaction accepts records, state bindings
   and artifacts. Only then does the registry make exposed handles available.
6. Scheduler readiness can now launch a downstream worker while the producer
   continues. No producer HDF5 file is required.
7. Finalization eventually joins exports and publishes the accepted model with
   `h5_publication.py`.

The worker may continue after receipt because neither worker exit nor its
resource cleanup removes the run-owned snapshots. Authoritative acceptance is
serialized on the coordinator thread, never in the receiver. No lock is held
across coordinator calculations. Receipt remains responsive during those
calculations; acceptance resumes when the coordinator returns to scheduling.

Output revision cursors send changed accepted records only. State is a binding
delta; artifacts are sent once. Final snapshots do not reapply prior datasets.
Identical commit retries are idempotent; conflicting retries fail. Journals
merge through `StageJournal`, not an independent transport implementation.

## Failure and lifetime

Accepted products/records survive stage failure, later acceptance rejection and
worker crash. Uncommitted edits remain invisible. A failed acceptance restores
its record/state transaction and removes files created by that acceptance.
The producer reports failure; its accepted product consumers may still succeed.
Whole-invocation consumers retain their stricter requirements.

`ProductStorage` is owned through `PendingOutput`: backing stays alive until
consumers and artifact jobs finish and publication completes or fails. Workers
close mappings and reclaim accounts but do not delete run-owned backing.
Cancellation terminates workers, finalizes pending products and releases
reservations. Publication uses local serialization, a destination staging file,
and atomic replacement; failure cannot expose a partial new final file.

One `ResourceBudget` supplies CPU/RAM admission and process-tree RSS observation
for calculations, workers, renderers and overlapping file finalizers. There
are no fixed per-pipeline RAM declarations. Numerical stages remain serial.
Native thread policy remains in resource/process initialization.

The input service opens each HD/DV source once, stages requested datasets in
bounded reads, and reuses source identities/hard links across workers. No new
NAS cache or prefetch framework is introduced.

## Verification

`test_product_exchange_processes.py` uses real spawned processes and actual
storage/transport/acceptance. Event proxies and process-exit waits prove early
consumption, crash recovery, delayed acceptance, duplicates, collisions, invalid
references, acceptance exceptions and cancellation without timing sleeps.
Existing graph, resource, workflow, displacement and writer tests preserve the
scientific schemas. `benchmark_execution_ownership.py` records copy/storage,
manifest/IPC, acceptance, publication and RSS counters. Representative NAS and
acquisition comparisons remain separate measurements.

The [implementation report](execution_ownership_report.md) lists actual owners,
source organization, measurements, verification and limits.
