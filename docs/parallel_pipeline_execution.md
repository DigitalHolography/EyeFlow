# Pipeline execution and product exchange

`runtime.py` owns an input run and its variants. `CoordinatorScheduler` chooses
ready invocations; `ProcessExecutor` owns workers and reaping. Both process and
coordinator invocations use the same checkpoints and
[stage executor](pipeline_stage_graphs.md).

## Readiness

A descriptor's `graph_factory(options=...)` supplies the same product declarations
for preflight and execution. A `ProductRequirement` adds producer closure to the
global DAG without requiring successful whole-invocation completion:

```python
ProductRequirement(
    alias="registration", producer_pipeline="displacement_map",
    product="registered_fields", variant="shared",
)
```

Product identity is `(run_id, pipeline, variant, product)`. `shared` selects the
unqualified producer; `same` selects the consuming workflow. Variant identities
are bound when the workflow producer publishes them.

The registry tracks `PENDING`, `AVAILABLE` and `FAILED`. Admission waits for all
required products to resolve, whole dependencies, and CPU/RAM availability.
Failed products skip dependent stages while independent branches run. Producer
termination fails its unproduced products. Whole-pipeline dependencies still
require successful completion. Displacement profiles use early registration
products while retaining whole dependencies on topology and velocity workflows.

## Delivery and acceptance

1. A stage locally commits its product and output.
2. The worker closes an incremental manifest in run-owned storage, then queues
   its reference and small metadata. Large arrays stay outside IPC queues.
3. The receiver queues the event and returns a correlated `TransportReceipt`.
   Only this receipt advances the worker's delivery cursor; it does not enable
   downstream readiness.
4. The coordinator's `ProductAcceptance` validates producer/commit identity,
   storage, output paths and artifacts. One existing output transaction accepts
   records, changed state bindings and files before publishing registry handles.
5. Ready consumers may start while the producer continues. Finalization joins
   exports and publishes HDF5 through `h5_publication.py`.

Serialization, delivery or receipt failure ends the worker explicitly. Already
accepted commits survive, including a commit whose receipt was lost. Incomplete
manifests are removed; staged backing stays run-owned until cleanup. Replayed
commit IDs are idempotent; conflicting replays fail. Status-only events avoid
repeated scientific snapshots. State deltas compare bindings against their last
acknowledged value, including rebindings to original input values.

## Resource and output lifetime

Accepted products survive later failures; rejected acceptance rolls back its
records, state and newly copied artifacts. `PendingOutput` retains backing until
consumers, exports and publication finish. Worker cleanup closes its mappings
without deleting shared backing. Publication uses a destination staging file and
atomic replacement, preserving the prior output on failure.

Logical completion retains CPU/RAM reservations until actual process exit and
nonblocking reaping. Each reservation releases once. Cancellation terminates
workers; cleanup exceptions report failure; repeated close is safe. One
`ResourceBudget` covers workers, calculations, renderers and file finalizers.
There are no fixed per-pipeline RAM declarations.

The input service opens each HD/DV source once and reuses bounded local replay
across workers and source aliases. Native thread policy belongs to process and
resource initialization.

`test_product_exchange_processes.py` covers actual transport, crash recovery,
receipt loss, acceptance rollback and exit/resource ordering. Reproduce synthetic
copy, manifest, IPC and RSS measurements with
`python scripts/benchmark_execution_ownership.py --processes`.

## Performance diagnostics

Run logs use `[PERF]` records with `perf_counter()` timestamps and durations in
seconds. Pipeline and stage records identify the pipeline, execution variant,
and stage; stage timings distinguish callback work from CPU admission,
dependency dispatch, and artifact waits. Handoff records split snapshots,
packing, and output commits. Aggregate storage metrics report source staging,
array copying, serialization, and final HDF5 publication.

Background artifact records separate snapshot preparation, admission/renderer
queueing, rendering, file copies, and publication. `stage_timings` and
`artifact_timings` remain the authoritative persisted timing data. The run-end
summary sorts the slowest pipeline/stage/artifact spans and estimates accounted
wall time as the interval union of stage and artifact spans; nested and
concurrent spans are therefore not added together. The remainder is reported
as unaccounted wall time, commonly including setup and finalization.

The shared `Logger` shows `DEBUG` and `PERF` records by default in development,
and only `INFO`/`WARNING`/`ERROR` records in frozen builds. Set
`EYEFLOW_DIAGNOSTICS=1` to enable diagnostic records in a packaged build;
`EYEFLOW_PERF_LOGGING=1` independently enables performance records;
`EYEFLOW_LOG_LEVEL` accepts `DEBUG`, `INFO`, `WARNING`, or `ERROR` as a minimum
severity. Startup configuration is inherited by spawned workers and artifact
renderers, so filtering is consistent across the run and UI.
