# Execution ownership implementation report

## Scope and baseline

Implemented on `refactor/pipeline-execution-ownership`, starting at
`01b2223fc087624b2b3246f3a524f040f3ddbe79`. The
[initial audit](execution_ownership_audit.md) records the original flow, owners,
failure guarantees, measurements and baseline failures before implementation.
Existing untracked HDF5 investigation/benchmark files were left untouched.

## A. Final responsibility matrix

| Operation | Actual owner and entry points |
|---|---|
| Scientific declarations | Pipeline `graph.py` factories and `stages.py` callbacks |
| Local graph validation/planning | `StageGraphPlanner.plan`, immutable `StagePlan` |
| Stage lifecycle/readiness | `StageGraphExecutor.run`, `_run_stage` |
| Logical local availability/access | `ProductStore.record_accepted`, `validate_read`, `require` |
| Stable handles/snapshot ownership | `DataManager.accept`, `read`, `writable`, `acquire` |
| Physical persistence/mapping primitive | One retained backend: `ProductStorage.store`, `store_array`, `_write_array`, `read`, `close` |
| Handoff coordination | `StageGraphExecutor.handoff` |
| Stage/descriptor transactions | `invocation._pipeline_checkpoint`, `OutputModel.checkpoint`, `ProductCheckpoint.commit/rollback` |
| Artifact readiness/barriers | `ArtifactDispatcher.dispatch`, `ready`, `wait` |
| Rendering/publication jobs | Existing `BackgroundArtifacts` / `ArtifactWorkers` |
| Invocation/thread/stage scope | `GraphScope.activate`, `executing`, `validate` |
| Stage diagnostic representation/merge | `StageJournal.mark`, `committed`, `persist`, `reconcile`, `observe`, `finalize` |
| Processes/receipt/supervision | `ProcessExecutor.start`, `_launch`, `_receive`, `collect`, `_finish`, `close` |
| Authoritative cross-process acceptance | `ProductAcceptance.accept_event`; `accept_result` implements its output/artifact transaction |
| Local exposed-product acceptance | `ProductAcceptance.accept_local`, after local output acceptance |
| Global product availability | `ProductExchangeRegistry.declare`, `accept`, `fail`, `finalize_producer` |
| Invocation readiness/admission | `CoordinatorScheduler._products_ready`, `_start_processes`, `coordinator_jobs`; shared `ResourceBudget` |
| Per-input/variant coordination | `runtime._execute_input`, `_bind_product_requirements` |
| One descriptor/outer checkpoint | `invocation._run_pipeline_descriptor`, `_execute_pipeline` |
| Managed records/revision cursors | `OutputModel`, `committed_delta`, `ProductCheckpoint` |
| Final lifetime | `PendingOutput.finish/close` |
| Physical HDF5 publication | `h5_publication.publish_output_model` |

`ProductStorage` is the existing single physical backend, reused by DataManager,
managed output, input replay and legacy transport. DataManager delegates array
snapshots to it; there is no second array-copy or `.npy` writer implementation.
Writable scratch allocation remains distinct from accepted immutable storage.
Acceptance and journal reconciliation were removed from `ProcessExecutor`.

## B. Resulting source organization

```text
src/data_manager/
  handles.py          immutable handles, metadata and explicit local resources
  ownership.py        structured snapshots, read access, adoption, resource cleanup
src/pipeline_engine/
  stage_graph/
    model.py          immutable executable definitions
    planner.py        dedicated local Kahn planner and immutable plan
    executor.py       stage lifecycle and coordinated handoffs
    products.py       logical local product state/access validation
    artifacts.py      early readiness, deduplication and completion barriers
    journal.py        statuses, product diagnostics and scoped reconciliation
    scope.py          thread/context/stage lifetime enforcement
    exceptions.py     expected scientific skips
  product_exchange/
    identity.py       run/pipeline/variant/product identities and requirements
    registry.py       global terminal dependency states and dependents
    acceptance.py     authoritative incremental commit acceptance
  metrics.py          thread-safe ownership/transport counters
  storage.py          retained physical payload implementation
  output_model.py     existing transaction model plus accepted revision cursors
  protocol.py         typed product events, manifests and receipt acknowledgments
  worker.py           spawned invocation and incremental manifest creation
  process_executor.py transport and process supervision
  scheduler.py        whole-invocation and product-level readiness
  runtime.py          input lifecycle, service assembly and variant qualification
src/pipelines/{low_freq_spectra,displacement_map}/
  graph.py            authoritative declarations used in preflight and execution
  stages.py           scientific callbacks
  runner.py           pipeline entry adapters
```

Model/planner imports do not depend on the global DAG or process scheduler.
DataManager imports no stage execution logic; product exchange imports no
scientific pipeline packages. Static import-direction checks and a test that
forbids construction of global descriptors during local planning pass.

## C. Real stage walkthrough

`spectral_graph()` resolves to:

```text
source -> sampling -> calibration -> spatial -> beats -> qc -> figures
```

That is the deterministic execution order, not a claim that every adjacent pair
is a dependency. Spatial requires source; beats requires spatial and optionally
calibration; QC requires spatial and optionally beats; figures requires spatial
and optionally calibration, beats and QC. The declaration retains all three
dependency categories. A product survives a later producer-stage exception;
an `after` dependency still requires that stage to succeed.

Registration declares prepare, independent artery/vein stages, and finalize.
Finalize requires at least one available field and exposes `registered_fields`.
Profiles consume that explicit product and retain topology/timing dependencies.
Analysis, profile packing and optional segment-map packing have explicit
scientific products. Replay/readers use an invocation-owned resource instead of
closure dictionaries or `True` placeholders.

## D. Cross-process product walkthrough

1. A stage calculates a writable array and calls `ctx.dag.handoff`.
2. Scope, declaration and scientific validation run. DataManager asks the
   existing storage backend for one immutable snapshot. Large arrays become
   read-only `.npy` references; compatible managed backing is reused.
3. The scientific packer produces the existing output payload. Existing nested
   and outer checkpoints accept its paths before local availability is recorded.
4. ProductStore associates the stable handle; StageJournal records the handoff.
5. A worker closes an incremental manifest and sends `ProductAcceptedEvent`.
   Exposed structured results use a stable metadata-file reference, so local
   Python types or extensive object trees do not enter the IPC queue.
6. The receiver acknowledges receipt and forwards logs. It never merges outputs
   or makes products available. Backing already resides in the coordinator's
   run-owned directory and survives worker exit.
7. On the coordinator thread, ProductAcceptance validates identity, commit,
   storage ownership/header, record merge, diagnostics and artifacts. Existing
   transactions apply the delta. Only then does the registry mark availability.
8. Scheduler admission can launch the dependent invocation before the producer
   finishes. The consumer maps the same backing read-only, with no full-array
   copy and no dependency on the producer's final HDF5 file.
9. PendingOutput retains backing until consumers, exports and final publication
   finish. Cleanup closes mappings and removes the run-owned temporary tree.

Explicit writable access makes one private copy. Mutable dataclass envelopes
are reconstructed without rerunning scientific constructors or copying arrays.
Repeated aliases within one product share one snapshot and one mapping.
Generic read-only flags alone do not authorize storage reuse: owning arrays can
re-enable writes, and external read-only maps may have writable file aliases.

## E. Failure recovery

The real spawned-process acceptance test uses manager events and process-exit
waits, without arbitrary sleeps to establish overlap. A produces an 8 MiB
array; B must start and signal A before A can crash with exit code 23. B checks
contents before and after that actual exit and completes successfully.

Assertions cover early launch, correct contents, final-HDF5 independence,
retained accepted records, absent unfinished records, no double application,
accurate stage/invocation outcomes, backing lifetime through publication,
eventual cleanup, released CPU/RAM reservations and consistent worker/coordinator
logs. Additional real-process cases cover pre-emission crash, delayed acceptance
after receipt, small exposed arrays, identical/conflicting retries, invalid
references, failed record merges, actual artifact publication failure,
acceptance exceptions, identity mismatch and cancellation. Optional missing
products and unproduced variants reach terminal dependency states.

Both production displacement pipelines also execute together in real workers,
including actual registration, profile calculations, optional segment maps and
final schema checks.

## F. Measured storage/performance comparison

One 16 MiB float32 product and three local stage consumers:

| Metric | Before | After |
|---|---:|---:|
| Full-array in-memory ownership/read copies | 4 | 0 |
| Bytes copied by those in-memory operations | 67,108,864 | 0 |
| Full-array staging copies | 1 | 1 |
| `.npy` writes | 1 | 1 |
| `.npy` bytes | 16,777,344 | 16,777,344 |
| Consumer data buffers | Separate copies | One shared buffer |
| Three manifest sizes, bytes | 929 / 929 / 929 | 982 / 662 / 662 |
| Sampled process-tree RSS, bytes | Not measured | 94,498,816 |
| Final publication | Not measured | 0.0654 s |
| NAS/source reads | No source in synthetic workload | No source in synthetic workload |

The staging copy is necessary to separate mutable producer aliases from accepted
data. Three consumers add no array copy or `.npy` write. Array-copy counters
cover ownership/storage boundaries, not arbitrary NumPy/library internals.

Separate real-process crash workload, after implementation:

| Metric | Observed |
|---|---:|
| Product size | 8,388,608 bytes |
| Array staging writes | 1, totaling 8,388,736 bytes |
| Worker-to-coordinator events | 22 |
| Pickled event payload volume | 10,660 bytes |
| Scientific worker processes | 2 |
| Accepted worker manifest serializations | 8, totaling 9,228 bytes |
| Initial coordinator state serializations | 2 |
| Exposed product emission-to-availability | 0.00902 s |
| Sampled process-tree peak RSS | 349,270,016 bytes |
| Final publication | 0.03993 s |
| Final CPU slots / reserved bytes | 0 / 0 |

The process benchmark includes one synchronization-manager process. IPC volume
describes worker events, not every queue framing byte or synchronization-proxy
message. RSS is sampled by the shared monitor, not an OS allocation trace.
Times are single diagnostic observations, not speedup claims.

Raw results: [execution_ownership_benchmark.json](execution_ownership_benchmark.json).
Reproduce with:

```powershell
.\.venv\Scripts\python.exe scripts/benchmark_execution_ownership.py --processes
```

Representative acquisition/NAS before/after throughput, baseline RSS, baseline
acceptance/publication latency and measured artifact-render latency are **not
measured**. No representative acquisition was used. Instrumentation is present
at snapshot/storage, source-read, transfer, acceptance, artifact-dispatch and
publication boundaries; `ResourceBudget` records process-tree RSS. Existing
source-replay tests assert one read per canonical dataset/hard-link identity.
Early artifact behavior is established by synchronization tests, without a
claimed rendering speedup.

## G. Removed/replaced code

- Removed `graphs.py`, `GraphExecution`, and temporary descriptor adaptation to
  the global DAG. There is no permanent legacy graph adapter.
- Removed universal product/read `deepcopy` operations.
- Removed cumulative worker result snapshots and invocation-end-only acceptance.
- Removed acceptance, artifact-copy transaction and journal merge code from
  ProcessExecutor; its accepted deltas now use one acceptance component.
- Removed displacement closure dictionaries, boolean placeholder handoffs,
  temporary-directory state exchange and duplicate direct-run execution paths.
- Restored the spectral hidden/automatic scheduling declarations whose removal
  caused all 20 focused baseline failures.

## H. Validation

All commands use `.venv\Scripts\python.exe` because default Python lacks pytest.

- Baseline requested 13 regression files: **124 passed, 20 failed**, 6 subtests.
- Same regression group after migration: **144 passed**, 6 subtests.
- Full available suite, allowing the known collection error: **1,030 passed**,
  **1 collection error**, 42 subtests. Final focused checks cover subsequent
  envelope/callback hardening: **60 passed** in the final four-file check.
- Real process exchange/failure group: **12 passed**.
- Production registration/profile process integration passes with segment maps.
- Ownership/transport/output/profile focused group after snapshot consolidation:
  **53 passed**. Earlier graph/exchange/process/profile group: **80 passed**.
- `git diff --check` and static import-direction checks pass.

The remaining collection error is pre-existing: `test_input_selection.py`
imports `_Var` from `test_stem_list_inputs.py`, where it is absent in both the
initial commit and the working tree. This UI test dependency was not changed.
No behavioral assertion was weakened: mutation tests request explicit writable
access; legacy read-only-cache tests now establish actual immutable ownership
and add checks for re-enabled writes/external writable aliases.

Primary validation commands (from the repository root):

```powershell
.\.venv\Scripts\python.exe -m pytest -q test/test_pipeline_graphs.py test/test_pipeline_failure_isolation.py test/test_parallel_coordinator.py test/test_worker_results.py test/test_output_model.py test/test_process_admission.py test/test_resource_budget.py test/test_low_freq_spectra_pipeline.py test/test_displacement_map_pipeline.py test/test_displacement_profile_outputs.py test/test_dual_velocity_workflows.py test/test_pipeline_context.py test/test_run_service.py
.\.venv\Scripts\python.exe -m pytest -q --continue-on-collection-errors
.\.venv\Scripts\python.exe -m pytest -xq test/test_product_exchange_processes.py
.\.venv\Scripts\python.exe -m pytest -q test/test_stage_graph_ownership.py test/test_worker_results.py test/test_low_freq_spectra_pipeline.py test/test_displacement_profile_outputs.py
```

## I. Remaining limitations

- No representative NAS/acquisition performance claim is made.
- Variant descendants retain the existing coordinator execution policy; their
  qualified product identities and same/named/shared resolution are tested.
- Long coordinator calculations defer authoritative worker acceptance until the
  coordinator returns to scheduling. Receipt/log forwarding remains responsive.
- Unmigrated pipelines retain whole-invocation dependencies and legacy state
  transport. They do not gain inferred product-level readiness.
- Additional pipeline migrations, NAS-aware scheduling, shared memory, automatic
  numerical-stage concurrency and refined admission prediction remain future work.

## J. Intentional implementation choices

- Retained ProductStorage as the single physical backend instead of duplicating
  or cosmetically relocating it. DataManager adds structured ownership and
  handles; existing transaction, scratch, replay and publication systems remain.
- Graph dataclasses permit positional declarations as well as keyword arguments.
  Names and execution semantics follow the new API.
- Workers use explicit receipt-only acknowledgments and run-owned backing;
  global acceptance is a separate serialized coordinator operation. A second
  acknowledgment is unnecessary for durability and is not simulated by receipt.
- Structured exposed products use an immutable metadata-file reference in
  addition to array references. This supports Python scientific types while
  keeping IPC bounded to references and small metadata.
- No shared-memory manager or new rendering pool was introduced.
