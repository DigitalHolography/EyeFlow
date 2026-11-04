# Execution ownership audit

## Baseline, before implementation

Checkout: `refactor/pipeline-execution-ownership`,
`01b2223fc087624b2b3246f3a524f040f3ddbe79`.
Four existing untracked HDF5 investigation/benchmark files are outside this change.

Registration (`low_freq_spectra.__init__`) creates a descriptor. `PipelineDAG`
selects it, `CoordinatorScheduler` admits the invocation, and `invocation.py`
constructs a context/checkpoint. The runner returns `PipelineGraph`.
`GraphExecution` adapts stages to temporary global descriptors, executes stages,
deep-copies handoffs, packs scientific payloads, writes the managed output,
advances checkpoints, and recursively dispatches ready artifacts. Workers call
the same invocation code. Worker events serialize cumulative accepted snapshots;
`ProcessExecutor._finish` accepts only the last snapshot, then finalization
joins exports and publishes through `h5_publication.py`.

| Existing operation | Existing authoritative implementation |
|---|---|
| Stage definitions | Pipeline `stages.py`; closures for displacement |
| Stage planning | `PipelineGraph.plan` adapting to `PipelineDAG` |
| Stage execution/readiness | `GraphExecution.run`, `_run_node`, `_dispatch_artifacts` |
| Logical products, snapshots, commits, diagnostics | `GraphExecution.handoff`, `_record`, `require` |
| Physical payload storage | `ProductStorage`; separate writable scratch allocator |
| Descriptor transaction | `invocation._pipeline_checkpoint`, `OutputModel.checkpoint` |
| Transport/lifetime | `worker.send`, `ProcessExecutor._receive`, `collect`, `close` |
| Coordinator acceptance | `process_executor.accept_result`, `_finish` |
| Invocation scheduling | `CoordinatorScheduler` |
| Per-input/variant workflow | `runtime._execute_input` |
| Artifact rendering/publication | `BackgroundArtifacts`, `ArtifactWorkers` |
| Final output lifetime/publication | `PendingOutput`, `publish_output_model` |

Array ownership: calculation owns a writable NumPy array; handoff deep-copies
it; the managed writer snapshots the array to `.npy` for large payloads; each
stage `require` deep-copies the handoff again. Worker snapshots reference
managed arrays but repeatedly pickle cumulative records/state. Writable state
arrays can incur another `.npy` per snapshot. Coordinator reads references from
the run-owned parent directory; downstream processes receive state snapshots
only after producer completion. Final serialization reads retained references.
The input service stages each requested dataset once per canonical source ID.

Failure traces:

- A stage failure restores uncommitted records/bindings/artifacts. Earlier
  handoffs advance nested and outer checkpoints and survive. Product consumers
  can run after their producing stage fails.
- A worker crash preserves its latest received snapshot on run-owned disk, but
  it becomes authoritative only in `_finish`; no early global readiness exists.
- Acceptance validates record/artifact collisions, then uses the existing model
  checkpoint and restores state/files on failure. No earlier snapshots have
  been independently accepted, so a rejected final snapshot loses that worker's
  otherwise valid earlier work.
- Publication failure leaves the old final file intact and cleans destination
  staging. `PendingOutput.finish` releases storage in `finally`; it occurs after
  invocation consumers, not while scientific consumers are active.

Baseline command (repository `.venv`): the 13 files requested in section 22,
from `test_pipeline_graphs.py` through `test_run_service.py`: **124 passed,
20 failed, 6 subtests passed, 136.11 s**. All 20 failures are spectral scheduling
regressions: `visibility="hidden"` and `auto_run_after="velocity"` are commented
out in the checked-out spectral declaration. This explains missing spectral
outputs and missing expected spectral failures. No graph/process/resource
baseline failures were observed. A separate parallel coordinator run passed
14 tests; the resource-budget file passed 2 tests.

Synthetic baseline: one writable 16 MiB float32 array, one persisted dataset,
three declared stage consumers. Instrumented `graphs.deepcopy` and enumerated
actual staging files: **4 array deep copies, 67,108,864 copied bytes, 1 `.npy`
write, 16,777,344 `.npy` bytes**. Three repeated cumulative transport snapshots
were 929 bytes each (array data already referenced). Stage elapsed 0.0411 s;
this is diagnostic, not a speed claim. Representative acquisitions/NAS and
before-refactor process-tree RSS/acceptance latency were not measured.

## Intended corrections

Dedicated local planner; separated execution, products, scope, journal and
artifact readiness; DataManager over the existing physical storage; explicit
product identities/requirements; incremental worker deltas accepted on the
coordinator thread before invocation completion; run-owned backing lifetime
through consumers/publication. Existing transaction and scientific writers
remain authoritative.
