# Executable stage graphs

A pipeline returns an immutable `StageGraph`. Its package declares calculations,
products, packers and artifacts; `StageGraphPlanner` resolves the local graph
without constructing global pipeline descriptors. `StageGraphExecutor` runs
numerical stages serially and continues independent branches after a failure.

```python
from pipeline_engine.stage_graph import ProductSpec, StageGraph, StageSpec

def calculate(ctx):
    values = compute(ctx.dag.require("source"))
    ctx.dag.handoff("metric", values)

graph = StageGraph((
    StageSpec("source", load_source, products=(ProductSpec("source"),)),
    StageSpec("metric", calculate, requires=("source",), products=(
        ProductSpec("metric", pack=pack_metric, exposed=True),
    )),
))
```

Callbacks return `None`. A product may be internal, an explicit
`InvocationResource`, one direct `path` with metadata, a `pack(value, ctx)`
mapping, or an exposed product for another pipeline. `path` and `pack` are
mutually exclusive. Logical IDs are independent of HDF5 paths.

`graph.plan()` exposes `order`, `stage_by_id`, `producer_by_product`,
`dependencies(id)`, `dependents(id)` and `ancestors(id)`. Kahn ordering preserves
declaration order among ready stages. Invalid IDs, producers, dependencies,
cycles and conflicting direct output paths fail preflight.

Dependency meanings remain distinct:

- `requires`: an accepted product, even if its producer subsequently fails.
- `optional_requires`: wait for acceptance or definitive producer termination;
  absence is allowed and `ctx.dag.get(name, default)` supplies a fallback.
- `after`: successful completion of the named stage.

`ProductSpec(optional=True)` permits omission. Missing mandatory products fail
their producer; missing required inputs skip their consumers. `ctx.dag.skip`
records expected missing data. Every read must be declared by that stage.

## Acceptance and ownership

`ctx.dag.handoff` calls the executor through a narrow stage interface. The
executor validates the scope/declaration, calls scientific validation, asks
`DataManager` for a stable handle, calls the scientific packer, and uses existing
output checkpoints/writers. Only after the output payload is accepted does
`ProductStore` record availability and `StageJournal` record the commit.
Artifact dispatch and exposed-product notification follow.

`ProductStore` holds identities, producer associations, availability and handles.
It performs no copying, storage, HDF5 writes, IPC or artifact dispatch.
`DataManager` recursively snapshots scientific values using `ProductStorage`.
Large arrays receive one immutable `.npy` snapshot; consumers reuse read-only
maps. Compatible managed backing is reused. Read-only views of writable bases
still need a snapshot. `require_mutable` explicitly creates a private writable
copy. Mutable dataclass envelopes are private per read without copying arrays.
Small string-array metadata retains its original dtype and shape.

Stage checkpoints and the invocation's outer checkpoint advance at handoffs.
Later failures roll back unfinished edits while preserving accepted records,
state bindings and owned artifacts. Packing/writing failures expose no product.
Graph stages use handoff rather than direct managed HDF5 writes.

Explicit local resources use `InvocationResource(value)`. Their `close` method
runs at graph completion. They are never transported; use them for readers and
replay buffers, not cross-pipeline data. Dense displacement fields instead live
in the run-owned product directory and survive replay and worker termination.

## Artifacts, scope and diagnostics

`ArtifactDispatcher` launches ready `artifact=True` stages during handoff,
including before a numerical producer returns. Dispatch is deduplicated and
reentrancy controlled. It reuses `BackgroundArtifacts`; there is no new render
pool. `wait_for_artifacts=True` waits for declared ancestors' jobs, releases
compute slots during the wait, and skips a consumer if an export failed.

`GraphScope` validates invocation identity, active stage, thread ownership and
lifetime, and restores context activation on exit. Copied contexts, other
threads and late handoffs cannot borrow an active graph.

`StageJournal` owns `stage_status`, `stage_failures`, `stage_skips` and
`stage_products`. Persisted states remain `pending`, `running`, `succeeded`,
`failed`, `skipped`. Reconciliation merges qualified invocation scopes and
individual stage/product records, preserving other pipelines and variants.

The spectral graph in `pipelines/low_freq_spectra/graph.py` resolves to
`source, sampling, calibration, spatial, beats, qc, figures`. Registration and
profile declarations live in `pipelines/displacement_map/graph.py`.

See [cross-process execution](parallel_pipeline_execution.md), the
[baseline audit](execution_ownership_audit.md), and the
[implementation report](execution_ownership_report.md).
