# Pipeline-owned stage graphs

A pipeline can return a `PipelineGraph` instead of a result mapping. Its own
package defines the graph, calculations, products, packers and export callbacks.
The engine supplies execution and uses the existing `PipelineDAG` planner;
scientific definitions do not belong in the engine or the global catalog.
Legacy function pipelines remain supported with their whole-invocation checkpoint.

`pipelines/low_freq_spectra/stages.py` is the first production example. Its
source, sampling, frequency calibration, spatial reduction, beat products, QC
and figure nodes have separate dependency and commit boundaries.

## Handoff instead of large return objects

```python
from pipeline_engine.imports import PipelineGraph, PipelineNode, PipelineProduct

def calculate(ctx):
    metric = compute_metric(ctx.dag.require("spectrum"))
    ctx.dag.handoff("metric", metric)
    # Additional work can follow. A later exception does not undo metric.

def pack_metric(value, ctx):
    paths = EyeFlowOutputPaths.active()
    return {
        f"{paths.spectrograms_masked_metrics_root}/F50":
            (value, {"unit": "kHz", "dimDesc": []}),
    }

node = PipelineNode(
    "metric", calculate,
    requires=("spectrum",),
    outputs=(PipelineProduct("metric", pack=pack_metric),),
)
```

The pipeline returns `PipelineGraph((source_node, node, figure_node))`.
Stage functions return `None`; they publish products explicitly. A deeply nested
helper can import `handoff` or `current_pipeline_dag()` instead of accepting
the whole context. These resolve an invocation-local `ContextVar`, not a singleton.

A product may specify `path` and `attrs` for a single dataset, or `pack(value,
ctx)` for a mapping of schema paths to existing `DatasetValue`/`(data, attrs)`
payloads. Packing remains scientific pipeline code. Shared HDF5 writers retain
their serialization, units, axis metadata, compression and provenance policies.
An optional `validate(value)` callback owns product-specific numerical checks;
the engine does not change NaN, sign or boundary policies.

On handoff the engine validates the declared producer, snapshots the value,
packs and atomically writes its datasets, advances the checkpoint for those
paths, flushes the HDF5, and makes the logical product available. Duplicate
handoffs and overlapping dataset paths within a graph are rejected. Failed
packing leaves no partial product, including when the caller catches its error.
Managed direct HDF5 writes from graph stages are rejected; use handoff.

Handed values and values returned by `require()` are copied so subsequent
caller/consumer mutations cannot change another stage's product. Use handoff
for reduced calculated products, not lazy HD datasets or full source videos;
typed input handles remain on the owning execution context.

## Dependencies and partial success

- `requires` names products that must have been handed off.
- `optional_requires` orders a node after those producers, but allows it to run
  if they fail or omit a product. Read them with `ctx.dag.get(name, default)`.
- `after` requires successful completion of another named node, useful when
  completion itself matters rather than an earlier partial handoff.
- `PipelineProduct(optional=True)` permits a stage to omit that product.
- `ctx.dag.skip(reason)` records expected missing input without failing the run.

Read access is restricted to a node's declared products. The planner validates
missing producers, duplicate names, self-dependencies and cycles before running
the graph. A failed producer can still have committed products. Consumers of
those products continue; consumers of missing products skip. Independent nodes
also continue. A failed stage reports an incomplete parent pipeline after safe
nodes finish, preserving committed data and exports.

Root JSON attributes `stage_status`, `stage_failures` and `stage_skips` record
outcomes per pipeline and execution variant. Existing external pipeline
dependencies remain whole-pipeline dependencies; this change does not silently
relax an unmigrated consumer's requirements.
`stage_products` records each handed-off name, producer and actual HDF5 paths,
including products committed by a stage that subsequently failed.

## Artifact nodes and ownership

Mark a node `artifact=True` to run it as soon as its declared products are ready,
including during the producing function's handoff. It can call shared figure
writers or `ctx.defer_artifacts(...)` with a pipeline-specific snapshot exporter.
Only immutable snapshots go to worker processes, never live HDF5 handles.
Exports start after the corresponding commit. A later calculation failure
cannot cancel them or remove their files.

Shared PNG, EPS, paired article figures and image writers capture exports
automatically during managed execution. They preserve existing paths and
formatting. Figure construction stays with the producer; rendering/encoding
runs in isolated Agg processes. Spectral exports use a reduced HDF5 snapshot
so their figure construction also happens in the worker.

Consumers such as reports declare `wait_for_artifacts=True`. Within a stage
graph, declare the relevant product/`after` dependencies as well; the consumer
joins those upstream export jobs. Top-level `pdf_report` joins the file's jobs.
Other numerical nodes do not wait for exports.

OutputManager and shared writers record explicit artifact ownership. Checkpoint
rollback removes only files reserved by its uncommitted invocation/segment.
There is no directory scan inferring ownership from creation time. Direct file
writes that bypass managed writers/reservations cannot be automatically owned;
custom writers must call `register_artifact_path(path, writing=True)`.

The batch reuses at most two rendering processes and permits analysis of the
next file while exports of the preceding file finish. At most two unfinished
inputs are retained to bound snapshot accumulation. The asynchronous completion
log identifies the input and producer. A file's final completion log and the
batch return wait for its exports and final metadata. Stop-between-files stops
new analysis and drains exports for inputs already started.

## Velocity and execution scope

Every file and pipeline invocation has a new graph execution. The moments and
band-ratio variants also get separate contexts, product dictionaries, status
records and artifact namespaces. Logical HDF5 paths pass through the existing
`PipelineH5Output` namespace mapping, so a `ProcessingMoments/...` declaration
is mapped to `ProcessingBandRatio/...` in that workflow.

Copying a context clears its graph handle. A copied context cannot publish into
the original graph, and a handle from an inactive/different invocation is
rejected. Handoffs execute on the owning pipeline thread; compute threads return
their results to it. HDF5 commit and scheduling are never performed by rendering
processes. Tests cover partial handoffs in one velocity workflow while the
other workflow succeeds with the same product names.

## Verification

Focused tests cover atomic and duplicate handoffs, invalid dependency graphs,
committed data surviving later failures, early export dispatch, report barriers,
artifact ownership, mutation isolation and velocity-variant separation. The
full dual-workflow run with PNG/EPS exports and PDF reports also passes.

The representative synthetic `(640,64,128,128)` HD fixture was replayed through
the graph. All **149 scientific datasets**, including 127 metrics, were exactly
equal in dtype and values to the prior reduced reference. The engine additionally
writes its ordinary `app_versions` provenance dataset. One local run took
11.34 seconds without figures, including engine initialization and checkpoints;
that is not directly comparable to earlier isolated calculation-stage timings.
The result and summary are under
`.tmp/spectral_performance/stage_graph_640_final/`.

```powershell
.\.venv\Scripts\python.exe benchmarks/benchmark_spectral_stage_graph.py `
  --hd .tmp/spectral_performance/representative_640/synthetic_HD.h5 `
  --reference .tmp/spectral_performance/cpu_parallel_640/profile.h5 `
  --output .tmp/spectral_graph_new
```

The benchmark uses a synthetic mask/cardiac timing and suppresses figure export;
it does not measure a patient's full workflow or NAS performance. The frozen
launcher supports multiprocessing dispatch, but an installer rebuild was not
part of this verification.
