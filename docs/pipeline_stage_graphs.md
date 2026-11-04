# Executable stage graphs

Pipeline packages declare an immutable `StageGraph`; the engine plans and runs
it through the managed invocation path. Graph factories can be inspected without
execution. Spectral, displacement-map and displacement-profile runners reject
direct calls with `Pipeline 'NAME' requires the EyeFlow managed execution runtime.`

```python
from pipeline_engine.stage_graph import ProductSpec, StageGraph, StageSpec

def calculate(ctx):
    ctx.dag.handoff("metric", compute(ctx.dag.require("source")))

graph = StageGraph((
    StageSpec("source", load_source, products=(ProductSpec("source"),)),
    StageSpec("metric", calculate, requires=("source",), products=(
        ProductSpec("metric", pack=pack_metric, exposed=True),
    )),
))
```

Callbacks return `None` and persist products through `ctx.dag.handoff`.
`ProductSpec.path` and `pack(value, ctx)` are alternative output declarations;
logical product IDs are independent of HDF5 paths. `graph.plan()` exposes ordered
stages, producers and dependency queries. Invalid definitions fail preflight.

| Dependency | Required behavior |
|---|---|
| `requires` | Product must be accepted; later producer failure does not invalidate it. |
| `optional_requires` | Wait for producer resolution; absence permits `ctx.dag.get` fallback. |
| `after` | Named stage must finish successfully. |

Missing mandatory handoffs fail the producer. Missing required inputs skip the
consumer; `ProductSpec(optional=True)` allows omission and `ctx.dag.skip(reason)`
records expected absence. Independent stages continue after failure.

## Ownership and outcomes

`DataManager` snapshots products into `ProductStorage`; readers share immutable
arrays, while `require_mutable` makes a private copy. Mutable result envelopes
are private per read. `InvocationResource` holds local readers/buffers that close
at graph completion and cannot cross process boundaries.

A handoff validates and packs the product, accepts its output checkpoint, then
records availability and diagnostics. Later failures roll back pending output,
state bindings and owned artifacts while retaining earlier commits. State
changes must rebind through `ctx.state.set`; in-place mutations are not tracked.
See [output ownership](data-contracts.md#output-ownership-and-schema) for dataset
immutability, shared groups and metadata owners.

| Invocation outcome | Meaning |
|---|---|
| `succeeded` | Required work completed; optional absence or legitimate no-output callbacks are allowed. |
| `partial` | Required work was blocked, but independent products were accepted. |
| `skipped` | Required work was blocked without an own accepted product. |
| `failed` | A stage raised or omitted a mandatory handoff. Earlier products remain usable. |

## Artifacts and diagnostics

Ready `artifact=True` stages may start during handoff. `wait_for_artifacts=True`
waits for declared ancestors' exports and skips the consumer if an export fails.
Existing `BackgroundArtifacts` owns rendering; numerical stages remain serial.

Graph scope enforces context identity, owning thread and lifetime. `StageJournal`
owns root `stage_status`, `stage_failures`, `stage_skips` and `stage_products`,
merging qualified invocation scopes without overwriting other pipelines.
Invocation-local reader adapters and scratch intermediates use `InvocationResource`;
they remain local and do not serialize full working arrays at each stage.
Durable datasets still pass through the authoritative output transaction. Only
explicitly exposed scientific products may cross process boundaries.

Production factories and callbacks live together in each package's `runner.py`:

| Pipeline | Calculation flow and optional branches |
|---|---|
| `low_freq_spectra` | source, sampling, calibration, spatial reductions, beats, QC, figures |
| `displacement_map` | prepare replay; independent artery/vein registration; finalize shared fields |
| `displacement_profiles` | analysis, profiles; optional segment maps |
| `topology_core` | inputs, gap support, shared sampling topology |
| `spatial_gradient_moment0` | inputs, cached gradient profiles, lumen figures; persist only when targeted |
| `absolute_waveform_metrics` | shared waveform inputs, global/selected segment metrics; optional quadrants reuse those metrics |
| `waveform_shape_metrics` | shared waveform inputs, global/selected segment metrics; optional quadrants reuse those metrics |
| `lowrank_waveform_decomposition` | shared waveform inputs, decomposition; optional independent quadrant decompositions |
| `blood_volume_rate` | inputs, selected gradient/mask branches, corresponding figures |
| `velocity` | selected inputs, all estimates, method-specific cardiac cycles, workflow publication |
| `velocity_analysis` | inputs, optional segments, per-beat analysis, continuous/shared waveform publication; selected maps/profiles/fits/FFTs/quadrants and artifacts |
| `pdf_report` | one atomic report render after upstream artifact completion |

Profile FFT accumulation stays fused with segment extraction to avoid replaying
the velocity volume; its selected output stage publishes the retained result.
Velocity method failures retain their established per-workflow outcome policy.
See [process execution](parallel_pipeline_execution.md) for global readiness,
transport and resource lifetime.
