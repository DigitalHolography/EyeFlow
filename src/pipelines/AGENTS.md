# Pipeline context

## Boundary

Pipelines adapt typed inputs and shared scientific calculations to the runtime.
They declare dependencies/options, orchestrate calculations, own run-local state
keys, pack HDF5 datasets, and export product-specific artifacts. Generic DAG and
execution behavior belongs in `pipeline_engine`; reusable formulas belong in
`calculations`; external path aliases belong in `input_output/schema`.

## Discovery and package shape

`pipelines/__init__.py` scans top-level modules/packages and imports them.
Registration occurs through `@pipeline` or `@registerPipeline`; no registry list
needs editing. A helper-only top-level module must not register a pipeline.
Tutorial pipelines live under `pipelines/tutorials/` and are hidden from normal
selection.

Small pipelines may be one module. Larger pipelines use:

```text
pipeline_name/
  __init__.py       declaration only
  runner.py         orchestration and state ownership
  sources.py        typed input assembly, when needed
  outputs.py        schema packing/artifact export
  calculator.py     product-specific calculation, when not broadly reusable
```

Keep imports in `__init__.py` light enough that catalog discovery remains
predictable and reports missing optional dependencies as descriptor
availability rather than unrelated import crashes.

## Current product groups

| Group | Responsibility and important upstream state |
|---|---|
| `heartbeat_core` (hidden) | Velocity-derived cycle boundaries; caches compatible estimator work; produces `heartbeat` |
| `topology_core` (hidden) | Canonical artery/vein prepared topology; produces `prepared_topology` |
| `waveform_velocity_core` (hidden) | Shared velocity/per-beat/segment/profile state; requires both hidden cores |
| `waveform_velocity` | Selectable persisted continuous, per-beat, segment, quadrant, profile, FFT, map, and movie products |
| `spatial_gradient_moment0` | Independently selected moment0 gradient/lumen profiles; requires heartbeat and topology, not waveform core |
| `blood_volume_rate` | Two option families with different DAG dependencies: gradient edges and mask-derived geometry |
| `waveform_shape_metrics`, `absolute_waveform_metrics`, `lowrank_waveform_decomposition` | Downstream waveform metrics requiring `waveform_velocity` |
| `velocity_profile_analysis` | Weighted fits of persisted artery and vein transverse masked profiles |
| `pdf_report` | Report assembly after waveform and shape outputs |
| `displacement_map` | Separate image-registration/displacement and cross-section path |

The declarations in each package are authoritative if this table drifts.

## Shared scientific flow

`vessel_inputs.py` builds canonical `RetinalSourceData`: active HD volumes,
aligned vessel masks, optic-disc geometry, timing/pitch, DV background settings,
and velocity method. Both heartbeat and waveform core use compatible sources.

Heartbeat detection runs first because per-beat consumers need its zero-based
cycle boundaries. Topology is prepared once and cached in run state. Waveform
core may reuse heartbeat's velocity estimate when the conservative cache key
matches; it then produces common context and per-beat outputs for visible
pipelines. This reuse is an optimization, not an alternate contract.

`spatial_gradient_moment0` shares heartbeat/topology but owns its gradient
profiles and lumen edges. Do not put gradient work back into waveform core.
It requires the optional HD flat-field moment at `/moment0ff` or `/M0FF`, and
publishes under `/Processing/SpatialGradientProfiles` and
`/Processing/SpatialGradientMetrics`; its branch lumen plots come from the
persisted `(time, branch)` metric, not a separate plot-only calculation.
`blood_volume_rate` retrieves the products required by each enabled option and
must keep segment labels, centers, branch/radius dimensions, and pixel scale
aligned.

## Velocity semantics

All targets remain available for both configured methods. `doppler_moments`
and `frequency_bands` represent physical velocity (`mm/s`); band mode first
converts `HF / LF` to frequency using the recorded
`band_ratio_frequency_scale_hz`. `waveform_velocity_core/velocity_semantics.py`
centralizes display and dataset interpretation. New velocity-derived outputs or
plots must resolve semantics from payload/provenance instead of hard-coding
`mm/s`.

Velocity is never dimensionless. Reject legacy unit-`1` or relative-index
velocity metadata instead of adapting labels or output units for it.

Historically physical downstream pipeline names are not scheduling guards.
Retain clear method, calibration, and unit provenance for their results.

## State and output rules

- Give shared state keys one owning runner and a declared producer DAG key.
- Use typed objects for multi-field state and validate type/shape in consumers.
- Put durable products in the HDF5 result, not only in `ctx.state` or figures.
- Use `EyeFlowOutputPaths.active()` for centralized output families.
- Pack data with units and `dimDesc`; method-dependent velocity datasets also
  need velocity semantics.
- Options should control only their described product families. Express their
  intra-pipeline prerequisites with `PipelineOption.requires` and upstream data
  with `PipelineOption.dag_requires`.
- Do not infer whether a sibling pipeline ran from target selection; use
  `ctx.pipeline_scheduled` or retrieve declared state.

## Change navigation

| Change | Read/inspect | Usually skip |
|---|---|---|
| New pipeline or option | `CONTRIBUTING.md`, package `__init__.py`, runner, DAG tests | GUI views; discovery is automatic |
| Velocity source/semantics | `vessel_inputs.py`, heartbeat and waveform-core sources/runners, estimator, data contracts | displacement internals |
| Per-beat/segment output | waveform core plus visible waveform packer, output schema, profile/segment tests | settings UI |
| Spatial-gradient lumen metric | spatial-gradient package, topology preparation/chunks, image filters, spatial-gradient/topology tests | CLI and report code unless output is exposed there |
| Blood-volume-rate formula | option declaration, runner, outputs, `calculations/blood_volume_rate.py`, BVR tests | unrelated waveform metric calculators |
| Metric family | that pipeline's runner/calculator/outputs plus waveform input contracts and matching test file | pipeline engine unless dependencies/options change |
| Displacement | only `displacement_map/`, topology functions it imports, and displacement tests | heartbeat/waveform metrics |
| Profile fit | `velocity_profile_analysis/`, profile schema producer, dedicated doc/test | GUI and settings |
