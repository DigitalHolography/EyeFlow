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
| `velocity` (hidden) | Physical retinal-velocity estimation, whole-vessel signals, average maps, and cardiac-cycle timing; produces `velocity` and `cardiac_cycles` |
| `topology_core` (hidden) | Canonical artery/vein prepared topology; produces `prepared_topology` |
| `velocity_analysis` | Shared per-beat/segment/profile state plus selectable persisted continuous, quadrant, profile-fit, FFT, map, and movie products; requires velocity and topology |
| `spatial_gradient_moment0` | Independently selected moment0 gradient/lumen profiles; requires cardiac cycles and topology, not velocity analysis |
| `blood_volume_rate` | Two option families with different DAG dependencies: gradient edges and mask-derived geometry |
| `waveform_shape_metrics`, `absolute_waveform_metrics`, `lowrank_waveform_decomposition` | Downstream waveform metrics requiring `velocity_analysis` |
| `pdf_report` | Report assembly after waveform and shape outputs |
| `displacement_map` | Separate image-registration path producing magnitude videos and run-local dense-field artifacts |

The declarations in each package are authoritative if this table drifts.

## Shared scientific flow

`vessel_inputs.py` builds canonical `RetinalSourceData`: active HD volumes,
aligned vessel masks, optic-disc geometry, timing/pitch, DV background settings,
and velocity method. The velocity and velocity-analysis pipelines use compatible
sources.

The velocity pipeline estimates the physical velocity field and detects
zero-based cardiac-cycle boundaries. Native geometry lives in
`calculations/topology`; the spatial `SegmentSamplingPlan` is prepared by
`calculations/vessel_segments/sampling` and cached by `topology_core` in run
state. The existing `prepared_topology` DAG key is unchanged.
Velocity analysis consumes both results and produces the common
per-beat, segment, and profile state used by visible products.

The velocity estimator is product-specific and therefore remains in
`velocity/estimation.py` rather than `calculations`. It processes `(frame, y,
x)` volumes in bounded temporal chunks, uses moments or calibrated frequency
bands as its active source, applies the shared local-background difference, and
returns physical velocity in `mm/s`. Exact-zero LF samples map to zero in band
mode; finite/non-negative validation and calibration provenance are part of the
pipeline contract.

`spatial_gradient_moment0` shares cardiac-cycle/topology state but owns its
gradient profiles and lumen edges. Do not put gradient work into velocity
analysis.
It requires the optional HD flat-field moment at `/moment0ff` or `/M0FF`, and
publishes under `/Processing/SpatialGradientProfiles` and
`/Processing/SpatialGradientMetrics` only when it is a direct target. When it
runs solely as a blood-volume-rate dependency, keep the products in `ctx.state`
without persisting those intermediate HDF5 families. Its branch lumen plots
come from the computed `(time, branch)` metric, not a separate plot-only
calculation.
`blood_volume_rate` retrieves the products required by each enabled option and
must keep segment labels, centers, branch/radius dimensions, and pixel scale
aligned.

## Velocity semantics

All selected downstream targets run for both automatic velocity workflows. `doppler_moments`
and `frequency_bands` represent physical velocity (`mm/s`); band mode first
converts `HF / LF` to frequency using the recorded
`band_ratio_frequency_scale_hz`. `velocity/semantics.py`
centralizes display and dataset interpretation. Moments publish under
`/Processing`, bands under `/ProcessingAlt`, with isolated downstream state and
`moments`/`bandratio` artifact folders. The single velocity pipeline declares
itself as an execution-variant producer and publishes one typed
`ExecutionVariant` per successful method; the engine fans out its declared DAG
dependents. Each workflow detects its own cycles from raw RMS frequency before
estimation and uses them for its downstream analyses; segmentation remains
under `/Segmentation`. New velocity-derived outputs or plots must resolve
semantics from payload/provenance instead of hard-coding `mm/s`.

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
| Velocity source/semantics | `vessel_inputs.py`, velocity sources/runner/estimator, data contracts | displacement internals |
| Per-beat/segment output | velocity-analysis builder and packers, output schema, profile/segment tests | settings UI |
| Spatial-gradient lumen metric | spatial-gradient package, vessel-segment sampling/measurement, image filters, spatial-gradient/topology tests | CLI and report code unless output is exposed there |
| Blood-volume-rate formula | option declaration, runner, outputs, `calculations/blood_volume_rate.py`, BVR tests | unrelated waveform metric calculators |
| Metric family | that pipeline's runner/calculator/outputs plus waveform input contracts and matching test file | pipeline engine unless dependencies/options change |
| Displacement | `displacement_map/` and `test_displacement_map_pipeline.py` | topology and waveform metrics unless the dormant segment helper is activated |
| Profile fit | `calculations/vessel_segments/profiles/fits/quadratic.py`, then `velocity_analysis/outputs/profile_analysis.py`, dedicated doc/test | GUI and settings |
