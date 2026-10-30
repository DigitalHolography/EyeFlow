# Architecture and navigation

This document is the repository-wide map. Detailed data paths belong in
[data contracts](data-contracts.md); local change guidance belongs in the
nearest subsystem `AGENTS.md`.

## System boundary

EyeFlow starts from a `.holo` marker, resolves separately produced HD and DV
companion data, and performs downstream analysis. HD owns image sequences and
acquisition metadata. DV owns vessel/optic-disc segmentation and a small amount
of analysis configuration. EyeFlow owns pipeline selection, derived analysis,
figures/reports, and the `<stem>_EF` result directory.

The application deliberately has four layers:

```text
front ends                 `src/cli.py`, `src/ui/`, `src/eye_flow.py`
    -> orchestration       `src/pipeline_engine/`
        -> adapters        `src/input_output/`, pipeline source loaders
            -> science     `src/pipelines/` -> `src/calculations/`
        -> persistence     `src/input_output/schema/`, writers, OutputManager
```

Dependencies should point downward. `calculations` must not know about the UI,
settings files, or run-directory layout. Pipelines may combine calculations,
run state, and output paths. Front ends may select work but must not implement
scientific behavior.

## Entry points and convergence

`pyproject.toml` publishes two commands:

- `eyeflow` calls `launcher.main`. With arguments it dispatches to `cli.main`;
  without arguments it calls `eye_flow.main` and creates the desktop app.
- `eyeflow-cli` calls `launcher.cli_main` and always uses `cli.main`.

CLI flow:

```text
argparse in `cli.main`
  -> `expand_run_inputs` (file, list, recursive folder, or ZIP)
  -> settings/catalog target and option selection
  -> `resolve_run_spec`
  -> `execute_run`
```

GUI flow:

```text
`eye_flow.main` -> `EyeFlowApp`
  -> input/pipeline/settings controllers build selection
  -> `RunController._build_run_spec`
  -> `resolve_run_spec`
  -> worker calls `execute_run`
```

`resolve_run_spec` is the shared validation boundary: it validates calibration,
selectable targets, option closure, pipeline availability, DAG order,
input layouts, and unique destinations. `execute_run` owns batch iteration,
per-file cleanup, failure collection, and stop-between-files behavior.

GUI folder selection/drop and CLI folder expansion share
`input_output.runs.discovery`. It finds `.holo` markers and directories containing
matching `<stem>_HD` and `<stem>_DV` companions, including the selected directory
itself. Marker-free inputs retain their actual folder identity for progress and
errors; HDF5 validation and execution use the same run-service path.

## Pipeline model

`pipelines.load_pipeline_catalog()` discovers top-level modules from the
runtime pipeline directory followed by the built-in package path. Importing a
module executes `@pipeline` or `@registerPipeline`, which places a
`PipelineDescriptor` in `PIPELINE_REGISTRY`. There is no hand-maintained module
list.

Each descriptor distinguishes:

- `requires`: importable optional Python packages and availability;
- `dag_requires` / `dag_produces`: semantic data dependencies;
- `options`: user-selectable product families, including option-to-option and
  option-to-DAG dependencies;
- `visibility`: visible target versus hidden shared implementation;
- `input_slot`: preferred source for merged attributes, not permission to omit
  companion data from a normal `.holo` run;
- `produces_execution_variants`: whether the pipeline publishes isolated
  downstream executions;
- `auto_run_after`: automatically include a companion after a scheduled pipeline;
- `run_once`: execute once with shared state/output despite variant ancestry.
- `wait_for_artifacts`: join completed sidecars before an artifact consumer runs.

Every pipeline implicitly produces its own name. A required DAG key with no
producer is treated as external, so a misspelled key does not fail DAG
construction; consumers must still validate the external value. The resolved
plan preserves discovery order subject to topological constraints.

The important shared scientific chain is:

```text
velocity      -> `velocity` + `cardiac_cycles` + `velocity_workflows`
  -> low_freq_spectra (hidden, automatic, once per run; optional HD spectra)
topology_core -> `prepared_topology` + shared `capillary_analysis_zones`
        \          /
         velocity_analysis
          -> continuous, per-beat, segment, quadrant,
             profile, FFT, map, and artifact products
                    |
                    +-> waveform metric and report pipelines

spatial_gradient_moment0
  requires cardiac_cycles + prepared_topology
  -> spatial_gradient_edges

blood_volume_rate options
  gradient_edges -> velocity_analysis + spatial_gradient_edges
  masked_edges   -> velocity_analysis + prepared_topology
```

The hidden core pipelines are shared work, not UI targets. The declarations in
each pipeline package remain the authoritative graph; the diagram is a
navigation aid.

## Runtime ownership

For each input, `run_pipelines_to_output` creates a coordinator-owned input
service and in-memory output structure. HD/DV datasets are loaded centrally and
staged once into local read-only storage when requested. Pipeline contexts use
these dataset adapters; worker processes do not reopen the NAS source files.

The scheduler starts dependency-ready process pipelines alongside coordinator
work. `displacement_map` is the first process pipeline; small waveform metrics
remain sequential. Existing workflow variants keep isolated state and output
namespaces. Process results arrive through a bounded queue as stable local-file
references. Private stage journals preserve completed products after unrelated
failures, without restoring snapshots of the shared output file.

Shared figure writers retain their bounded renderer pool and explicit artifact
ownership. Renderers and calculations share CPU/RAM admission accounting. Reports
run last and read the assembled output view and published images. The final HDF5
is serialized locally and published after required exports complete; batch
analysis can overlap the preceding input's finalization.

See [parallel pipeline execution](parallel_pipeline_execution.md) for input
caching, queue events, stage commits, resource policy, local scratch selection
and the initial process rollout. See [stage graphs and handoff](pipeline_stage_graphs.md)
for pipeline-owned product/dependency declarations.

## Scientific processing boundaries

- `velocity` estimates physical retinal velocity, prepares whole-vessel signals,
  detects each workflow's cardiac cycles from its physical velocity after estimation,
  attempts moments and band-ratio estimation independently, and retains each
  velocity video when downstream analysis requires spatial products.
- `low_freq_spectra` runs once after velocity when HD `/spectrograms` is an
  available frequency-bin video. It owns shared spatial/gap reductions, per-beat
  spectra and filters, and spectral metrics under `/ProcessingSpectrograms`.
  Detected moments cardiac timing takes priority, with band-ratio timing as a
  fallback; metadata identifies the selected workflow. Missing timing preserves
  raw/spatial spectra. Spectral failures do not disable independent velocity
  products.
- `topology_core` aligns masks, creates annular branch/segment topology, and
  publishes reusable prepared artery and vein topology.
- `velocity_analysis` consumes velocity and topology state, performs per-beat
  and segment/profile analysis, and publishes the selected user-facing velocity
  datasets and artifacts. Shape, absolute, low-rank, blood-volume-rate, and
  report pipelines consume its declared state. The existing DAG executes these
  downstream descriptors for each successful `ExecutionVariant` using isolated
  state, `/ProcessingMoments` or `/ProcessingBandRatio`, and separate artifact folders.
  Topology and segmentation remain shared.
- `spatial_gradient_moment0` is independent of velocity analysis computation;
  it reuses cardiac-cycle timing and topology, applies its own ordered image-processing
  chain, and owns gradient-derived lumen products. Dependency-only execution
  keeps those products in run state; their HDF5 families are published only
  when `spatial_gradient_moment0` is selected directly.
- `blood_volume_rate` composes either gradient-derived edges or mask-derived
  geometry with velocity products according to its enabled option families.
- `displacement_map` is a separate registration analysis. It currently emits
  per-vessel magnitude videos and retains temporary dense-field paths in run
  state; the topology-based segment helper is not part of its runner. Do not
  assume a waveform dependency merely because it uses the same inputs.

See [pipeline context](../src/pipelines/AGENTS.md) and
[calculation context](../src/calculations/AGENTS.md) before scientific changes.

## Settings lifecycle

`default_settings.json` is the fresh-settings template. `AppSettingsStore` in
`src/app_settings.py` owns the platform-specific persisted path, loading,
normalization, import validation, and saving. Catalog-aware normalization adds
new default-selected pipelines/options to older settings without using old
files as a complete schema.

The GUI exposes settings through `ui/controllers/settings.py` and pipeline
selection through `ui/controllers/pipeline_library.py`. The CLI reads the same
store; `--pipelines` changes target selection. The JSON `velocity_estimation_method` setting selects `doppler_moments`,
`bandratio`, or `both` (the default when absent). The legacy `frequency_bands`
value is accepted as an alias for `bandratio`. Only selected workflows load
inputs, estimate velocity, and execute their downstream pipelines. The calibration comes from persisted settings and passes through
`resolve_run_spec` into `PipelineContext.band_ratio_frequency_scale_hz`.
`PipelineContext.velocity_estimation_method` identifies a workflow internally;
there is no GUI selector. Workflow selection is configured in JSON.

For a new CLI option controlling an existing setting, inspect only
`src/cli.py`, the relevant `AppSettingsStore` methods/default, the field passed
to `resolve_run_spec`, and its tests. The GUI matters only if the shared setting
contract or its presentation must also change.

Both velocity workflows reconstruct their local background with OpenCV Telea
on float32 frequency images, using a fixed 5-pixel radius. The velocity pipeline
records the method and radius in output provenance.

## Source-of-truth map

| Concept | Source of truth |
|---|---|
| Installed commands and package dependencies | `pyproject.toml` |
| GUI and CLI dispatch | `src/launcher.py` |
| Pipeline discovery | `src/pipelines/__init__.py` |
| Pipeline metadata and registration | decorators in pipeline `__init__.py` files; types in `src/pipeline_engine/base.py` |
| Dependency resolution | `src/pipeline_engine/dag.py` |
| Run validation and batch policy | `src/pipeline_engine/run_service.py` |
| Per-pipeline runtime/context | `src/pipeline_engine/runtime.py`, `src/pipeline_engine/context.py` |
| Default and persisted settings | `default_settings.json`, `src/app_settings.py` |
| Companion input layout | `src/input_output/runs/layout.py` and source layouts in `src/input_output/schema/` |
| HD/DV logical fields | `src/input_output/schema/holodoppler.py`, `src/input_output/schema/doppler_view.py` |
| Canonical aligned source models | `src/retinal_source_data.py`, `src/pipelines/vessel_inputs.py` |
| Reusable figure styling and layout | `src/plotting/article.py`, `src/plotting/formatting.py` |
| Current EyeFlow HDF5 paths | `src/input_output/schema/eyeflow_output.py` |
| HDF5 source access and output serialization | `src/input_output/readers/h5.py`, `src/input_output/writers/h5_output.py`, `src/input_output/writers/h5.py`, `src/input_output/writers/h5_storage.py`, `src/input_output/writers/h5_initialization.py` |
| Existing artifact migration | `scripts/migrate_artifacts.py` |
| Native topology and centerline geometry | `src/calculations/topology/` |
| Segment sampling, measurement, and profile calculations | `src/calculations/vessel_segments/` |
| Velocity estimator semantics | `src/pipelines/velocity/estimation.py` and `src/pipelines/velocity/semantics.py` |
| Release behavior | `.github/workflows/release.yml`, `build_installer.ps1` |

## Resource and release notes

Cross-section work is temporally chunked. CPU/GPU backend selection lives in
`calculations/compute_backend.py`; vessel-segment sampling enforces a shared
scratch-memory budget, while retained outputs are outside that budget. Preserve
the staged gradient path's temporal halo and operation order when optimizing.

The checked-in GitHub workflow is release-only: a newly created `vX.Y.Z` tag at
the current `dev` head is validated, the tree is squash-published to `main`, an
installer is built, and refs are updated atomically. Do not infer continuous
test coverage from that workflow.
