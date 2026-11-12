# Architecture

EyeFlow starts from a `.holo` run or matching run folder, reads HoloDoppler
(HD) and DopplerView (DV) companions, and writes one `<stem>_EF` result. DV
supplies segmentation; EyeFlow does not segment vessels. See
[data contracts](data-contracts.md) for input and output boundaries.

## Entry points and ownership

`pyproject.toml` exposes `eyeflow` and `eyeflow-cli`. `src/launcher.py` sends
arguments to `src/cli.py` and a no-argument desktop launch to
`src/eye_flow.py`. Both front ends select inputs and targets through
`pipeline_engine.run_service.resolve_run_spec`, then call `execute_run`.
Folder and marker discovery is shared in `input_output/runs/`.

| Layer | Owns | Start with |
|---|---|---|
| Front ends | Input, selection, progress, presentation | `src/cli.py`, `src/ui/` |
| Run service and engine | Validation, DAG order, execution, failures, resource admission | `src/pipeline_engine/run_service.py`, `dag.py`, `runtime.py` |
| Pipelines | Declarations, source adaptation, run state, product packing | `src/pipelines/<name>/` |
| Calculations | Reusable numerical operations | `src/calculations/` |
| I/O and schemas | HD/DV adapters, HDF5 paths, output serialization, artifacts | `src/input_output/` |
| Settings | Defaults, persisted settings, validation | `default_settings.json`, `src/app_settings.py` |

Keep calculations independent of settings, Tkinter, and result directories.
Pipelines may combine calculations and runtime state; front ends must not
implement scheduling or science. External path aliases belong in
`input_output/schema/`, and current EyeFlow paths come from
`EyeFlowOutputPaths.active()`.

## Pipeline execution

`pipelines.load_pipeline_catalog()` discovers top-level pipeline packages.
Their `@pipeline` or `@registerPipeline` declarations define options,
visibility, and dependencies. `requires` names importable packages;
`dag_requires` names produced data. A pipeline implicitly produces its name.
A missing producer key is treated as external, so consumers must validate
external values. The resolved DAG and option closure are shared by GUI and CLI.

`resolve_run_spec` validates settings, target availability, dependencies,
companion layouts, and unique destinations before execution. `execute_run`
replaces the EyeFlow-owned output directory per attempt, runs the plan, and
collects partial failures. Cancellation is checked between files. Runtime
contexts share run-local `ctx.state`; durable products go through
`ctx.output.h5`. Variant producers create isolated velocity workflow state and
HDF5/artifact namespaces. See [parallel execution](parallel_pipeline_execution.md)
for process scheduling and publication, and [stage graphs](pipeline_stage_graphs.md)
for product-level handoffs.

The main scientific flow is:

```text
HD + DV -> velocity (moments and/or calibrated bands) -> cardiac timing
        -> topology_core -----------------------------+-> velocity_analysis
        -> low_freq_spectra (optional HD video)       |       -> metrics/reports
        -> spatial_gradient_moment0 (timing + topology)
        -> displacement_map -> displacement_profiles (optional)
```

`blood_volume_rate` chooses gradient or mask geometry through its declared
options. The pipeline declarations, not this sketch, are the current graph.

## Settings and outputs

`AppSettingsStore` owns persisted settings and import validation; the GUI and
CLI use the same store. `velocity_estimation_method` selects
`doppler_moments`, `bandratio`, or `both` (default). The legacy
`frequency_bands` value is accepted as a setting alias. Calibration passes
through the run spec to the velocity context. The GUI has no method selector.

The coordinator opens HD/DV inputs, stages requested datasets in local read-only
storage, accepts pipeline products, and writes the final HDF5 locally before
publishing it to the destination. Reports read the assembled output view after
required artifacts complete. `input_output` owns the serialized schema and
writers; pipelines own scientific values and report page selection.

## Sources of truth

| Concern | File or package |
|---|---|
| Commands and dependencies | `pyproject.toml` |
| Pipeline metadata and discovery | pipeline `__init__.py` files; `src/pipelines/__init__.py` |
| DAG and run validation | `src/pipeline_engine/dag.py`, `run_service.py` |
| Execution and context | `src/pipeline_engine/runtime.py`, `context.py` |
| HD/DV paths and alignment | `src/input_output/schema/`, `readers/canonical.py` |
| Output HDF5 paths and provenance | `src/input_output/schema/eyeflow_output.py`, `provenance.py` |
| Velocity semantics | `src/pipelines/velocity/estimation.py`, `semantics.py` |
| Segment geometry and profiles | `src/calculations/topology/`, `vessel_segments/` |
| Figure styling | `src/plotting/` |
| Installer and release | `build_installer.ps1`, `installer/`, `.github/workflows/release.yml` |
