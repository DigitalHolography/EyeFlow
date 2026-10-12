# Input/output context

## Boundary

This package owns external run discovery, typed HD/DV access, output-folder and
HDF5 serialization policy, archives, media writers, and PDF assembly. It does
not decide which pipelines run or implement scientific formulas.

Read [data contracts](../../docs/data-contracts.md) before changing a path,
shape, unit, layout, or writer behavior.

## Navigation

| Area | Start with | Follow only if needed |
|---|---|---|
| `.holo` companion resolution | `holo_run_layout.py` | `schema/base.py`, `run_service.py` |
| Raw access and merged attributes | `h5_access.py`, `input_access.py`, `inputs.py` | the typed source adapter |
| HD/DV datasets and sidecars | `schema/holodoppler.py`, `schema/doppler_view.py` | `schema/source_data.py`, pipeline source loader |
| Current output paths | `schema/eyeflow_output.py` | producer and consumer pipelines |
| HDF5 writes | `h5_access.py`, `writers/h5.py` | output-schema tests and producer tests |
| Output folders/artifacts | `output_manager.py`, `writers/artifact_names.py` | the specific writer or report module; `artifact_migration.py` for existing results |
| ZIP inputs | `archives/zip_archive.py` | input-expansion tests in `test_run_service.py` |
| Profile dataset serialization | `writers/h5.py` | profile producer/consumer tests |

## Contracts and ownership

Typed source adapters translate stable logical operations into external HDF5
paths. Pipeline code should use these adapters instead of repeating path
fallbacks. Canonical DV-to-HD alignment happens in
`pipelines/vessel_inputs.py`, because it must move masks and optic-disc geometry
atomically.

`EyeFlowOutputPaths.active()` owns reusable output path families. Before
changing a path, search all source and tests for both the field name and current
literal. If a pipeline owns a specialized path constant, keep its producer,
consumers, tests, and documentation synchronized.

`PipelineH5Output` is the pipeline-facing API. `writers/h5.py` owns actual
serialization, downcasting, bool representation, attributes, initialization,
and selected HD pass-through data. Avoid opening the primary HDF5 independently
inside a pipeline.

Per-beat profile interpolation and dataset packing live in
`pipelines/profile_outputs.py`; they are pipeline adaptation rather than I/O
serialization policy.

`OutputManager` creates type directories lazily. Supported types are H5, PNG,
MP4, AVI, PDF, and EPS; there is no generic JSON output helper. A pipeline may
reserve a path or open an auxiliary HDF5, but required machine-readable results
belong in the primary file.

Artifact writers prefix each basename with the acquisition `<stem>_` exactly
once, preserving subfolders. `OutputManager.path_for()` applies the same rule
to reserved paths and auxiliary outputs. The primary `<stem>_EF.h5` name is
unchanged. Keep the naming policy in the shared writer helper.

## Invariants

- Input volumes stay lazy until an algorithm requests a bounded slice.
- A schema adapter must reject malformed rank/type/shape close to the boundary
  and name the missing external path in its error.
- HDF5 logical paths use `/` regardless of host path separators.
- Dataset axes and units travel in attributes such as `dimDesc` and `unit`.
- Existing target datasets are replaced atomically at the dataset level; parent
  groups are reused only if they are groups.
- An output-directory replacement is owned by `execute_run`, not by individual
  pipelines or writers.
- Sidecar artifacts must not be the sole copy of a required analysis result.

## Tests

Use `test_dopplerview_compat.py` for DV/source compatibility,
`test_frequency_band_velocity.py` for exact HD band contracts,
`test_h5_writer.py` for serialization, `test_avi_writer.py` and
`test_eps_writer.py` for artifacts, `test_stem_list_inputs.py` and
`test_run_service.py` for input/layout behavior, and the PDF tests for report
path lookup. Producer-specific output schemas are generally asserted in that
pipeline's test file rather than one global schema test.
