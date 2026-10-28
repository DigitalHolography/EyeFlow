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
| `.holo` companion resolution | `runs/layout.py` | `schema/base.py`, `run_service.py` |
| Recursive run-folder discovery and batch expansion | `runs/discovery.py`, `runs/inputs.py` | GUI input controller and `run_service.py` |
| External frame and mask readers | `readers/frame_sequences.py`, `readers/canonical.py` | `pipelines/displacement_map/`, `pipelines/vessel_inputs.py`, and focused tests |
| Raw access and merged attributes | `readers/h5.py` | the typed source adapter; context conveniences live in `pipeline_engine/input_access.py` |
| HD/DV datasets and sidecars | `schema/holodoppler.py`, `schema/doppler_view.py`, `schema/base.py` | `src/retinal_source_data.py`, `readers/canonical.py` |
| Current output paths | `schema/eyeflow_output.py` | producer and consumer pipelines |
| HDF5 writes | `writers/h5_output.py`, `writers/h5.py`, `writers/h5_storage.py`, `writers/h5_initialization.py`, `schema/output_payload.py` | output-schema tests and producer tests |
| Output folders/artifacts | `output_manager.py`, `writers/artifact_names.py` | the specific writer or report module; `scripts/migrate_artifacts.py` at the repository root for existing results |
| Article figure styling | `src/plotting/article.py` | the producer's figure builder; writers save already formatted figures |
| General ticks/legend/colorbar layout | `src/plotting/formatting.py` | article styles or non-article producer; no font/export policy is imposed |
| Article PNG/EPS pairs | `writers/article.py` | `src/plotting/article.py` and producer artifact tests; general PNG exports stay opt-in |
| ZIP inputs | `archives/zip_archive.py` | input-expansion tests in `test_run_service.py` |
| Profile dataset serialization | `writers/h5_storage.py` | profile producer/consumer tests |
| Displacement frame and artifact I/O | `readers/frame_sequences.py`, `writers/displacement_arrays.py`, `writers/magnitude_avi.py` | `pipelines/displacement_map/` and its focused test |
| JSON and text sidecars | `writers/sidecars.py` | displacement metadata and run logging |
| PDF measurements and page layout | `reports/pdf_data.py`, `reports/pdf_report.py` | report image lookup and runner tests |

## Contracts and ownership

Typed source adapters translate stable logical operations into external HDF5
paths. Pipeline code should use these adapters instead of repeating path
fallbacks. `readers/canonical.py` owns canonical DV-to-HD alignment, moving
masks and optic-disc geometry atomically. Pipeline-specific fallback policy
belongs to `pipelines/vessel_inputs.py`. `schema/serialization.py` owns
published output axis order.

`runs/discovery.py` shares marker and marker-free run discovery between GUI and
CLI. A run folder must contain its own matching HD and DV directories; EF/AE
are optional. Discovery checks structure, while normal layout resolution checks
the HDF5 files. Keep actual `HoloRunLayout.input_path` distinct from its canonical
`holo_path` naming anchor; never create marker files during discovery.
`runs/inputs.py` validates selected markers, folders, and text lists, expands
ZIP/folder selections, and computes batch roots. `output_manager.py` maps
validated runs to output destinations and rejects destination collisions.

`EyeFlowOutputPaths.active()` owns reusable output path families. Before
changing a path, search all source and tests for both the field name and current
literal. If a pipeline owns a specialized path constant, keep its producer,
consumers, tests, and documentation synchronized.

`PipelineH5Output` in `writers/h5_output.py` is the pipeline-facing API.
`writers/h5.py` owns generic file and dataset writes; `writers/h5_storage.py`
owns payload conversion, compression, and bool representation;
`writers/h5_initialization.py` owns root initialization and selected HD
pass-through data. Avoid opening the primary HDF5 independently inside a pipeline.
Large float32 map/profile payloads may request sparse NaN chunk storage; the
published shape and values remain dense on read. Symmetric velocity-profile
FFT magnitudes may use a same-file virtual dataset over compressed unique
frequency bins. The writer falls back to ordinary storage when exact symmetry
or a useful size reduction is absent. HD registration and Zernike pass-through
copies are selected by the run-level `copy_hd_passthrough` setting.

Per-beat profile interpolation and dataset packing live in
`pipelines/profile_outputs.py`; they are pipeline adaptation rather than I/O
serialization policy.

`OutputManager` creates type directories lazily. Supported types are H5, PNG,
MP4, AVI, PDF, and EPS; JSON and text sidecars use dedicated writer helpers. A pipeline may
reserve a path or open an auxiliary HDF5, but required machine-readable results
belong in the primary file.
PNG figure writers record report image paths in run state so PDF assembly can
use the actual variant artifacts; legacy filename lookup remains available.

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
- Processing parents retain only the fixed summary attributes defined in
  `schema.provenance.processing_summary_keys()`. Other parent metadata is discarded.
  Use `schema.provenance.read_dataset_attrs()` or `PipelineH5Output.read_attrs()`
  for effective attributes; dataset-local overrides win. Keep workflows separate.
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
