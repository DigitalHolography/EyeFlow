# Input/output map

This package owns run discovery, HD/DV adapters, HDF5 schemas and writers,
result folders, artifacts, archives, and PDF file I/O. Pipelines own scientific
values and report page selection. Read [data contracts](../../docs/data-contracts.md)
when changing a public path, shape, unit, or layout.

| Task | Start with |
|---|---|
| `.holo` companions and folder/list/ZIP input | `runs/layout.py`, `runs/discovery.py`, `runs/inputs.py` |
| External HD/DV paths and logical fields | `schema/holodoppler.py`, `schema/doppler_view.py` |
| Canonical mask/image alignment | `readers/canonical.py`, `schema/serialization.py` |
| Current EyeFlow output paths | `schema/eyeflow_output.py` |
| Effective attributes and provenance | `schema/provenance.py` |
| HDF5 output and storage policy | `writers/h5_output.py`, `writers/h5.py`, `writers/h5_storage.py` |
| Output folders and artifact names | `output_manager.py`, `writers/artifact_names.py` |
| Figures and reports | `writers/article.py`, `reports/`; scientific layout stays in pipelines and `src/plotting/` |

Use typed adapters for external aliases and `EyeFlowOutputPaths.active()` for
shared output families. `PipelineH5Output` is the pipeline-facing writer.
Native HD/DV arrays remain lazy until bounded reads; canonical alignment moves
related DV masks and optic-disc geometry together. Preserve `dimDesc`, units,
and workflow provenance. Read effective attributes through the provenance
helper because raw HDF5 dataset attributes may omit inherited values.

Search the producer, every consumer, and a focused schema or writer test before
changing a path or representation. `execute_run` owns replacement of an
existing result directory; writers and pipelines do not.
