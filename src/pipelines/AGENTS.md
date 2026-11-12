# Pipeline map

Pipelines declare their dependencies/options, adapt typed sources to
calculations, own run-local state keys, and pack products. The engine owns DAG
execution; `calculations` owns reusable formulas; `input_output/schema` owns
shared HDF5 paths. Discovery imports top-level packages automatically through
`src/pipelines/__init__.py`. Start with the affected package's `__init__.py`
declaration, then its runner and output module.

| Product or concern | Start with |
|---|---|
| Canonical HD/DV source assembly | `vessel_inputs.py`, `src/retinal_source_data.py` |
| Physical velocity and cardiac timing | `velocity/` |
| Shared topology and gap masks | `topology_core/` |
| Per-beat, segment, profile, map, and movie products | `velocity_analysis/` |
| Shared HD spectra, metrics, and QC | `low_freq_spectra/` |
| Gradient-derived lumen products | `spatial_gradient_moment0/` |
| Flow rate option families | `blood_volume_rate/` |
| Displacement registration and profiles | `displacement_map/` |
| Waveform descriptors | `waveform_shape_metrics/`, `absolute_waveform_metrics/`, `lowrank_waveform_decomposition/` |
| PDF measurements and page composition | `pdf_report/` |

Hidden `velocity`, `topology_core`, and `low_freq_spectra` publish shared work.
Visible products declare their upstream requirements; never infer availability
from UI selection alone. State keys need one owning producer and a declared DAG
edge. Durable values go through `ctx.output.h5` with units, axes, and provenance.

Velocity moments and calibrated frequency bands produce separate physical
`mm/s` workflows under `/ProcessingMoments` and `/ProcessingBandRatio`.
Shared segmentation and spectra are written once. Preserve method,
calibration, source, unit, and timing information in derived outputs. See
[data contracts](../../docs/data-contracts.md) for input/output boundaries;
[spectral metrics](../../docs/spectral_metric_selection.md) and
[profile analysis](../../docs/velocity_profile_analysis.md) for those
specialized numerical contracts.

For new options and examples, see [CONTRIBUTING](../../CONTRIBUTING.md).
Search producer, consumers, and focused tests before changing a shared state
key or output family.
