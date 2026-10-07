# EyeFlow

EyeFlow is a Python analysis engine and desktop application for retinal Doppler
holography. It consumes HoloDoppler image data plus DopplerView segmentation,
runs dependency-aware scientific pipelines, and produces HDF5 results, plots,
movies, and reports. Vessel segmentation and AI inference are intentionally
outside its scope.

## Install and run

EyeFlow requires Python 3.10 or newer.

```powershell
python -m venv .venv
.\.venv\Scripts\activate
pip install -e .
```

Launch the desktop application:

```powershell
eyeflow
```

Run the CLI over a `.holo` file, a text list of `.holo` files, a folder tree,
or a ZIP archive:

```powershell
eyeflow --data path\to\input.holo
eyeflow-cli --data path\to\folder --output path\to\results
```

`eyeflow` opens the GUI when invoked without arguments and uses the CLI when
arguments are supplied. `eyeflow-cli` always uses the CLI. Run
`eyeflow-cli --help` for pipeline-list and ZIP-output options.

## Input runs

A `.holo` file anchors separately generated HD and DV companions. For
`sample.holo`, the standard input locations are:

```text
sample/
  sample_HD/h5/sample_HD_output.h5
  sample_DV/h5/sample_DV.h5
```

Normal runs require both companions. HD volumes use `(frame, y, x)` axes. DV
provides `(y, x)` artery and vein masks under
`/segmentation/Retina/artery_mask` and
`/segmentation/Retina/vein_mask`, plus optional labels and optic-disc geometry.
Timing and pixel-pitch metadata are also required by calibrated products.

See [data contracts](docs/data-contracts.md) for exact paths, fallback names,
sidecar configuration, validation, axes, and units.

## Velocity estimation

Every velocity run attempts both methods:

- `doppler_moments` uses HD `moment0` and `moment2`, producing calibrated
  physical velocity in `mm/s` under `/Processing`.
- `frequency_bands` (band ratio) uses exact HD datasets
  `/band_0_3000_9000` (LF) and `/band_1_9000_18000` (HF), both `(frame, y, x)`.
  It converts `HF / LF` to RMS frequency using the positive finite persisted
  `band_ratio_frequency_scale_hz` setting, then applies the same local-background,
  background-difference, and physical velocity conversion. Outputs use
  `/ProcessingAlt`; the calibration default is `1 Hz` per ratio unit.

Cardiac cycles are detected once before either estimate, from raw moment-derived
RMS frequency, falling back to calibrated HF/LF frequency when moments are
unavailable or unusable. Both workflows use identical cycle timing and shared
`/Segmentation` geometry. Selected downstream products are recalculated for
each method, and PNGs, videos, EPS files, and PDFs use separate `moments/` and
`bandratio/` folders beneath their artifact-type directories.

Missing inputs, invalid data, or a downstream failure skip only the affected
workflow. A run fails if neither workflow completes. The HDF5 root records
completed methods and failure reasons; method-specific provenance is attached
to each processing group and its datasets. Exact-zero LF values map to ratio
zero without an epsilon, and LF quality counts remain diagnostic.

There is no method selector. Legacy `velocity_estimation_method` settings are
ignored. Fresh-install defaults are in `default_settings.json`; velocity is
always physical in `mm/s`, and legacy dimensionless velocity metadata is rejected.


## Pipelines and outputs

The GUI and CLI use the same catalog, dependency resolver, and execution
service. Visible pipelines select products; hidden core pipelines automatically
prepare shared heartbeat, topology, and velocity state. Options such as
per-beat signals, segments, profiles, quadrants, and blood-volume-rate families
add their own dependencies.

For `sample.holo`, EyeFlow owns `sample/sample_EF/`. A new attempt replaces an
existing directory before running. The primary result is
`sample_EF/h5/sample_EF.h5`; artifact directories such as `png`, `avi`, `eps`,
and `pdf` are created when needed. The output HDF5 records source files,
selected targets, resolved execution order/options, version data, and velocity
semantics.

Artifact filenames use the acquisition prefix, for example
`sample_lumen_size_by_branch_artery.png`. Writers preserve subfolders and
avoid adding the prefix twice. To rename artifacts in existing result folders,
run `python -m input_output.artifact_migration path\to\results --dry-run`
to inspect the changes, then omit `--dry-run` to apply them.

Major result families include continuous and per-beat velocity, heartbeat,
topology/segmentation, cross-section profiles, waveform metrics, spatial-gradient
lumen metrics, blood-volume-rate products, displacement products, and reports.
The executable output-path schema is centralized in
`src/input_output/schema/eyeflow_output.py`.

## Development

Run the test suite with:

```powershell
pip install pytest
python -m pytest
```

Repository navigation for contributors and coding agents starts at
[AGENTS.md](AGENTS.md). Deeper references are:

- [architecture and source-of-truth map](docs/architecture.md)
- [data contracts](docs/data-contracts.md)
- [test-to-production map](test/AGENTS.md)
- [pipeline contribution guide](CONTRIBUTING.md)
- [weighted velocity-profile analysis](docs/velocity_profile_analysis.md)

The tag-triggered GitHub workflow builds Windows releases from `dev`; it is not
a continuous test workflow. `benchmarks/rtx4090_cross_section.json` is a
hardware-specific cross-section parity/performance snapshot, not an automated
benchmark harness.
