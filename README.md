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

The persisted `velocity_estimation_method` setting has two values:

- `doppler_moments` is the default. It uses HD `moment0` and `moment2` and
  produces calibrated physical velocity in `mm/s`.
- `frequency_bands` uses exact HD datasets `/band_0_3000_9000` (LF) and
  `/band_1_9000_18000` (HF), both `(frame, y, x)`. It begins with `HF / LF`,
  converts the ratio to RMS frequency using a provisional `1 Hz` per ratio-unit
  calibration, and then applies the same mask-based local-background,
  background-difference, and physical velocity conversion to produce `mm/s`.

Band mode never falls back to moments. Missing bands, invalid values, or
mismatched shapes produce explicit errors. Exact-zero LF values map to ratio
zero without an epsilon. All pipelines remain selectable in band mode,
including `blood_volume_rate` and `absolute_waveform_metrics`; their outputs
carry the stored velocity method, calibration, quantity, and unit provenance.

The fresh-install default is in `default_settings.json`. Existing user settings
are loaded and normalized by `AppSettingsStore`.

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
