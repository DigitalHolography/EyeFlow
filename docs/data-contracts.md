# Data contracts

This is the authoritative documentation layer for data moving into, through,
and out of EyeFlow. Executable path constants in `src/input_output/schema/`
remain the final source of truth.

## External run layout

A selected `sample.holo` anchors this default tree:

```text
<parent>/sample.holo
<parent>/sample/
  sample_HD/h5/sample_HD_output.h5
  sample_HD/json/parameters.json          # optional sidecar values
  sample_DV/h5/sample_DV.h5
  sample_DV/json/DV_params.json           # optional sidecar values
  sample_EF/                               # EyeFlow-owned output
```

For each HD/DV `h5` directory, the preferred filename wins. If it is absent,
exactly one valid `.h5`/`.hdf5` file is accepted; zero or multiple ambiguous
files fail. Normal production runs require both HD and DV companions. Folder,
text-list, and ZIP inputs are expanded into `.holo` selections before layout
resolution.

## HD input

All image volumes are lazy HDF5 datasets with axes `(frame, y, x)`.

| Logical value | HDF5 path(s) | Requirement |
|---|---|---|
| Zeroth moment | `/moment0`, legacy `/M0` | Required by the moments workflow and moment-based products |
| Second moment | `/moment2`, legacy `/M2` | Required by the moments workflow |
| Low-frequency band | `/band_0_3000_9000` | Required by the band-ratio workflow; no heuristic aliases |
| High-frequency band | `/band_1_9000_18000` | Required by the band-ratio workflow; no heuristic aliases |
| Flat-field zeroth moment | `/moment0ff`, legacy `/M0FF` | Optional |
| Registration | `/registration` | Optional; copied to output `/Meta/registration` |
| Zernike coefficients | `/zernike_coefs_radians` | Optional; copied through unchanged |

Moment and band datasets must be 3-D. The two bands must be numeric and have
identical shapes; their chunk values must be finite and non-negative.

Timing uses scalar `sampling_freq` and `batch_stride`, resolved from HDF5 or HD
sidecar configuration. The frame interval is
`batch_stride / sampling_freq` seconds. Spatial calibration is a scalar JSON
value at `/HD_parameters` containing `pixel_pitch: [x, y]` in metres. Current
profile consumers require the two pitches to be approximately equal before
using the isotropic value.

## DV input

DV spatial arrays are native `(y, x)` values:

| Value | HDF5 path | Requirement |
|---|---|---|
| Artery mask | `/segmentation/Retina/artery_mask` | Required |
| Vein mask | `/segmentation/Retina/vein_mask` | Required |
| Labeled vessels | `/segmentation/Retina/labeled_vessels` | Optional |
| Optic-disc mask | `/segmentation/OpticDisc/mask` | Optional if fallback/other measurements suffice |
| Optic-disc center | `/segmentation/OpticDisc/center` | Optional input; canonical source always resolves a center |
| Optic-disc width/height | `/segmentation/OpticDisc/width`, `/segmentation/OpticDisc/height` | Used to reconstruct an ellipse when no mask is available |

`VelocityEstimation/LocalBackgroundDist` comes from the DV sidecar configuration
and defaults to `2`. `pipelines.vessel_inputs.load_retinal_source_data` is the
only canonical alignment boundary: it validates masks against the active HD
spatial shape and may transpose all DV geometry together when axes are swapped.
Do not add isolated transposes in downstream algorithms.

Missing optic-disc geometry falls back to a centered circle. That fallback
retains the original vein mask as velocity-background support but disables vein
measurement in the canonical segmentation.

## Velocity methods

Every run that schedules velocity attempts both workflows, independently:

- `doppler_moments`: derives RMS frequency from moments, estimates local
  background, subtracts it with the established signed RMS rule, and converts
  frequency to calibrated velocity in `mm/s`. Its outputs use `/Processing`.
- `frequency_bands` (band ratio): computes `HF / LF` and converts it to RMS
  frequency with `fRMS_Hz = band_ratio_frequency_scale_hz * (HF / LF)`, then
  applies the same vessel-mask dilation, biharmonic background inpainting,
  signed background-difference, and frequency-to-velocity conversion. Its
  outputs use `/ProcessingAlt`. The persisted calibration default is `1 Hz`
  per ratio unit; quantity is `physical_velocity`, unit `mm/s`.

`velocity_estimation_method` is no longer a selectable or persisted setting.
Legacy configuration values are ignored. Method identifiers remain in workflow
provenance and in internal method-specific contexts.

Cardiac cycles are detected independently for each workflow **before** velocity
estimation, from raw RMS frequency averaged over artery/vein support: moment-derived
frequency for moments, and calibrated HF/LF frequency for bands. Each workflow
uses the existing artery, then vein, then full-record fallback. Its own cycle
boundaries, durations, and spectral timing are persisted under its processing
group's `CardiacCycle` family and used throughout its downstream analyses.
Velocity filtering always filters each workflow's physical velocity signal
independently of the detector's signal.
The workflows must have identical `(frame, y, x)` dimensions.

All selected downstream analyses are recalculated in isolated run state for
each successful estimate, using shared segmentation/topology and workflow-specific
cycle boundaries. A failed input, cycle detection, estimate, or downstream workflow is skipped; its
processing group and method artifacts are removed, while the other workflow
continues. If both fail, the run fails. Root `velocity_workflow_failures`
records a JSON map of method to error; `velocity_workflows_completed` lists
successful methods when the run finishes.

For band mode, an exactly zero LF sample maps to ratio zero. No epsilon is
added. A nonzero ratio beyond finite `float32` range raises a clear error rather
than emitting infinity. Missing exact band paths, mismatched shapes, nonnumeric
data, NaN/Inf, and negative power values also fail explicitly. A failed band
workflow is skipped independently of the moments workflow.

Band outputs report the LF quality threshold (`1e-6` of each frame maximum) and
exact-zero/near-zero sample counts for the full volume, vessel pixels, and the
unmasked neighborhood used by inpainting. These counts are diagnostic only and
do not change the zero rule or discard low samples.

All pipeline targets remain schedulable in band mode, including
`blood_volume_rate` and `absolute_waveform_metrics`. Consumers must inspect the
processing-group method, calibration, quantity, and unit provenance when
interpreting derived values.

## Internal contracts

`RetinalSourceData` in `schema/source_data.py` is the shared typed input:

- lazy `ImageMaps` for moments or bands;
- aligned boolean artery/vein/background masks and optional labels;
- resolved optic-disc geometry;
- timing, pixel pitch, background distance, alignment flag, and velocity method.

Common axes used by pipeline products include:

| Shape | Meaning |
|---|---|
| `(frame, y, x)` | HD images and velocity video |
| `(y, x)` | masks and spatial maps |
| `(time, beat)` | global per-beat waveform |
| `(time, beat, branch, radius)` | segment per-beat signal |
| `(x, time, beat, branch, radius)` | transverse velocity profile before profile fitting |

Specific datasets carry `dimDesc`; treat it as part of the contract. Beat and
frame indexes are zero-based unless a dataset explicitly says otherwise.

Run-local `ctx.state` contains objects that should not be serialized merely to
connect pipelines. The producer must declare a DAG key and the consumer must
retrieve the expected typed state. Stable or externally consumed results belong
in the output HDF5.

## Output ownership and schema

The output root is `<stem>_EF`. `execute_run` removes an existing output root
before a new attempt, so a failed run leaves only that attempt's partial files.
Subdirectories are created lazily for `h5`, `png`, `mp4`, `avi`, `pdf`, and
`eps`. Velocity-derived artifacts are separated beneath each type directory
into `moments/` and `bandratio/`, including PDF reports. Each report reads its
corresponding processing group. The primary file is `h5/<stem>_EF.h5`.

All artifact basenames start with `<stem>_`, including PNG, EPS, AVI, MP4,
PDF, auxiliary HDF5, and archive outputs. The shared artifact naming helper
and `OutputManager.path_for()` enforce this at the writer boundary: existing
prefixes are retained once and subdirectories are preserved. The primary
`<stem>_EF.h5` filename is unchanged.

Existing EyeFlow result folders can be migrated with
`python -m input_output.artifact_migration RESULTS_PARENT --dry-run`, then
the same command without `--dry-run`. The utility only renames files inside
identified `<stem>_EF` folders, leaves directories and contents unchanged,
and checks for collisions before renaming.

`EyeFlowOutputPaths.active()` defines the current `eyeflow_v2` paths. Important
families are:

- `/Processing/Velocity/global`, `/segments`, and `/VelocityPerBeat`;
- `/Processing/Heartbeat` and `/Processing/FrequencyMaps`;
- `/Processing/VelocityProfiles` and `/VelocityProfilesFFT`;
- `/Processing/Metrics/{waveform_shape_metrics,absolute_waveform_metrics,lowrank_waveform_decomposition}`;
- `/Processing/SpatialGradientMetrics` for a directly selected spatial-gradient
  target, and `/Processing/BloodVolumeRate` for blood-volume-rate products;
- `/Segmentation` for aligned masks, topology, areas, and lumen geometry;
- `/Meta` for provenance and selected pass-through data.

The `/ProcessingAlt` tree mirrors the selected `/Processing` families for the
band-ratio workflow. `/Segmentation` and acquisition `/Meta` remain shared.

Per-beat durations are published once per workflow at
`/Processing/CardiacCycle/Systole/CycleDurationSeconds/value` as a `(beat,)`
vector. The
`/Processing/VelocityPerBeat` family contains beat-aligned velocity waveforms,
not a duplicate duration dataset.

`/Processing/Maps/VelocityAverageMasked/value` is the temporal mean of the
background-subtracted velocity used by the vessel waveform analysis. The
unmasked velocity average and intermediate delta-frequency average are not
persisted.
All two-dimensional datasets under `/Processing/Maps` are serialized in the
same lower-left `(x, y)` image frame as
`/Segmentation/Artery/BranchLabelMap/value`: the internal `(y, x)` array is
flipped vertically and transposed before persistence.

Do not hand-copy path strings into new consumers. Obtain them from the schema or
from a pipeline-owned constant when the path family has not yet been centralized.
Writers normalize leading/trailing separators, replace an existing dataset at
the target path, downcast float/complex payloads to 32-bit, preserve integer
arrays, serialize boolean payloads as `uint8` with `original_class="bool"`, and
attach `nameID` unless supplied.

The output root records source files, selected targets, actual pipeline order,
selected options, requested methods, and workflow status. Processing groups and
their datasets record method-specific physical velocity provenance. The
band-ratio group additionally records the
exact LF/HF source paths, factor, linear-through-origin calibration model,
calibration source/version, wavelength, numerical aperture, and LF quality
counts. Persisted RMS-frequency maps use `Hz`; velocity datasets use `mm/s` and
carry method/calibration provenance. Legacy dimensionless velocity metadata is
not accepted.

When spatial-gradient processing runs only to supply the `gradient_edges`
blood-volume-rate option, its profiles and metrics remain transient run state:
`/Processing/SpatialGradientProfiles` and
`/Processing/SpatialGradientMetrics` are not written. Selecting
`spatial_gradient_moment0` directly publishes both families. Gradient-derived
blood-volume-rate datasets identify these inputs as `transient_run_state` and
do not expose dangling HDF5 source-path attributes.

The release-default selection disables `gradient_edges`, velocity profiles,
and the waveform `velocity_profile_analysis` option. Its HDF5 therefore omits
`/Processing/SpatialGradientProfiles`,
`/Processing/SpatialGradientMetrics`,
`/Processing/BloodVolumeRate/{Artery,Vein}/{dynamicEdges,staticEdges}`,
`/Processing/VelocityProfiles`, and
`/Processing/VelocityProfileAnalysis`. Mask-derived `maskedEdges` and
`totalMaskedEdges` blood-volume-rate datasets remain selected.

For the weighted profile-analysis schema, see
[velocity-profile analysis](velocity_profile_analysis.md).
