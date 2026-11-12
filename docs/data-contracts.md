# Data contracts

This page describes the stable boundaries between HD, DV, and EyeFlow. For an
exact path or attribute, use the executable schemas in
`src/input_output/schema/`; for a numerical definition, use its calculation or
pipeline implementation and focused tests. Do not copy path strings into new
consumers.

## Run layout

A `sample.holo` marker normally identifies these companions and output:

```text
<parent>/sample.holo
<parent>/sample/sample_HD/h5/sample_HD_output.h5
<parent>/sample/sample_DV/h5/sample_DV.h5
<parent>/sample/sample_EF/h5/sample_EF.h5
```

Optional HD and DV JSON sidecars live beside their `h5` directories. The
preferred HDF5 filename wins; if absent, exactly one valid `.h5` or `.hdf5`
file is accepted. Normal runs require both companions. Folder discovery also
accepts a matching `sample/` directory without a marker, including the selected
directory itself. It never creates a `.holo` file. Marker and folder references
to one run are deduplicated. Text lists, recursive folders, and ZIP inputs are
expanded before shared run validation. See `src/input_output/runs/` for the
exact discovery and list syntax.

`HoloRunLayout.input_path` is the actual selection shown in progress and errors;
`holo_path` is the naming anchor even for marker-free runs. Validation rejects
missing or ambiguous companions before `execute_run` replaces an existing
EyeFlow output directory.

## HD source

HD moments and frequency-band videos are lazy numeric `(frame, y, x)` datasets.
The external aliases belong to `schema/holodoppler.py`:

| Source | HDF5 path | Use |
|---|---|---|
| Zeroth and second moments | `/moment0`, `/moment2` (legacy `/M0`, `/M2`) | Moments velocity |
| Flat-field zeroth moment | `/moment0ff` (legacy `/M0FF`) | Optional gap and spatial-gradient analyses |
| LF and HF bands | `/band_<i>_<n>_<m>` | Band-ratio velocity; first and last numeric ranges |
| Spectral video | `/spectrograms` | Optional shared spectral analysis |
| Registration and Zernike metadata | `/registration`, `/zernike_coefs_radians` | Optional pass-through |

Band videos must have the same 3-D shape and finite numeric samples. The
velocity estimator replaces finite negative LF/HF samples with zero in each
working chunk and warns; source data is unchanged. Zero LF gives ratio zero.

A usable `/spectrograms` is a real numeric 4-D dataset with exported axes
`(time, freq, y, x)`. Older HD exports may have a group of that name; the group
means the optional video is unavailable. A present dataset with invalid rank or
dtype fails spectral validation. Spectral values may be negative after HD
registration and are not clamped at input. Spectral spatial dimensions may
differ from moment and DV dimensions. Gap support is resized to the spectral
`(y, x)` grid with nearest-neighbor pixel-center sampling. Spectral reductions
apply this support in native HD coordinates; published segmentation uses the
EyeFlow orientation below.

Frequency calibration uses HD `batch_size` as the FFT length, `sampling_freq`,
and HD's selected frequency range and binning. The signed bin table remains in
export order, positive bins followed by negative bins. If selection bounds are
missing, use the full FFT range `(start=0, end=0.5)`; do not infer a crop from
output bin count. `calculations/math/frequency.py` owns the conversion. Invalid
calibration omits calibrated tables, metrics, and frequency figures while
retaining available power spectra.

Frame interval is `batch_stride / sampling_freq` seconds. The HD `pixel_pitch`
JSON value is `[x, y]` in metres; consumers requiring isotropic pitch validate
it before conversion.

## DV source and spatial alignment

DV arrays are native `(y, x)`. `readers/canonical.py` aligns all masks and
optic-disc geometry to the HD frame together; downstream code must not add
independent transposes. Required artery and vein masks are
`/segmentation/Retina/artery_mask` and `/segmentation/Retina/vein_mask`.
Optional sources include retinal logits, choroid vessel mask, labeled vessels,
and optic-disc mask, center, width, and height. The retinal logits use `> 0.1`
for vessel support; choroid mask uses nonzero. Missing or malformed optional
sources have different policies: absence may skip a dependent product, while a
present malformed source fails validation.

When optic-disc geometry cannot be used, the canonical source creates a circle
centered on the image with physical diameter 1.91 mm. This requires calibrated
pixel pitch and disables vein measurement while retaining vein background
support. `pipelines/vessel_inputs.py` owns that policy.

Calculations use native `(y, x)` masks alongside HD images. Output
segmentation and 2-D maps flip Y once and transpose to EyeFlow's lower-left
`(x, y)` frame. The output `dimDesc` records the published axes. Do not infer
alignment from equal or swapped dimensions alone.

## Gap support and spectra

`topology_core` can produce three gap masks:

| Output under `/Segmentation/Gaps` | Meaning |
|---|---|
| `Legacy` | Retinal/choroid vessel and optic-disc exclusion within the central cutoff |
| `GapsFromM0ff` | White gaps after flat-field temporal median, contrast enhancement, bright-vessel threshold, and dilation |
| `Mask` | Intersection of available gap supports, with central cutoff and optic-disc exclusion |

`Mask` is the downstream gap support. Missing `moment0ff` uses `Legacy`; missing
retinal or choroid support may use `GapsFromM0ff`. If both are unavailable,
`Mask` is omitted. `GapsFromM0ff` itself remains full-field. All published masks
are `uint8` 0/1 with `original_class="bool"`, `unit="1"`, and `(x, y)` axes.
The numerical operations and fallback parameters live in
`calculations/topology/` and `pipelines/topology_core/`.

The hidden `low_freq_spectra` pipeline runs once after scheduled velocity work
when `/spectrograms` is available. Shared products use
`/ProcessingSpectrograms`, without workflow copies:

| Product | Axes and meaning |
|---|---|
| `Global/SpatialAverage`, `Global/SpatialMedian` | `(time, freq)` whole-frame reductions |
| `FrequencyBinTable` | `(freq,)` signed calibrated centers in kHz |
| `Gaps/SpatialMedian` | `(time, freq)` median over included gap pixels |
| `Gaps/PerBeat/Raw`, `BandLimited` | `(phase, freq, beat)` aligned and optionally Fourier-filtered spectra |
| `Gaps/PerBeat/SVDContinuous` | `(beat * phase, freq)` concatenated filtered beats |
| `Gaps/TimeAggregated/MedianPhaseBeat` | `(freq,)` median over valid phase and beat samples |
| `Gaps/TimeAggregated/MedianBeat` | `(phase, freq)` median across beats |
| `Gaps/TimeAggregated/MedianPhaseBeatEarly`, `MedianPhaseBeatLate` | `(freq,)` phase-window medians of `MedianBeat` |
| `Gaps/TimeAggregated/StdBeat`, `StdBeatEarly`, `StdBeatLate` | `(freq,)` population STD across per-beat phase medians |
| `Gaps/Metrics`, `Gaps/QC` | Calibrated descriptors and descriptive quality measures |

A nonfinite included pixel makes its spatial reduction NaN; excluded pixels do
not matter. Negative finite spatial reductions are clamped to zero afterward.
An empty gap mask skips gap products while whole-frame products remain. Detected
moments cardiac timing is preferred, then detected band-ratio timing. Without
detected timing, spatial products remain and beat products are omitted. Beats
are resampled to 128 phases in `[0, 1)` without interpolation across NaNs.
Early and late windows are closed `[0, 1/3]` and `[1/3, 1]`; a boundary sample
can belong to both. Aggregations ignore NaNs and retain NaN where no valid
samples exist. `BandLimited` uses signed reconstruction and the configured
harmonic count; it differs from velocity's magnitude reconstruction.

Metric definitions, their selection evidence, invalid-value rules, and units
are in [spectral metric selection](spectral_metric_selection.md). QC is
descriptive: it does not reject beats or change metrics. The source pipeline's
`metrics.py`, `qc.py`, and focused tests own exact dataset lists and metadata.

## Velocity workflows

The persisted `velocity_estimation_method` setting selects `doppler_moments`,
`bandratio`, or `both` (default). `frequency_bands` is an accepted legacy setting
alias and the internal band workflow identifier. Only selected workflows load
inputs and run downstream products. The two methods execute independently:

- Moments derive RMS frequency as `sqrt(M2 / spatial_mean(M0))`, flooring the
  denominator at float32 epsilon and preserving NaNs.
- Bands convert `HF / LF` to RMS Hz using the explicit
  `band_ratio_frequency_scale_hz` calibration. Zero LF maps to zero; finite
  negative inputs are clamped only in the estimator.

Both use the shared local-background difference and conversion to **physical
velocity in mm/s**. Each successful method detects its own cardiac cycles from
its physical velocity, then runs downstream analyses with isolated state and
`/ProcessingMoments` or `/ProcessingBandRatio` output. The workflows share
segmentation and spectral results. Preserve method, calibration, quantity,
unit, and timing provenance. Legacy dimensionless velocity metadata is invalid.
A failed method does not erase the other's accepted products; run status and
failure attributes identify partial results.

## Output and serialization

`EyeFlowOutputPaths.active()` in `schema/eyeflow_output.py` owns current HDF5
paths. The primary file is `h5/<stem>_EF.h5` under `<stem>_EF/`. Major roots are
`/ProcessingMoments`, `/ProcessingBandRatio`, `/ProcessingSpectrograms`,
`/ProcessingWidth`, `/Segmentation`, and `/Meta`. Every artifact basename starts
with `<stem>_`; the shared naming helper applies this once. Type directories
(`png`, `eps`, `avi`, `mp4`, `pdf`) are created when used. Reports and figures
are selected by their producing pipelines; use their code for current page
layout and artifact lists.

Common array contracts are `(frame, y, x)` for HD/velocity, `(y, x)` for
native masks, `(x, y)` for published 2-D maps, `(time, beat, branch, radius)`
for per-beat segments, and `(x, time, beat, branch, radius)` for transverse
profiles before fitting. Beat and frame indexes are zero-based unless the
dataset says otherwise. Each dataset's `dimDesc`, `unit`, and provenance are
part of the contract. [Velocity-profile analysis](velocity_profile_analysis.md)
documents the weighted quadratic fit and its outputs.

`PipelineH5Output` is the pipeline-facing HDF5 API. Writers replace target
datasets, downcast floating and complex payloads to 32-bit, preserve integers,
and store booleans as `uint8`. Numeric arrays of at least 4 KiB use gzip level
9 with shuffle by default; explicit storage options override that choice.
Some profile/map outputs use sparse NaN chunks or virtual datasets without
changing their published dense values. `copy_hd_passthrough` controls optional
HD registration/Zernike copies; the default is false.

The root records sources, selections, order, and status. Workflow groups carry
method-specific summary attributes; dataset-specific axes, units, and
calibration remain local or inherit through the schema provenance policy.
Consumers should use `schema.provenance.read_dataset_attrs()` or
`PipelineH5Output.read_attrs()` for effective attributes. Raw HDF5
`dataset.attrs` may omit inherited values. Run-local `ctx.state` carries
transient products only; stable results belong in HDF5 and state consumers
require declared DAG dependencies.
