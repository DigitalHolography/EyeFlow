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
text-list, and ZIP inputs are expanded into run selections before layout resolution.

Folder discovery accepts the selected directory itself or any descendant whose
name is `<stem>` and which contains matching `<stem>_HD/` and `<stem>_DV/`
directories. `<stem>_EF/` and `<stem>_AE/` are optional; EF/AE alone, a missing
HD/DV directory, or mismatched prefixes do not qualify. This permits marker-free
input via Add folder, drag-and-drop, and CLI folder inputs. Discovery checks the
directory structure; run validation still rejects missing or ambiguous HD/DV
HDF5 files before output cleanup. No `.holo` file is created. Existing real
markers take precedence when discovered alongside their data folders.

`HoloRunLayout.input_path` identifies the actual input (marker or run folder)
for progress and errors. `holo_path` remains a canonical sibling naming anchor
for output mapping even when that marker does not exist. Default outputs stay
inside the discovered root at `<stem>/<stem>_EF/`; a selected run directory
does not introduce an extra `<stem>` nesting level. Folder and marker references
to the same run share a discovery/deduplication key. TXT-list syntax is unchanged.

Text-list entries are independent: absolute `.holo` paths may come from different
directories, drives, or UNC shares. Matching surrounding double or single quotes
are removed; blank lines and a UTF-8 BOM are accepted. Relative paths resolve
against the text list's parent directory. Legacy bare stems still resolve there.
Each entry retains its own companion and default output directory. A text list
must be selected by itself. The GUI uses separate native multi-file and folder
dialogs. Both add to the current file selection, allowing files and recursively
expanded folders from different locations to coexist; Clear resets it. Selecting
a text list replaces the current selection, and adding files/folders replaces
a selected list. Drag-and-drop can combine files and folders too. Overlapping
paths are deduplicated before shared run validation.

## HD input

All image volumes are lazy HDF5 datasets. Moments and LF/HF bands have axes
`(frame, y, x)`; the frequency-resolved capillary source has the separate layout
described below.

| Logical value | HDF5 path(s) | Requirement |
|---|---|---|
| Zeroth moment | `/moment0`, legacy `/M0` | Required by the moments workflow and moment-based products |
| Second moment | `/moment2`, legacy `/M2` | Required by the moments workflow |
| Low-frequency band | `/band_0_3000_9000` | Required by the band-ratio workflow; no heuristic aliases |
| High-frequency band | `/band_1_9000_18000` | Required by the band-ratio workflow; no heuristic aliases |
| Frequency-resolved capillary source | `/spectrograms` | Optional gap-spectrum preprocessing source; exact HD root dataset name |
| Flat-field zeroth moment | `/moment0ff`, legacy `/M0FF` | Optional |
| Registration | `/registration` | Optional; copied to output `/Meta/registration` |
| Zernike coefficients | `/zernike_coefs_radians` | Optional; copied through unchanged |

Moment and band datasets must be 3-D. The two bands must be numeric and have
identical shapes; their chunk values must be finite and non-negative.

`/spectrograms` contains quantitative band-power values and must be a real
numeric 4-D dataset with exported axes `(time, freq, height, width)`. The spelling
`spectrograms` follows the actual HD export; the earlier supplied spelling
`spectograms` was incorrect. It replaces the previous stand-in `band_2_0_6000`
source. `ImageMaps.spectograms` retains its existing internal field name and
can hold the lazy HDF5 dataset for direct source-adapter callers. Runtime
velocity/topology source assembly omits it; `low_freq_spectra` loads it
independently. Legacy HD exports may instead have a `/spectrograms` group containing
`longtimes`, `singlebeat`, and `spectral_cube_parameters`; that group is not the
4-D video and is treated as an unavailable optional source. Existing moments
and band-ratio workflows continue in both cases. Present datasets with invalid
rank or dtype still fail validation. Source assembly
checks rank and dtype without reading the full volume. Spatial Y is dimension 2
and X is dimension 3, with Y inverted relative to EF's published lower-left
`(x, y)` gaps mask. Source assembly preserves the raw ordering and laziness.
`pipelines.vessel_inputs.align_spectrogram_frame` explicitly flips spatial Y,
then swaps the spatial axes on each read frame to yield `(freq, x, y)`.
Never infer this transform from shape or apply the DV mask-alignment flag to
the spectrogram video.

The checked `260907_DEY_6` export (HD commit `340b1b3`) contains
`(640, 32, 256, 256)` float32 spectra while its moment/band volumes and DV masks
use `(640, 512, 512)` and `(512, 512)` respectively. HD preserves the spatial
axis order when resizing spectra; comparison with resized LF maps confirms
the same orientation for this example. `/Segmentation/Gaps/value` and shared
topology state retain the original 512x512 grid. Gap support is interpolated
to the actual last two `/spectrograms` dimensions `(height, width)`, using
nearest-neighbor sampling of pixel centers (`order=0`, `grid_mode=True`).
Arbitrary positive target dimensions and noninteger scale ratios are supported;
matching grids use a copy. The mask stays Boolean and is not averaged or reduced
using all-true blocks. `calculations.topology.downsample_gap_zones` owns this
array transform and receives the target shape explicitly; it does not read HDF5.
`pipelines.low_freq_spectra.masks.gap_velocity_mask` adapts it to the video's spatial
shape. Gap preprocessing now uses this application helper, then flips and
transposes the reduced native mask into the same published `(x, y)` frame as
the aligned spectra. The persisted mask is unchanged.

This HD exporter uses `sampling_freq=37037` Hz and `batch_size=256`, retains
FFT indexes `0:64` followed by `192:256`, and averages each consecutive group
of four bins. The 32 frequency planes are signed and remain in FFT order,
with center spacing 578.703125 Hz within each sign block. Centers run from
217.013671875 to 8897.560546875 Hz, then from -9042.236328125 to -361.689453125 Hz.
No frequency-coordinate dataset is stored. These centers are reconstructed
from the matching exporter and parameters, not universally assumed from the
plane count. They represent averages over four frequency samples, not single
unbinned FFT frequencies. Sampled registered spectra include negative values;
the input adapter permits them and does not clamp or apply the LF/HF bands'
non-negative validation to this source. HD starts from `abs(spectrum_f)**2`,
but its subsequent subpixel Fourier registration of these spatial power maps
can ring and introduce negative undershoots. Linear resizing and frequency-bin
averaging occur afterward and do not enforce nonnegativity.

`pipelines/low_freq_spectra/spectrograms.py` reads each HD spectral frame once and
publishes `/ProcessingSpectrograms/SpatialAverage` and
`/ProcessingSpectrograms/SpatialMedian` as `(time, freq)` arrays. These are
the arithmetic mean and median over every exported spatial pixel, independent
of gap support. Any nonfinite pixel produces NaN for that reduction's
time/frequency pair. Finite negative reductions are clamped to zero after
aggregation. Frequency order and source power units are retained. Outputs are
written after velocity-workflow failure cleanup, so shared spectra survive a
failed moments workflow when the band workflow succeeds. They are produced
once, with no `/ProcessingBandRatio/Spectrograms` duplicate.

`pipelines/low_freq_spectra/preprocessing.py` prepares optional gap support for
the same streaming reduction pass, which computes the arithmetic spatial mean
of the exported band-power samples included by the application mask, for every
time/frequency pair. Only included pixels enter the denominator; excluded
pixels are not zero-filled. The `low_freq_spectra` runner invokes this once using
native gap support from `topology_core`, and publishes `/ProcessingSpectrograms/Gaps/Raw`, with
`dimDesc=["time", "freq"]`, preserving HD frequency-bin order. Finite negative
spatial means are replaced with zero after averaging, before beat subdivision;
individual pixel samples are not clipped before averaging.
Units are copied from the source dataset when supplied,
otherwise `a.u.` with uncalibrated-power provenance. Float64 accumulation is
serialized as float32 by the standard writer. Any nonfinite included sample
produces NaN for its time/frequency pair; nonfinite excluded pixels have no
effect. If native or interpolated gap support is empty, log a
warning and skip the entire gap-analysis stage without reading spectral
gap samples or writing `/ProcessingSpectrograms/Gaps/Raw`. Full-frame spatial
outputs, segmentation, and other pipeline
work remain available. An unavailable/legacy spectra source omits this
product with a warning. Physical capillary velocity remains deferred.

The hidden `low_freq_spectra` pipeline owns all spectral stages. Its
`auto_run_after="velocity"` declaration automatically adds it to any plan that
already schedules velocity, without making it a selectable target or pulling
velocity into unrelated runs. It declares `velocity_workflows` and
`capillary_analysis_zones` dependencies, owns the `gap_spectrograms` and
`gap_per_beat_spectrograms` products, and executes once in the shared output
namespace (`run_once=True`). Missing or legacy-group `/spectrograms` skips all
spectral processing. Invalid datasets fail this pipeline without disabling
velocity/topology source assembly or independent downstream velocity analyses.
Velocity publishes typed successful workflow results for shared consumers.
Detected moments cycles take priority; unavailable moments results, fallback-only
moments timing, or fewer than two moments boundaries cause a retry with detected
band-ratio cycles. If neither workflow has detected intervals, spatial/raw gap
spectra remain available and per-beat products are omitted with a warning.

`/ProcessingSpectrograms/Gaps/PerBeat/Raw` contains `(phase, freq, beat)` spatial mean band
power using consecutive zero-based cycle boundaries from the selected workflow.
Each beat is linearly resampled to 128 phases in `[0, 1)`. The next boundary is included as
an interpolation anchor, but its duplicated phase endpoint is not stored.
Intervals outside the detected boundary sequence are omitted. Invalid samples
remain NaN, with no interpolation across NaN gaps. Frequency-bin order and
power units are preserved. Attributes identify the selected timing workflow and
its HDF5 source path, boundary indexes, original beat durations, interpolation,
and phase samples.
Spectra and the selected timing signal must have identical frame counts. The
detector's full-record fallback (`cardiac_cycle_source="none"`) is never used as
detected timing. Missing detected timing in both workflows skips per-beat gap
analysis with a warning, preserving shared preprocessing.

`/ProcessingSpectrograms/Gaps/TimeAggregated/Median` contains a float32 `(freq,)`
pooled median of `/ProcessingSpectrograms/Gaps/PerBeat/Raw` over its phase and
beat axes. Ignore NaNs; a frequency with no valid samples remains NaN without
an all-NaN warning. Every valid phase/beat sample has equal weight, with no beat
duration weighting and no median-of-beat-medians reduction. Preserve frequency
order and power units. Attributes record the source path, aggregation axes,
weighting, source dimensions, and valid sample counts per frequency. Produce
this shared dataset whenever raw per-beat gaps are available, with no
`ProcessingBandRatio` duplicate or fallback to continuous/filtered gap spectra.

When the HD spectral video is available, `low_freq_spectra` also exports each
available gap product to both `png/spectrograms/` and `eps/spectrograms/`:
`<stem>_spectrograms_gaps` displays `/ProcessingSpectrograms/Gaps/Raw` as a
linear viridis heatmap with horizontal frame tick labels every 50 samples, every
fourth frequency-bin label, and a vertical power colorbar;
`<stem>_spectrogram_gaps_median` plots
`/ProcessingSpectrograms/Gaps/TimeAggregated/Median` against integer bins in
the original HD positive-then-negative order, labeling every second bin vertically. Both use
an 8 cm wide, 1.618:1 canvas, 8 pt Times New Roman axis captions, 5 pt tick labels
with horizontal Y-axis tick labels, 0.4 pt black borders, and 320 DPI.
`<stem>_gaps_mask` displays the native `/Segmentation/Gaps/value` mask as white
gaps on black in its published lower-left spatial frame, on an 8 cm square
canvas without ticks, labels, borders, or surrounding white margins.
Missing products omit only their own figures;
the median figure requires detected cardiac timing, while the mask does not
require nonempty gap support. Figures do not change any persisted data.

`/ProcessingSpectrograms/Gaps/TimeAggregated/Metrics` analyzes that median
spectrum in stored bin order and writes fourteen scalar datasets directly,
without `/value` suffixes: `index_max1`, `max1`, `index_max2`, `max2`,
`index_min1`, `min1`, `index_min2`, `min2`, `h0`, `h1`, `rho_h`, `A0`, `A1`, and `rho_A`.
Peak 1 is the maximum over zero-based bins 0–15; peak 2 is the maximum over
16–31. Include endpoints, ignore NaNs for peak selection, and choose the
lowest index on ties. An all-NaN half has a NaN peak index and value. Valley
indexes are fixed at 15 and 16; valley values are read directly from those bins.
`h0=(max1+max2)/2`, `h1=(min1+min2)/2`, and `rho_h=h1/h0`, with NaN for a zero
denominator. Means propagate missing inputs. `A0=h0*(index_max2-index_max1+1)`;
`A1=sum(spectrum[index_max1:index_max2+1])`, propagating any NaN in that inclusive
interval. `rho_A=A1/A0`, with NaN for a zero or missing denominator; a missing
`A1` also propagates NaN. These are unweighted bin sums without wraparound or Hz-spacing
weights. All values are float32; valid index values are integer-valued floats
so missing indexes can be represented by NaN. Indexes, `rho_h`, and `rho_A` have unit `1`;
other values retain the spectrum's power unit. Attributes record definitions,
source, zero-based indexes, selection rules, and interval inclusion. A spectrum
whose shape is not exactly `(32,)` produces all fourteen NaN metrics and logs a
warning; the median and other gap products remain available. No alternative
workflow copy is generated.

`/ProcessingSpectrograms/Gaps/PerBeat/BandLimited` independently filters the
interpolated raw beats with the shared `calculations.math.bandlimiting.band_limit`
API in `signed` reconstruction mode. Retain DC plus the configured number of
non-DC harmonics, default 12, capped at Nyquist. For a beat of duration
`T`, phase sample spacing is `T / 128` seconds; retaining harmonic `H` gives
an effective cutoff of `H / T` Hz. The count does not depend on duration. Each
frequency bin is filtered separately within each periodic beat. Incomplete
phase series remain unchanged with a warning; do not fill NaNs or let them
poison unrelated series. Signed inverse-FFT results are retained without abs
or post-filter clipping, preserving DC. This differs from continuous velocity's
fourth-order 15 Hz Butterworth filter and the existing velocity per-beat
MATLAB magnitude reconstruction. The persisted `band_limited_harmonics`
setting controls both velocity and gaps and excludes DC.

`/ProcessingSpectrograms/Gaps/PerBeat/SVDContinuous` independently appends the raw
interpolated beats in chronological order to form `(beat * phase, freq)`,
subtracts each bin's temporal mean, removes the first two modes in descending
singular-value order, and restores the mean. Persist the concatenated 2-D
result with `dimDesc=["time", "freq"]`, shape `(beat * phase, freq)`, without
unfolding back into beats. The matrix uses normalized beat phase rather than
uniform acquisition time; SVD does not require a sampling interval. Entirely
invalid bins are excluded; fit rows must be finite across all remaining bins.
Excluded rows/bins remain unchanged and NaNs are never filled. Require at least
three usable matrix dimensions to remove two modes while retaining one;
otherwise publish unchanged raw data with a warning. Filtered undershoots
remain signed to preserve the baseline. Settings, singular values, fitted mean,
support counts, matrix ordering, and filtering availability are recorded in
attributes. The mode count is configurable in `gaps/signal_processing.py`.
Relative sample time is `sum(previous beat durations) + phase * current beat
duration`; spacing is uniform within each beat but may differ between beats.
Source beat count, samples per beat, and within-beat phase coordinates are
recorded to make this continuous ordering explicit.

These paths replace the former root gap datasets; no legacy aliases or removed
MAD/Hampel/distribution/streak filter products are emitted.

Timing uses scalar `sampling_freq` and `batch_stride`, resolved from HDF5 or HD
sidecar configuration. The frame interval is
`batch_stride / sampling_freq` seconds. Spatial calibration is a scalar JSON
value at `/HD_parameters` containing `pixel_pitch: [x, y]` in metres. Current
profile consumers require the two pitches to be approximately equal before
using the isotropic value.

## DV input

DV spatial arrays are native `(y, x)` values: dimension 0 is Y and dimension 1
is X, including `/segmentation/Retina/vessel_segmentation_logits`. Their row
order has Y inverted relative to the published lower-left EF frame. Preserve
that native order when reading and binarizing masks for calculations alongside
raw HD images. EF segmentation serialization flips spatial Y once and then
transposes to published `(x, y)`; do not add a separate input-side Y flip.

| Value | HDF5 path | Requirement |
|---|---|---|
| Artery mask | `/segmentation/Retina/artery_mask` | Required |
| Vein mask | `/segmentation/Retina/vein_mask` | Required |
| All retinal vessels | `/segmentation/Retina/vessel_segmentation_logits` | Optional for existing analyses; required for capillary analysis zones, binarized with `logits > 0.1` |
| Choroid vessels | `/segmentation/Choroid/vessel_mask` | Optional for existing analyses; required for capillary analysis zones, binarized with `value != 0` |
| Labeled vessels | `/segmentation/Retina/labeled_vessels` | Optional |
| Optic-disc mask | `/segmentation/OpticDisc/mask` | Optional if fallback/other measurements suffice |
| Optic-disc center | `/segmentation/OpticDisc/center` | Optional input; canonical source always resolves a center |
| Optic-disc width/height | `/segmentation/OpticDisc/width`, `/segmentation/OpticDisc/height` | Used to reconstruct an ellipse when no mask is available |

`VelocityEstimation/LocalBackgroundDist` comes from the DV sidecar configuration
and defaults to `2`. `pipelines.vessel_inputs.load_retinal_source_data` is the
only canonical alignment boundary: it validates masks against the active HD
spatial shape and may transpose all DV geometry together when axes are swapped.
Do not add isolated transposes in downstream algorithms.

Missing or invalid optic-disc geometry falls back to an image-centered circle
of assumed physical diameter **1.91 mm**. Its center is `(nx / 2, ny / 2)`;
its diameter in pixels is `1.91 / pixel_pitch_mm`, using HD's native isotropic
pitch in mm/pixel. Rasterization uses the physical radius and clips the circle
at image boundaries without resizing it. The fallback requires valid physical
calibration; it no longer scales with the image half-diagonal. That fallback
retains the original vein mask as velocity-background support but disables vein
measurement in the canonical segmentation.

The two additional vessel arrays must be finite, numeric, two-dimensional and
aligned with the other DV geometry. `pipelines.vessel_inputs` applies the same
axis swap to these masks as to the artery, vein and optic disc. The combined
retinal logits include vessel support beyond the artery/vein masks; never
substitute the artery/vein union for this input.

`topology_core` calls `calculations.topology.capillary_analysis_zones` to form
`cutoff_disc AND NOT (dilate(retinal_vessels OR choroid_vessels) OR optic_disc)`
in native `(y, x)` coordinates. The cutoff disc is centered at `(nx / 2, ny / 2)`
and has diameter `0.8 * nx` pixels (production mask images are square). Vessel
support is dilated three pixels in every direction using one pass of a 7x7
square footprint, including diagonal expansion (Chebyshev distance three).
Dilation is clipped at image borders;
the exterior is treated as background. Dilate the original vessel support
before applying the cutoff, so vessels just outside the circle still exclude
their nearby interior pixels. Neither the optic-disc mask nor the cutoff disc
is dilated. For optic-disc exclusion, a usable DV mask takes precedence;
if absent, the physical centered fallback is used even when DV exported only
center/width/height. True pixels are included in capillary analysis. If either
retinal logits or the choroid mask is absent, omit the capillary output and
warn; continue the existing vessel analyses. Malformed masks fail validation
rather than being silently replaced with empty masks.

## Velocity methods

Every run that schedules velocity attempts both workflows, independently:

- `doppler_moments`: derives RMS frequency from moments, estimates local
  background, subtracts it with the established signed RMS rule, and converts
  frequency to calibrated velocity in `mm/s`. Its outputs use `/ProcessingMoments`.
- `frequency_bands` (band ratio): computes `HF / LF` and converts it to RMS
  frequency with `fRMS_Hz = band_ratio_frequency_scale_hz * (HF / LF)`, then
  applies the same vessel-mask dilation, biharmonic background inpainting,
  signed background-difference, and frequency-to-velocity conversion. Its
  outputs use `/ProcessingBandRatio`. The persisted calibration default is `1 Hz`
  per ratio unit; quantity is `physical_velocity`, unit `mm/s`.

`velocity_estimation_method` is no longer a selectable or persisted setting.
Legacy configuration values are ignored. Method identifiers remain in workflow
provenance and in internal method-specific contexts.

Cardiac cycles are detected independently for each workflow **after** velocity
estimation, from its signed, background-corrected whole-vessel velocity in `mm/s`.
The detector filters the raw physical signal at 15 Hz, then detects derivative
upstrokes using both a 95th-percentile height gate and MAD-based prominence.
Both gates apply before cadence estimation and to final peak selection; adaptive
spacing preserves genuine interval variation. The strict global height gate may
miss lower-amplitude beats; resulting long intervals are retained, and suspected
two- or three-period gaps are reported rather than filled with synthesized beats.
Each workflow uses the existing artery, then vein, then full-record fallback.
Its own cycle boundaries, durations, and spectral timing are persisted under its processing
group's `CardiacCycle` family and used throughout its downstream analyses.
Timing provenance records `cardiac_cycle_detection_quantity="physical_velocity"`
and `cardiac_cycle_detection_unit="mm/s"`, plus the workflow and vessel source.
Continuous velocity outputs filter the raw physical signal separately, without
passing the already filtered detector signal through a second filter.
The workflows must have identical `(frame, y, x)` dimensions.

All selected downstream analyses are recalculated in isolated run state for
each successful estimate, using shared segmentation/topology and workflow-specific
cycle boundaries. Failed inputs, cycle detection, or estimates prevent that
method from producing a workflow. A downstream failure rolls back only the
failed pipeline's writes and skips its dependent consumers. Earlier results,
including velocities, remain in the HDF5; independent pipelines continue in
both workflows. If both workflows fail, the run reports failure after safe
pipelines finish and exposes its partial HDF5 path in `RunResult.outputs`.
Root `execution_variant_failures` records a JSON map of method to error;
`execution_variants_completed` lists methods without failures. Root
`pipeline_status`, `pipeline_failures`, and `pipeline_skips` record shared and
per-workflow outcomes. Velocity-analysis data is flushed before optional pulse
PNG, segment-map AVI, and velocity-signal figure exports; export failures are
warnings and do not invalidate the scientific products.

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

- lazy `ImageMaps` for moments, LF/HF bands, and optional exported spectograms;
- aligned boolean artery/vein/background masks and optional labels;
- resolved optic-disc geometry;
- optional aligned retinal-logit and choroid exclusion masks;
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

- `/ProcessingMoments/Velocity/{Continuous,PerBeat}/{Global,Segments}/{Raw,Filtered}`;
- `/ProcessingMoments/Heartbeat` and `/ProcessingMoments/FrequencyMaps`;
- `/ProcessingMoments/VelocityProfiles` and `/VelocityProfilesFFT`;
- `/ProcessingMoments/Metrics/{waveform_shape_metrics,absolute_waveform_metrics,lowrank_waveform_decomposition}`;
- `/ProcessingMoments/SpatialGradientMetrics` for a directly selected spatial-gradient
  target, and `/ProcessingMoments/BloodVolumeRate` for blood-volume-rate products;
- `/Segmentation` for aligned masks, topology, areas, and lumen geometry;
- `/Segmentation/Gaps/value` for central support outside dilated vessels and
  the optic disc, when its inputs are available;
- `/ProcessingSpectrograms/SpatialAverage` and `SpatialMedian` for full-frame
  spatial power reductions, shaped `(time, freq)`;
- `/ProcessingSpectrograms/Gaps/Raw` for shared spatial mean gap band power, shaped
  `(time, freq)`, when spectra and nonempty application support are available;
- `/ProcessingSpectrograms/Gaps/PerBeat/Raw` and `BandLimited` for moments-timed
  beat spectra and independent FFT filtering, shaped `(phase, freq, beat)`;
- `/ProcessingSpectrograms/Gaps/TimeAggregated/Median` for the NaN-aware median
  pooled over raw phase/beat samples, shaped `(freq,)`;
- `/ProcessingSpectrograms/Gaps/PerBeat/SVDContinuous` for SVD-filtered
  concatenated beats, shaped `(beat * phase, freq)` with `(time, freq)` axes;
- `/Meta` for provenance and selected pass-through data.

The `/ProcessingBandRatio` tree mirrors the selected `/ProcessingMoments` families for the
band-ratio workflow. Shared spectra have their own `/ProcessingSpectrograms` root
and are not duplicated inside either workflow. `/Segmentation` and acquisition
`/Meta` remain shared. These roots replace `/Processing` and `/ProcessingAlt`;
new outputs do not emit aliases at the former paths.

Per-beat durations are published once per workflow at
`/ProcessingMoments/CardiacCycle/Systole/CycleDurationSeconds/value` as a `(beat,)`
vector. The
`/ProcessingMoments/Velocity/PerBeat` family contains beat-aligned velocity waveforms,
not a duplicate duration dataset.

Velocity signal datasets are named `ArterySignal` and `VeinSignal` directly,
without a `/value` suffix. Both `/ProcessingMoments` and `/ProcessingBandRatio` use the same
structure. Continuous global signals have `(frame,)` axes and continuous
segment signals `(frame, branch, radius)`. `Filtered` is the fourth-order
15 Hz Butterworth forward/backward filter applied once to each continuous
recording, using MATLAB-style 12-sample odd reflection padding. Recordings
with at most 12 samples retain the established unchanged-signal fallback.
Both raw and filtered recordings are subdivided using the same zero-based
cardiac boundaries. Global per-beat axes remain `(beat, sample)` and segment
per-beat axes `(sample, beat, branch, radius)`.

Each per-beat `Raw` and `Filtered` group additionally has
`BandLimited/{ArterySignal,VeinSignal}`. Its reconstruction retains DC plus
12 harmonics by default, doubles retained non-DC coefficients, and takes
`abs(ifft(...))` to match EyeFlowMatlab. Magnitude reconstruction can change
the temporal mean and removes signs; it differs from the signed mode used
for gap power. Continuous velocity has no `BandLimited` output. Per-beat
global FFT magnitude and phase remain available at
`{Raw,Filtered}/{FFTAbs,FFTPhase}/{ArterySignal,VeinSignal}`.
Quadrants use `Quadrants/<region>/<signal>` below the corresponding signal
representation group. `VelocityPerBeat` and `VelocityPerBeatSafe` are no
longer emitted; no duplicate safe calculation remains. Waveform metric
consumers use the filtered per-beat branch, including its bandlimited
representation. Blood-volume rate uses raw per-beat segments.

The explicit calculation API is:

```python
from calculations.math.bandlimiting import band_limit

velocity = band_limit(beat, harmonics=12, reconstruction="matlab_magnitude")
power = band_limit(beat_power, harmonics=12, reconstruction="signed", axis=0)
```

`harmonics=0` retains only DC; counts above Nyquist clamp to Nyquist. Signed
reconstruction retains Nyquist once, while the MATLAB magnitude recipe doubles
every retained non-DC coefficient, including Nyquist. The app setting
`band_limited_harmonics` is validated before execution and shared by GUI/CLI,
both velocity workflows, and gap processing. API callers can also specify the
count directly. Filtering/reconstruction mode, count, timing, units, and axes
are recorded in output attributes.

The numerical API default and integer-scalar validation belong to
`calculations/math/bandlimiting.py`. Application defaults and JSON setting
validation belong to `app_settings.py`; the run service validates settings
before replacing an output directory and pipelines pass the count explicitly
to calculations. Neither layer imports the other to resolve its default.
Tests keep the scientific API, application fallback, and shipped JSON defaults
aligned, while cold-import tests enforce that settings do not load calculations
or SciPy and bandlimiting does not load application modules.

`/ProcessingMoments/Maps/VelocityAverageMasked/value` is the temporal mean of the
background-subtracted velocity used by the vessel waveform analysis. The
unmasked velocity average and intermediate delta-frequency average are not
persisted.
All two-dimensional datasets under `/ProcessingMoments/Maps` are serialized in the
same lower-left `(x, y)` image frame as
`/Segmentation/Artery/BranchLabelMap/value`: the internal `(y, x)` array is
flipped vertically and transposed before persistence.

The gaps mask used for capillary analysis uses the same lower-left `(x, y)`
frame and carries source/binarization/fallback, cutoff, and dilation provenance,
`unit="1"`, and `dimDesc=["x", "y"]`.
The HDF5 writer stores this boolean mask as `uint8` values 0/1 with
`original_class="bool"`. Native read-only support is available from
`topology_core` under its declared `capillary_analysis_zones` DAG product;
it remains shared across moments and bandratio workflows. Gap-spectrum
preprocessing runs once in the velocity runner alongside full-frame spatial
mean/median reductions, followed there by moments-timed per-beat spectra and
independent Fourier/SVD filtering under `/ProcessingSpectrograms`.
Physical capillary velocity conversion is not implemented yet. Its source is
`/spectrograms`, with exported `(time, freq, y, x)` axes and explicit spatial
alignment as described above.

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
`/ProcessingMoments/SpatialGradientProfiles` and
`/ProcessingMoments/SpatialGradientMetrics` are not written. Selecting
`spatial_gradient_moment0` directly publishes both families. Gradient-derived
blood-volume-rate datasets identify these inputs as `transient_run_state` and
do not expose dangling HDF5 source-path attributes.

The release-default selection disables `gradient_edges`, velocity profiles,
and the waveform `velocity_profile_analysis` option. Its HDF5 therefore omits
`/ProcessingMoments/SpatialGradientProfiles`,
`/ProcessingMoments/SpatialGradientMetrics`,
`/ProcessingMoments/BloodVolumeRate/{Artery,Vein}/{dynamicEdges,staticEdges}`,
`/ProcessingMoments/VelocityProfiles`, and
`/ProcessingMoments/VelocityProfileAnalysis`. Mask-derived `maskedEdges` and
`totalMaskedEdges` blood-volume-rate datasets remain selected.

For the weighted profile-analysis schema, see
[velocity-profile analysis](velocity_profile_analysis.md).
