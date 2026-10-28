# Data contracts

This is the authoritative documentation layer for data moving into, through,
and out of EyeFlow. Executable path constants in `src/input_output/schema/`
remain the final source of truth.

## Article figure formatting

`src/plotting/formatting.py` provides general-purpose, opt-in
Matplotlib helpers for article and non-article figures. It imposes no font,
backend, export format, or global settings:

- `TickFormat` and `apply_tick_format(axis, ...)` configure decimal or
  significant-digit precision, major interval density, integer ticks,
  scientific notation/offsets, and tick direction, length, width, and padding.
  Defaults retain existing locators and formatters, including logarithmic ticks
  and cardiac-phase fractions. For calibrated bin axes, pass `positions` and
  physical `values`; labels use the values without changing bin order or limits.
  `format_tick_labels()` provides the same precision rules without touching axes,
  including suppression of negative zero.
- `LegendLayout` and `add_legend()` share placement, columns, frame visibility,
  anchor, handle lengths, and spacing. Existing handles, labels, and theme fonts
  are preserved unless overridden. Axes-level and figure-level legends are supported.
- `ColorbarLayout` and `add_colorbar()` share orientation/location, relative
  padding/fraction/shrink/aspect/anchors, label padding, optional text sizes,
  outline appearance, and tick formatting. Both attached bars and explicit
  standalone `cax` are supported; explicit axes retain their geometry. Colormaps,
  normalization, physical units, and data limits remain owned by the caller.

Article and diagnostic producers use these helpers for their existing legends
and colorbars. Layout fields left as `None` delegate to Matplotlib; explicit
keyword overrides take precedence over reusable layouts. No additional EPS
exports or article styling are introduced by the general helpers.

`src/plotting/article.py` owns shared, opt-in article styling.
`ArticleStyle` configures Times New Roman text, Times-compatible STIX math,
separate axis-label/tick/legend/title/annotation sizes, borders, tick rotation,
label/tick padding in points, and optional X/Y data margins. `scaled_fonts()`
keeps text visually consistent across different canvas widths. The standard
style uses 10 pt labels/legends, 8 pt ticks, and 0.4 pt borders; the compact
style retains the 8 pt labels/legends used by velocity and flow signal plots.
Line widths for signals, secondary curves, references, event markers, and fits
also come from the style. Envelope color/opacity, analysis-region shading,
and histogram fill/edge appearance are shared settings. `envelope_kwargs()`,
`region_kwargs()`, and `histogram_kwargs()` supply Matplotlib artist options.
Callers assign these roles explicitly, preserving dashed/dotted signal styles,
scientific colors, map overlays, and dispersion values.

The spectrum preset retains darker translucent envelopes (gray 0.35, alpha 0.45)
and softer analysis-region shading (gray 0.5, alpha 0.1); compact velocity/flow
plots retain opaque gray 0.8 envelopes. Diagnostic presets retain their larger
line widths. Font scaling does not alter widths or opacity. The existing
spectrum fill raster layer continues to composite overlapping transparency
for EPS export. `prepare_article_fill_layer()` shares this behavior with
velocity/flow curves if a translucent envelope style is selected. Only the
fill layer is rasterized; signals and text remain vector in EPS. Opaque
velocity/flow envelopes do not require this layer.

`create_article_figure()` supplies fixed canvas dimensions, constrained layout,
and configurable outer padding in inches, or full-bleed square axes for masks.
Its `signal_padding=(horizontal, vertical)` tuple independently enables automatic
data margins on each axis. The default `(True, True)` retains Matplotlib's
defaults; `(False, True)` removes horizontal padding, and `(False, False)`
removes both. This is independent of outer `padding_inches`. Later explicit
axis limits or `ArticleStyle.x_margin`/`y_margin` settings take precedence.
`apply_article_style()` also accepts existing multi-panel figures or selected
axes. Apply it after adding labels, legends, annotations, and colorbars.
Neither helper changes global Matplotlib settings. Writers serialize the
already formatted figure without imposing scientific presentation choices.
`ArticleStyle` composes `x_ticks`, `y_ticks`, `legend_layout`, and
`colorbar_layout` from the general helpers. Build legends/colorbars with the
chosen layout before applying article typography. Article font sizes applied
last take precedence over general layout font sizes. Spectrum frequency labels
remain explicitly mapped, while its automatic power ticks use the shared
integer/precision policy. Scientific frequency mapping and peak-label collision
selection remain in the spectrum producer.

Article PNG/EPS pairs explicitly use
`input_output.writers.article.save_article_figure_pair(output, fig, name)`.
The caller builds and styles one figure; the exporter draws its layout once
at the requested DPI (320 by default), freezes axes positions, and writes both
formats with the same complete canvas bounds and background. Outer padding is
preserved even if global savefig settings request tight cropping or transparency.
Acquisition prefixes, namespaces, and legacy vessel labels still come from
the existing format writers. The exporter closes the figure after both saves,
including on failure. Gap/spectrum/metric pairs, velocity signals, flow signals,
and lumen-diameter histograms use this helper. General PNG writers and diagnostic
PNG-only outputs keep their own export/cropping policy; no EPS is added implicitly.

Consumers include gap spectrograms/spectra and phase metrics, velocity signal
plots, blood-volume-rate signals and lumen-diameter histograms, velocity
waveform/spectral/correlation diagnostics, map colorbars and transverse profile
maps, and spatial-gradient lumen-size curves. Pipeline code retains canvas
geometry for diagnostic panels, data limits, units, grids, scientific colors,
and overlays; styling preserves hidden map axes and borders.

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
| Low-frequency band | `/band_<i>_<n>_<m>` | First video after sorting numeric `(n, m)` ranges; index `i` is ignored |
| High-frequency band | `/band_<i>_<n>_<m>` | Last video in the same ordering; at least two matching datasets are required |
| Frequency-resolved capillary source | `/spectrograms` | Optional gap-spectrum preprocessing source; exact HD root dataset name |
| Flat-field zeroth moment | `/moment0ff`, legacy `/M0FF` | Optional |
| Registration | `/registration` | Optional HD metadata; copied to output `/Meta/registration` when `copy_hd_passthrough` is enabled |
| Zernike coefficients | `/zernike_coefs_radians` | Optional HD metadata; copied unchanged when `copy_hd_passthrough` is enabled |

Moment and band datasets must be 3-D. The two bands must be numeric and have
identical shapes; their chunk values must be finite. During velocity estimation,
negative samples in either LF or HF are replaced with zero, with a warning naming
the dataset and frame chunk, before ratio calculation; source videos are unchanged.

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
the same orientation for this example. `/Segmentation/Gaps/Mask` and shared
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

Spectral frequency calibration uses the original FFT length `B=batch_size`.
NumPy's `np.fft.fftfreq(B, d=1/fs)` supplies the signed FFT frequency grid;
HD metadata is still required to map that grid to the exported planes.
HD selects FFT indexes `[:int(B*spectrogram_frequency_end)]` followed by
`[int(B*spectrogram_frequency_start)+B//2:]`, then averages consecutive groups
of `spectrogram_frequency_binning` samples. Read the selection bounds from
`HD_parameters`; if either bound is missing or cannot be determined, use
`start=0`, `end=0.5`, meaning the complete signed FFT spectrum. Do not infer a
quarter-spectrum crop from the exported plane count. Resolve `batch_size`
from HD parameters with config fallback. Derive the averaging factor from the
selected sample count and exported plane count, checking any recorded binning.
For the default 256-to-64, four-sample full-spectrum export at 37037 Hz,
centers start at 217.013671875 Hz, with spacing 578.703125 Hz within each sign
block. DC belongs to the first average but is not its center; the negative
Nyquist sample likewise belongs to an average centered above `-fs/2`.
The frequency planes remain in exported FFT order. Sampled registered spectra include negative values;
the input adapter permits them and does not clamp or apply the LF/HF bands'
non-negative validation to this source. HD starts from `abs(spectrum_f)**2`,
but its subsequent subpixel Fourier registration of these spatial power maps
can ring and introduce negative undershoots. Linear resizing and frequency-bin
averaging occur afterward and do not enforce nonnegativity.

`pipelines/low_freq_spectra/spectrograms.py` reads each HD spectral frame once and
publishes `/ProcessingSpectrograms/Global/SpatialAverage` and
`/ProcessingSpectrograms/Global/SpatialMedian` as `(time, freq)` arrays. These are
the arithmetic mean and median over every exported spatial pixel, independent
of gap support. Any nonfinite pixel produces NaN for that reduction's
time/frequency pair. Finite negative reductions are clamped to zero after
aggregation. Frequency order and source power units are retained. Outputs are
written after velocity-workflow failure cleanup, so shared spectra survive a
failed moments workflow when the band workflow succeeds. They are produced
once, with no `/ProcessingBandRatio/Spectrograms` duplicate.

`pipelines/low_freq_spectra/preprocessing.py` prepares optional gap support for
the same streaming reduction pass, which computes the spatial median
of the exported band-power samples included by the application mask, for every
time/frequency pair. Only included pixels enter the median; excluded
pixels are not zero-filled. The `low_freq_spectra` runner invokes this once using
native gap support from `topology_core`, and publishes
`/ProcessingSpectrograms/Gaps/SpatialMedian` with `dimDesc=["time", "freq"]`,
preserving HD frequency-bin order. Finite negative
spatial medians are replaced with zero after reduction, before beat subdivision;
individual pixel samples are not clipped before taking the median. For an even
number of included pixels, use the mean of the two middle values. Gap products
record `quantity="spatial_median_band_power"` and `spatial_aggregation_method="median"`.
Units are copied from the source dataset when supplied,
otherwise `a.u.` with uncalibrated-power provenance. Float64 reduction is
serialized as float32 by the standard writer. Any nonfinite included sample
produces NaN for its time/frequency pair; nonfinite excluded pixels have no
effect. If native or interpolated gap support is empty, log a
warning and skip the entire gap-analysis stage without reading spectral
gap samples or writing `/ProcessingSpectrograms/Gaps/SpatialMedian`. Full-frame spatial
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
band-ratio cycles. If neither workflow has detected intervals, spatial/gap
spectra remain available and per-beat products are omitted with a warning.

`/ProcessingSpectrograms/Gaps/PerBeat/Raw` contains `(phase, freq, beat)` spatial median band
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

`/ProcessingSpectrograms/Gaps/TimeAggregated/MedianPhaseBeat` contains a float32 `(freq,)`
pooled median of `/ProcessingSpectrograms/Gaps/PerBeat/Raw` over its phase and
beat axes. Ignore NaNs; a frequency with no valid samples remains NaN without
an all-NaN warning. Every valid phase/beat sample has equal weight, with no beat
duration weighting and no median-of-beat-medians reduction. Preserve frequency
order and power units. Attributes record the source path, aggregation axes,
weighting, source dimensions, and valid sample counts per frequency. Produce
this shared dataset whenever raw per-beat gaps are available, with no
`ProcessingBandRatio` duplicate or fallback to continuous/filtered gap spectra.

`/ProcessingSpectrograms/Gaps/TimeAggregated/MedianBeat` contains float32
`(phase, freq)` power from the median of `/ProcessingSpectrograms/Gaps/PerBeat/Raw`
over its beat axis only. Each valid beat has equal weight at a given phase and
frequency, with no duration weighting. Ignore NaNs; a phase/frequency pair with
no valid beat samples stays NaN without an all-NaN warning. Preserve phase
coordinates, the excluded endpoint, frequency-bin order, power units, and timing
provenance. Attributes identify the raw source, `aggregation_axes=["beat"]`,
weighting, and source phase/beat counts. Both median products are emitted
whenever raw per-beat spectra are available.

`/ProcessingSpectrograms/Gaps/TimeAggregated/MedianPhaseBeatEarly` and
`MedianPhaseBeatLate` are float32 `(freq,)` spectra calculated from `MedianBeat`.
First take its beat median, then take a NaN-aware median across the early
phase interval `[0,1/3]` or late interval `[1/3,1]`, independently per frequency.
Use the same phase-window helper as the early/late metric summaries. A sample
at exactly one third belongs to both intervals; retain the source grid's
excluded endpoint without adding a sample at phase one. The normal 128-phase
grid selects 43 early and 85 late samples. Each valid phase sample has equal
weight regardless of how many valid beats contributed to that `MedianBeat`
cell. A frequency with no valid samples in a window remains NaN without a
RuntimeWarning. Preserve frequency order, power units, and cardiac timing
provenance. Attributes identify `MedianBeat` as the source, `dimDesc=["freq"]`,
`aggregation_axes=["phase"]`, closed interval endpoints, selected and total
source phase counts, and valid sample counts per frequency. These products
do not pool the raw phase/beat samples and do not retain phase-axis coordinates.

`/ProcessingSpectrograms/Gaps/TimeAggregated/StdBeat`, `StdBeatEarly`, and
`StdBeatLate` store float32 `(freq,)` population standard deviations across
beats for the whole cycle and the same closed early/late phase intervals.
Use the existing `/Gaps/PerBeat/Raw` array: first take a NaN-aware phase median
within each beat and frequency, then `nanstd` across those beat spectra with
`ddof=0`. Do not take STD over pooled phase/beat samples or over `MedianBeat`.
One valid beat gives zero STD; no valid beats gives NaN. Retain full frequency
order and raw power units, without plot normalization or peak cropping.
Attributes record the raw per-beat source, phase-median reduction followed by
beat STD, `ddof`, window endpoints, source phase/beat counts, and valid beat
counts per frequency. Reuse the existing beat division without repeating it.

When the HD spectral video is available, `low_freq_spectra` also exports each
available gap product to both `png/spectrograms/` and `eps/spectrograms/`:
`<stem>_spectrograms_gaps` displays `/ProcessingSpectrograms/Gaps/SpatialMedian` as a
linear viridis heatmap with horizontal elapsed-time ticks labeled `Time (s)`,
using the HD frame interval `batch_stride / sampling_freq`, signed frequency
ticks rounded to whole kHz values, and a vertical power colorbar. A matching
`<stem>_spectrograms_gaps_monochrome` uses white-to-black grayscale with a mild
gamma contrast adjustment. The final negative-frequency bin is labeled to make
the upper bound explicit; replace a neighboring tick when needed to avoid crowding.
`<stem>_spectrogram_gaps_median` plots
`/ProcessingSpectrograms/Gaps/TimeAggregated/MedianPhaseBeat` in the original HD
positive-then-negative order, with a thin black dotted horizontal line at zero
power, labeling frequency in signed kHz with vertical
tick labels. Tick selection keeps the sign boundary visible and avoids crowded
labels. Heatmaps use a 12 cm wide, 2.65:1 canvas; the median spectrum uses an
8 cm wide, 1.618:1 canvas. All use 10 pt Times New Roman axis captions and legends,
8 pt tick labels for readability at article column width,
with horizontal Y-axis tick labels, 0.4 pt black borders, and 320 DPI.
Metric phase plots also use Times New Roman text, with Times-compatible STIX
math lettering for subscripted metric symbols (`A_0`, `A_1`, `A_2`, `R_1`,
`R_2`, `Z_1`, `Z_2`) and the italic `t/T` in the cardiac-phase axis label.
`<stem>_gaps_mask` displays the combined `/Segmentation/Gaps/Mask` mask as white
gaps on black in its published lower-left spatial frame, on an 8 cm square
canvas without ticks, labels, borders, or surrounding white margins.
Missing products omit only their own figures;
the median figure requires detected cardiac timing, while the mask does not
require nonempty gap support. Figures do not change any persisted data.

The full-domain heatmap and median-spectrum axes read averaged signed centers
from `/ProcessingSpectrograms/FrequencyBinTable`, preserving source order.
The table averages the original signed FFT coordinates using exactly the HD
selection and binning above. For a full export with averaging factor `q`, centers
have spacing `q*fs/B` Hz and offset `(q-1)*fs/(2*B)` from the former DC-zero grid.
Read `fs` in Hz from `/Meta/SamplingFrequencyHz/value` and divide averaged
coordinates by 1000 for kHz. The spectral
runner ensures this metadata is available from HD `sampling_freq` (HDF5 parameters
preferred, HD config fallback) before figure export; frequency calibration does
not require `batch_stride`. Velocity analysis uses the same schema-owned metadata
path. Missing/invalid sampling or FFT calibration omits frequency figures and
metrics with a warning while retaining available power products, masks, and filters.

The shared spectral runner also writes `/ProcessingSpectrograms/FrequencyBinTable`
as a float32 `(N,)` array of frequencies. Array index `k` identifies source bin
`k`; no separate bin-index column is stored. `unit="kHz"` describes the values
and `dimDesc=["freq"]` identifies the axis.
Frequency values represent means of signed original FFT frequencies, including
the average containing DC and the negative-frequency averages. Attributes record
zero-based indexing, source bin count, HD source path, original FFT count and its
parameter source, averaging factor, selection bounds and fallback provenance,
sampling frequency and its shared metadata path, and bin order. The table is written once
when valid HD spectra and calibration are available, independently of gap
support, detected cardiac timing, and figure export. Missing/invalid calibration
omits the table with a warning while preserving available power products.
The conversion functions live in `calculations/math/frequency.py`, shared by
table publication, spectral metrics, and artifact axes.

Three additional monochrome PNG/EPS figures use these persisted spectra:
`<stem>_gaps_positive_freq_median_phase_early_late` overlays `MedianPhaseBeatEarly`
and `MedianPhaseBeatLate`; `<stem>_gaps_positive_freq_median_phase_beat` plots
`MedianPhaseBeat`; and `<stem>_gaps_positive_freq_median_phase_all` overlays all
three. These legacy `positive_freq` filenames now display the full signed
spectrum in the stored positive-then-negative HD bin order, with calibrated
kHz tick labels from `FrequencyBinTable`. Include ticks at both stored peak
frequencies, prioritizing these labels over overlapping regular or boundary
ticks. Use whole-number tick labels on both axes. The labeled axis ticks at
`/ProcessingSpectrograms/Gaps/Metrics/index_max1` and `index_max2` identify the
peaks without additional reference lines or border notches. Do not shade
the outer, low-interest regions. Require valid peak indexes in their respective
bin-domain halves; do not infer replacements for missing or invalid indexes.
Normalize each curve by its own `Hv`, using the same NaN-aware central-neighborhood
median as A0. Read `valley_neighborhood_length` from the stored metric metadata;
when absent, use the shared default central-quarter window (bins 24 through 39
for 64 source bins). Divide the matching STD by that same Hv. Shade exactly
those bins, from the first bin's lower edge to the last bin's upper edge, with
a faint gray background to identify the Hv estimation region. Reuse the metric
window selection and rounding; do not expand around NaNs. Figure export remains read-only.
Power remains linear, with the existing dotted power-zero reference.
The central curves remain the existing HDF5 medians; dispersion comes from the
corresponding `StdBeat`/`StdBeatEarly`/`StdBeatLate` datasets.
Use solid, dashed, and dotted lines for median, early, and late respectively.
Draw solid, translucent gray envelopes at normalized power plus/minus STD,
clearly darker than the Hv-region background;
overlaps darken gently, while the black signal curves retain their line styles.
Omit envelope outlines and interior patterns. Keep
lower envelopes below zero without clipping and draw a thin dotted zero line.
Remove left/right data padding while keeping vertical room for the zero reference,
frequency in kHz, and `Power (a.u.)`. Use lowercase signal names in each figure's
legend, with the same font size as the axis labels. Use a 16 cm wide canvas,
about 6.4 cm high (30% taller than the original full-spectrum plots), with
about 13.3 pt axis captions and legends and 10.7 pt ticks. Scale these sizes
with canvas width to match the `spectrograms_gaps` font sizes when both are
displayed at the same width. Retain the DPI, typeface, and border style.
Place the legend centrally above the valley to leave both peaks visible.
Composite background and envelope fills on white
in a single raster layer for matching PNG/EPS appearance; signal curves, axes,
and text remain vector graphics in EPS. Missing required spectra/STD
omit the corresponding figure; an all-NaN neighborhood or nonpositive normalization median
cannot be normalized and is omitted with a warning.

Normalization changes figure coordinates only; the source power and STD
datasets remain unnormalized.

`/ProcessingSpectrograms/Gaps/Metrics` groups all gap-spectrum metrics directly,
without `/value` suffixes or workflow copies. The former `TimeAggregated/Metrics`
group is no longer emitted. Twelve pooled scalar datasets analyze the
`MedianPhaseBeat` spectrum in stored bin order: `index_max1`, `max1`,
`index_max2`, `max2`, `index_min1`, `min1`, `index_min2`, `min2`, `A0`, `A`,
`A2`, and `R1`. Heights are internal calculations; `h0`, `h1`, `hv`,
and `rho_h` are not output datasets.
For a domain of `N` bins, peak 1 is the maximum over zero-based bins 0 through
`floor(N/2)-1`; peak 2 is the maximum over `floor(N/2)` through `N-1`.
Include endpoints, ignore NaNs for peak selection, and choose the
lowest index on ties. An all-NaN half has a NaN peak index and value. Valley
indexes are the two central bins, `floor(N/2)-1` and `floor(N/2)`; valley values are read directly from those bins.
Internally, `h0=(max1+max2)/2`, propagating missing inputs. All three areas use
the calibrated spacing `df_khz=q*fs/(B*1000)` and the common peak-to-peak width
`w_khz=(index_max2-index_max1)*df_khz`.
`A2=h0*w_khz`;
`A=trapezoid(spectrum[index_max1:index_max2+1], dx=df_khz)`, propagating any NaN
in that interval. Both peak samples are included with trapezoidal half weights
at the endpoints; the interval width counts the spaces between bins, without
an extra bin. Integrate with positive spacing in the stored HD bin order,
equivalent to an unwrapped kHz axis, without sorting the signed display axis
or wrapping around the array. The valley neighborhood is centered between the two central bins.
With `valley_neighborhood_length=None`, its length is computed as one quarter
of the actual bin domain, rounded to the nearest even length (ties upward), with
a minimum of two samples: 32 bins selects 8 samples (bins 12 through 19), and
64 bins selects 16 samples (bins 24 through 39). Any domain with at least two bins
is supported, including odd lengths; domains split at `floor(N/2)`.
Explicit neighborhood lengths
accept even integers from 2 through `N` in `run_gap_analysis`,
`pack_spectrum_metrics`, and `analyze_spectrum`. `hv` is the NaN-aware median of
that neighborhood; an all-NaN neighborhood returns NaN. With the same kHz
width used by `A2`, `A0=w_khz*hv`. The area names `A0`
and `A2` replace each other's former meanings. `A` is the former `A1` output.
The pooled dimensionless output is `R1=A/A0`. Ratios require a finite, strictly
positive reference area and finite numerator/result; otherwise they return NaN.
`SEI` is calculated from the phasewise `SE` curve below, not from pooled areas.
There is no pooled-spectrum SEI output. All values are float32; valid index values are integer-valued floats
so missing indexes can be represented by NaN. Indexes and `R1` have unit `1`;
`A0`, `A`, and `A2` have the source power unit times kHz (for example
`a.u.*kHz`), while peak/valley heights retain the spectrum's power unit.
Attributes record definitions,
source, zero-based indexes, selection rules, interval inclusion, neighborhood
length and exact bins, automatic versus explicit selection and rounding,
per-half valid (non-NaN) bin counts, valley aggregation/NaN policy, sampling
frequency and its shared source path, kHz bin spacing, integration method,
and the unwrapped area coordinate. Read calibration from
`/Meta/SamplingFrequencyHz/value`. Missing/invalid calibration skips the Metrics
family with a warning while preserving raw/aggregate spectra and filters;
never substitute a unit bin width for a physical kHz area.
An all-NaN half leaves its peak index and value as NaN and logs an explanatory warning.
A spectrum that is not one-dimensional with at least two bins
produces all twelve NaN pooled metrics and logs a
warning; the median and other gap products remain available. No alternative
workflow copy is generated.

Five phase curves are written in that same group as `A0_t`, `A_t`, `A2_t`,
`R1_t`, and `SE`, with `(phase,)` axes and phase coordinates from `MedianBeat`.
Calculate these curves from `MedianBeat` after its beat median, keeping the
peak indices selected from `MedianPhaseBeat` fixed at every phase. With fixed
width `w_khz=(index_max2-index_max1)*df_khz`, calculate
`A2(t)=w_khz*(S(t,index_max1)+S(t,index_max2))/2`,
`A(t)=trapezoid(S(t,index_max1:index_max2+1), dx=df_khz)`, and
`A0(t)=w_khz*nanmedian(S(t,Nv))` using the same valley neighborhood as the pooled
metrics. `R1(t)=A(t)/A0(t)` and `SE(t)=(A(t)-A0(t))/A0(t)`.
Missing fixed peak indices make every curve NaN. Missing peak heights invalidate
`A2(t)`, interval NaNs propagate through `A(t)`, and an all-NaN valley
neighborhood makes `A0(t)` NaN. Ratios and offsets require finite, strictly
positive `A0(t)` and finite numerator/result; invalid phases remain NaN.
Curve attributes identify their `MedianBeat` source, fixed peak values and
index-dataset paths, pooled peak source, timing, neighborhood, and units.

The dimensionless scalar `SEI=median(SE(t))` follows the temporal definition:
take the NaN-aware median over all valid phase samples of the registered cycle,
with equal weight per source phase sample and no added endpoint. This differs
from a ratio of pooled areas or a ratio of median phase areas. It is emitted
only when phase spectra are available, with source path `/ProcessingSpectrograms/Gaps/Metrics/SE`,
aggregation method, and source/valid/invalid phase counts. All-invalid cycles
produce NaN. There is no fallback to the pooled-spectrum ratio. These validity
checks do not impose an unconfigured low-signal or valley-stability threshold.

For the `A0_t`, `A_t`, `A2_t`, and `R1_t` curves, the same group also contains
scalars `<base>_early`, `<base>_late`, and `delta_<base>`. The `SE` curve uses
`SEI_early`, `SEI_late`, and `delta_SEI`. Early is the NaN-aware median of the
corresponding curve over a configurable closed `early_interval` (default `[0,1/3]`); late uses
`late_interval` (default `[1/3,1]`). Bounds are fractions of a cycle in `[0,1]`,
equivalent to 0–33⅓% and 33⅓–100% by default. A sample exactly at a shared
boundary belongs to both intervals. Preserve the source grid's excluded endpoint;
do not interpolate or add a sample at phase one. The standard 128-phase grid
therefore has 43 early and 85 late samples. Medians weight valid phase samples
equally; all-missing windows remain NaN. Aggregate the ratios themselves,
rather than dividing aggregated areas. The whole-cycle scalar corresponding
to `SE` is `SEI`. Every spectrogram metric uses late minus early, including
A0/A/A2, R1, and SEI. Preserve the signed
difference and propagate a missing early or late median. Window metadata records
phase intervals, included endpoints, source and valid sample counts, source
curve paths, and aggregation policy; deltas identify their late then early
scalar sources and both configured intervals.
Shared helpers in `calculations/math/phase_windows.py` provide interval masks,
NaN-aware window medians, and late-minus-early subtraction. `run_gap_analysis`
and `pack_spectrum_metrics` accept `early_interval` and `late_interval` keyword
arguments. The runner applies the same configuration to early/late spectra,
STD envelopes, metric summaries, and deltas. Empty/all-missing windows retain NaN.
These twenty phase/window datasets, the whole-cycle `SEI`, and twelve pooled scalars give 33 datasets
under `/ProcessingSpectrograms/Gaps/Metrics`. All persisted values are float32;
ratio and offset curves and their summaries/deltas have unit `1`, while area products
have the source power unit times kHz.

Each phase curve also produces `<stem>_gaps_{A0,A,A2,R1,SE}.{png,eps}` under
the shared `png/spectrograms/` and `eps/spectrograms/` folders. Use the same
8 cm wide golden-ratio canvas, thin black border, fonts, and resolution as the
spectral figures. Abbreviate `a.u.*kHz` to `a.u.` in artifact axis labels;
stored area units remain source power times kHz.
Plot normalized cardiac phase `t/T` horizontally with labels
at zero, one third, two thirds, and one. Preserve NaN gaps and the excluded
endpoint; area-axis labels include power units and ratio axes are dimensionless.
Each available curve exports independently.

`/ProcessingSpectrograms/Gaps/QC` contains descriptive diagnostics, separate
from `Metrics`. QC never filters samples, rejects beats, or changes SE/SEI.
There are no implicit pass/fail thresholds. A fully evaluated run writes only
six QC datasets:

- `DataCompleteness/{AcquisitionValidFraction,ValidSEPhaseFraction}`: overall
  finite fraction of reduced acquisition-time spectra and fraction of finite SE
  phases. Counts and complete-beat numbers can be read from the existing sources.
- `ValleyRepeatability/RelativeMAD`: unscaled across-beat valley-height MAD
  divided by absolute median height. Comparisons are at the same cardiac phase. Use the metrics'
  recorded valley neighborhood. Dispersion needs at least two finite beat
  heights; relative dispersion additionally needs a nonzero median. Undefined
  values stay NaN. This is dispersion, not an uncertainty interval.
- `SEIBeatSensitivity/DeltaSEI`: signed omission-minus-full SEI change for
  each complete beat in chronological order. For each omission,
  repeat the production raw-spectrum pooling, reselect pooled peaks, and
  recompute phasewise SE and temporal-median SEI. Match production float32
  reductions/scalar storage. At least two original beats and calibrated SEI
  metrics are required; otherwise estimates are NaN. Array beat indices are
  zero-based, figure beat labels one-based. This measures influence, not a
  confidence interval or proof of motion. Valid phase counts per omission are
  recorded as supporting metadata. Full SEI remains in `Metrics/SEI`; omitted-
  beat estimates and maximum changes are derived from that value and `DeltaSEI`.
- `PotentialMotionContamination/RelativeBroadbandChange` is the median
  across bins of `abs(S[t,j]-S[t-1,j])/median_time(S[:,j])`, using only finite
  differences and positive finite per-bin temporal medians. The score belongs
  to the later frame; the first frame and unsupported scores are NaN.
  It is invariant to positive
  per-bin multiplicative scaling and is a candidate transient indicator, not
  a validated saccade detector; physiological changes can contribute.
- `SpatialSupportQuality/RetainedFieldFraction`: reduced application-mask
  coverage relative to the entire rectangular spectrogram field. Retained pixel
  count is already recorded on the source spectra.
  Missing mask metadata yields NaN. Anatomical suitability requires visual review.

Acquisition and spatial QC remain available without cardiac registration;
phase-dependent outputs require registered beats and valley-dependent outputs
require at least two frequency bins. Missing SE calibration is explicitly
identified in SE-support and influence metadata. Dataset attributes record
definitions, units, axes, source paths, and the descriptive-only policy.
`png/spectrograms/<stem>_gaps_qc.png` shows data support, broadband changes,
per-beat valley heights, relative valley MAD, SEI beat influence, and the retained
mask. Plotting reconstructs frequency-bin counts and valley heights/dispersion
from existing spectra instead of persisting duplicate QC arrays.
It uses acquisition seconds when timing is available, otherwise frame
indices. This diagnostic has no automatic EPS export or article-format styling.

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
and defaults to `2`. `input_output.readers.canonical.load_retinal_source_data`
is the canonical alignment boundary: it validates masks against the active HD
spatial shape and may transpose all DV geometry together when axes are swapped.
`pipelines.vessel_inputs.load_retinal_source_data` then applies the venous
fallback policy. Do not add isolated transposes in downstream algorithms.

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
center/width/height. This original mask is preserved at
`/Segmentation/Gaps/Legacy`. If either retinal logits or the choroid mask is
absent, omit Legacy and warn; continue existing vessel analyses. Malformed
masks fail validation rather than being silently replaced with empty masks.

`topology_core` also produces `/Segmentation/Gaps/GapsFromM0ff` from the temporal
median of HD `/moment0ff` (legacy alias `/M0FF`). Read spatial tiles over all
frames to compute the exact median without loading the complete volume.
Ignore nonfinite temporal samples; pixels with none remain excluded. The
native `(y, x)` shape must match the HD reference, without inferred transposes.
Stretch the finite image's 1st–99th percentiles to `[0,1]` and clip; if these
percentiles coincide, use the finite minimum/maximum. Apply CLAHE with clip
limit 0.01, 256 bins, and kernel dimensions `max(1, dimension // 8)`, followed
by a sigma-1-pixel Gaussian with reflect borders and finite-support weight
normalization. Invalid pixels are median-filled only for CLAHE. Global Otsu
thresholding selects **bright vessels** (`intensity >= threshold`). Dilate
the vessel support once using the legacy 7x7 square footprint (Chebyshev
radius 3), clipping at image borders with exterior background, then invert
to obtain white gaps and exclude invalid pixels. Apply no anatomical or
field-of-view exclusions. Uniform/all-invalid images yield an empty mask.
Persist source, contrast, smoothing, threshold, and dilation parameters as
attributes.

`/Segmentation/Gaps/Mask` is
`NOT ((NOT Legacy) OR (NOT GapsFromM0ff)) AND cutoff_disc AND NOT optic_disc`
and is the sole
downstream gap support, including the `gaps_mask` figure. Reapply the same
image-centered 80%-width cutoff and DV optic-disc mask (or the centered
1.91 mm physical fallback) after combining. OR is applied to the vessel/exclusion
masks before inversion, equivalent to AND of the two saved white-gap masks.
A pixel remains a gap only if both sources include it. These constraints apply only to
Mask; GapsFromM0ff retains its full-field intensity-derived support. When
`moment0ff`/`M0FF` is absent, warn, omit GapsFromM0ff, and use Legacy unchanged.
If Legacy is unavailable but the HD flat-field moment exists, use
GapsFromM0ff with the same disc constraints. Omit Mask only when neither
source is available. Record the contributing paths, disc geometry, optic-disc
source, and legacy fallback status. The old
`/Segmentation/Gaps/value` dataset is no longer produced.

## Velocity methods

Every run that schedules velocity attempts both workflows, independently:

- `doppler_moments`: derives RMS frequency as
  `sqrt(M2 / max(spatial_mean(M0), eps('single')))` per frame, using float32
  epsilon (`1.1920929e-7`) as the denominator floor and preserving NaNs; estimates
  local background, subtracts it with the established signed RMS rule, and converts
  frequency to calibrated velocity in `mm/s`. Its outputs use `/ProcessingMoments`.
- `frequency_bands` (band ratio): computes `HF / LF` and converts it to RMS
  frequency with `fRMS_Hz = band_ratio_frequency_scale_hz * (HF / LF)`, then
  applies the same vessel-mask dilation, OpenCV Telea background inpainting
  with a 5-pixel radius on float32 frames,
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

For band mode, replace negative HF and LF samples with zero and log a warning.
A zero LF denominator, including a clamped negative sample, maps to ratio zero
without performing division. No epsilon is
added. A nonzero ratio beyond finite `float32` range raises a clear error rather
than emitting infinity. Missing exact band paths, mismatched shapes, nonnumeric
data, and NaN/Inf still fail explicitly. Negative finite samples do not skip a
workflow. A failed band
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

`RetinalSourceData` in `src/retinal_source_data.py` is the shared typed input:

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
corresponding processing group and the actual PNG paths registered by its
velocity-analysis variant. Existing filename lookup remains a fallback. The
primary file is `h5/<stem>_EF.h5`.
Displacement magnitude videos use MJPEG AVI files under
`avi/displacement_maps/`, with one file for each selected vessel mask.

All artifact basenames start with `<stem>_`, including PNG, EPS, AVI, MP4,
PDF, auxiliary HDF5, and archive outputs. The shared artifact naming helper
and `OutputManager.path_for()` enforce this at the writer boundary: existing
prefixes are retained once and subdirectories are preserved. The primary
`<stem>_EF.h5` filename is unchanged.

Existing EyeFlow result folders can be migrated with
`python scripts/migrate_artifacts.py RESULTS_PARENT --dry-run`, then
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
- `/Segmentation/Gaps/Legacy` for the original vessel/disc/cutoff support,
  `/Segmentation/Gaps/GapsFromM0ff` for dark regions of the enhanced HD temporal
  median after vessel dilation, and `/Segmentation/Gaps/Mask` for the intersection
  of their gap support (OR of vessels before inversion), restricted to the
  outer disc and excluding the optic disc, used downstream;
- `/ProcessingSpectrograms/FrequencyBinTable` for signed frequencies in kHz,
  shaped `(freq,)` and indexed by zero-based source bin;
- `/ProcessingSpectrograms/Global/SpatialAverage` and `SpatialMedian` for full-frame
  spatial power reductions, shaped `(time, freq)`;
- `/ProcessingSpectrograms/Gaps/SpatialMedian` for shared spatial median gap band power, shaped
  `(time, freq)`, when spectra and nonempty application support are available;
- `/ProcessingSpectrograms/Gaps/PerBeat/Raw` and `BandLimited` for moments-timed
  beat spectra and independent FFT filtering, shaped `(phase, freq, beat)`;
- `/ProcessingSpectrograms/Gaps/TimeAggregated/MedianPhaseBeat` for the NaN-aware median
  pooled over raw phase/beat samples, shaped `(freq,)`;
- `/ProcessingSpectrograms/Gaps/TimeAggregated/MedianBeat` for the NaN-aware median
  over raw beat samples at each phase/frequency, shaped `(phase, freq)`;
- `/ProcessingSpectrograms/Gaps/TimeAggregated/MedianPhaseBeatEarly` and
  `MedianPhaseBeatLate` for phase-window medians of `MedianBeat`, shaped `(freq,)`;
- `/ProcessingSpectrograms/Gaps/TimeAggregated/StdBeat`, `StdBeatEarly`, and
  `StdBeatLate` for population beat STD of whole/early/late per-beat phase medians;
- `/ProcessingSpectrograms/Gaps/Metrics` for pooled scalars, fixed-peak phase
  curves, configurable early/late medians, and signed late-minus-early differences
  for every metric;
- `/ProcessingSpectrograms/Gaps/PerBeat/SVDContinuous` for SVD-filtered
  concatenated beats, shaped `(beat * phase, freq)` with `(time, freq)` axes;
- `/Meta` for provenance and selected pass-through data.

The `/ProcessingBandRatio` tree mirrors the selected `/ProcessingMoments` families for the
band-ratio workflow. Shared spectra have their own `/ProcessingSpectrograms` root
and are not duplicated inside either workflow. `/Segmentation` and acquisition
`/Meta` remain shared. Pixel pitch is stored directly as the scalar dataset
`/Meta/PixelPitch_m` (unit `m`); `/Segmentation/PixelPitch_m/value` is no longer
written. These roots replace `/Processing` and `/ProcessingAlt`;
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
When velocity analysis needs the full `(frame, y, x)` float32 velocity videos,
the run shares a RAM budget of one quarter of available memory across both
workflows (1 GiB if availability cannot be measured). Arrays exceeding that
budget use temporary `.npy` scratch storage beside the EyeFlow output and are
removed when the run finishes, including after a failed pipeline. Scratch
videos are never published as HDF5 datasets or artifacts.
All two-dimensional datasets under `/ProcessingMoments/Maps` are serialized in the
same lower-left `(x, y)` image frame as
`/Segmentation/Artery/BranchLabelMap/value`: the internal `(y, x)` array is
flipped vertically and transposed before persistence.

All three gaps masks use the same lower-left `(x, y)`
frame and carry their source/binarization/fallback and processing provenance,
`unit="1"`, and `dimDesc=["x", "y"]`.
The HDF5 writer stores this boolean mask as `uint8` values 0/1 with
`original_class="bool"`. Native read-only support is available from
`topology_core` under its declared `capillary_analysis_zones` DAG product;
it remains shared across moments and bandratio workflows. Gap-spectrum
preprocessing runs once in `low_freq_spectra` alongside full-frame spatial
mean/median reductions, followed there by moments-timed per-beat spectra and
independent Fourier/SVD filtering under `/ProcessingSpectrograms`.
Physical capillary velocity conversion is not implemented yet. Its source is
`/spectrograms`, with exported `(time, freq, y, x)` axes and explicit spatial
alignment as described above.

Do not hand-copy path strings into new consumers. Obtain them from the schema or
from a pipeline-owned constant when the path family has not yet been centralized.
Writers normalize leading/trailing separators, replace an existing dataset at
the target path, downcast float/complex payloads to 32-bit, preserve integer
arrays, and serialize boolean payloads as `uint8` with `original_class="bool"`.
Writers no longer generate the redundant `nameID` attribute; HDF5 exposes the
absolute dataset path as `dataset.name`. Existing explicitly supplied attributes
remain supported for compatibility with older payloads.
Numeric arrays of at least 4 KiB (after dtype normalization) use lossless gzip
compression at level 9 with the shuffle filter by default. Smaller arrays stay
unfiltered to avoid chunk-index overhead exceeding compression savings. Explicit
dataset storage options override these defaults, including compression for small
arrays. Scalars and strings remain unfiltered; these writes preserve their
existing shapes and dtypes.
Per-beat velocity maps, profiles, and segment masks use gzip level 9. The
float32 map and profile writers leave chunks containing only the exact NaN fill
value unallocated; reading the published datasets returns the same dense shapes
and values. When the redundant bins represent at least 32 KiB of raw data,
the optional velocity-profile FFT stores the unique half of an exactly mirrored
magnitude spectrum in a gzip-compressed
backing dataset and exposes the original full shape through an HDF5 virtual
dataset at the published path. Smaller arrays or arrays whose two halves differ
at the bit level keep ordinary full-array storage.
The backing dataset is named `__<published-name>_unique_frequency_storage` in
the same group. A virtual dataset has no compression filter of its own; its
backing data has the gzip filter.

Fresh runs omit HD `/registration` and `/zernike_coefs_radians` pass-through
copies by default. Set `copy_hd_passthrough` to `true` in the imported settings
JSON, pass `--copy-hd-passthrough` to the CLI, or pass
`copy_hd_passthrough=True` to `resolve_run_spec` to retain them.
The root output attribute records this choice, and the HD source-file path
remains available when copies are omitted. When enabled, large numeric HD
pass-through datasets are copied in chunks and recompressed with gzip level 9;
their values, shapes, dtypes, and attributes remain intact. Unsupported or
small pass-through objects retain their source storage layout.

The output root records source files, selected targets, actual pipeline order,
selected options, requested methods, and workflow status. Each separate workflow
group (`/ProcessingMoments` and `/ProcessingBandRatio`) keeps only these summary
attributes: `velocity_estimation_method`, `velocity_inpainting_radius`, and
`velocity_unit`, plus workflow-specific sources. `/ProcessingMoments` records
`moment0_source_path` and `moment2_source_path` (including actual legacy aliases
such as `/M0` and `/M2`). `/ProcessingBandRatio` records `band_lf_source_path` and
`band_hf_source_path`. Source attributes from the other workflow are excluded.
Other processing-parent attributes are discarded, without relocating them to
new metadata groups. Duplicate heart-rate attributes are omitted in favor of
the cardiac-cycle datasets; sampling frequency, batch stride, and frame interval
remain in their existing acquisition `/Meta` datasets.
Only the fixed summary keys apply to workflow descendants. Dataset-specific
units, dimensions, calibration, metric definitions, and other existing payload
attributes remain local. Equal summary duplicates are omitted at write time;
dataset-local overrides take precedence. No implicit inheritance occurs from
the file root or between workflows. Older explicit inheritance declarations
remain readable. References to processing datasets use the owning namespace.
Persisted RMS-frequency maps use `Hz`; velocity datasets use `mm/s` and
inherit method/calibration provenance. Read effective attributes with
`input_output.schema.provenance.read_dataset_attrs(dataset)` or
`PipelineH5Output.read_attrs(path)`; raw `dataset.attrs` only contains local
attributes. The reader also accepts older files that store provenance locally
and have no inheritance declaration. Legacy dimensionless velocity metadata is
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
