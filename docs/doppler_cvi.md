# Frequency-resolved Doppler CVI proxy

The selectable `low_freq_spectra` pipeline publishes a functional linear-power
proxy. It does not establish choroidal specificity or anatomical OCT CVI equivalence.
It preserves the existing full-frame/gap spectra, filters and spectral-excess products.

## Code ownership

CVI follows the existing `low_freq_spectra` stages and artifact export path:

| Module | Responsibility |
|---|---|
| `runner.py` | Stage dependencies and handoffs for spectra, CVI and binary occupancy maps |
| `per_beat.py` | Reduced-spectrum and spatial-power phase alignment |
| `cvi.py` | CVI and binary occupancy analyses, including input requirements, calculations, result datasets and provenance |
| `spectrograms.py` | Source-spectrogram spatial preprocessing: full-frame mean/median and optional gap median |
| `artifacts.py` | Shared article heatmaps, spatial maps and AVI rendering |
| `figure_export.py` | Snapshotting and isolated background artifact export |

Like `gap_analysis.py`, `cvi.py` is a product owner. Its entry points are
`run_cvi_analysis` and `run_binary_occupancy_analysis`; its numerical helpers
accept arrays and remain independent of scheduling and rendering. The
`binary_occupancy` stage produces spatial maps used by both images and movies.
The shared spectral snapshot carries the persisted CVI products to the existing
artifact renderer; CVI does not have a separate export subsystem.

For each original signed HD frequency bin, spatially intact `/spectrograms`
power is linearly aligned between the same detected complete-beat boundaries
as the gap spectra (moments preferred, detected band-ratio timing as fallback).
The next boundary is an interpolation anchor; the 128 phases are `[0,1)`.
Exact acquisition anchors remain exact, and nonfinite temporal gaps are not bridged.
No cardiac timing is inferred when these workflows have no detected intervals.

Support is the intersection of the image with the closed disc centered at
`((W-1)/2, (H-1)/2)` with radius `0.4*W`. The fixed denominator is its actual
pixel count. Retinal vessels and optic disc remain included. There is no gap-mask
restriction, smoothing, power clipping, log transform or normalization.

At each phase/frequency, a beat image with any nonfinite included pixel is
excluded in its entirety. The threshold reference is the pixelwise median of
valid aligned beats. Pool its finite values over phases and disc pixels with
equal weight per valid phase/pixel. A single Otsu threshold per frequency uses
256 equal-width bins spanning the finite reference minimum and maximum.
`skimage.filters.threshold_otsu` receives the accumulated histogram and bin centers.
Constant reference values give that constant threshold and a degeneracy flag;
no reference gives NaN and an invalid flag.

Each valid individual beat image is compared using **power > threshold**.
Its white count divided by the geometric denominator gives the beat fraction.
The phase/frequency CVI is the median of valid beat fractions, with NaN when
none remain. This is distinct from binarizing the beat-median reference image.

## HDF5 contract

`EyeFlowOutputPaths.active().spectrograms_cvi_root` owns
`/ProcessingSpectrograms/DopplerCVI`. Its children are:

| Dataset | Axes / meaning |
|---|---|
| `CVI` | `(phase, freq)` fractions, unit `1` |
| `PerBeat` | `(phase, beat, freq)` fractions; invalid images are NaN |
| `Thresholds` | `(freq,)`, source power units |
| `Phase`, `Frequency` | normalized phase and original calibrated signed kHz centers |
| `DiscMask`, `PixelCount` | native HD `(y,x)` Boolean support, actual geometric count |
| `ValidBeatCounts` | `(phase, freq)` valid image counts |
| `ThresholdFlags` | `(freq,)`: 0 valid, 1 constant reference, 2 no finite reference |
| `ThresholdReferenceCounts` | `(freq,)` pooled finite reference sample counts |
| `Movies/PhaseMean` | `(phase,y,x)` binary occupancy averaged over valid beats, then the inclusive peak-bin interval |
| `Movies/FrequencyMean` | `(freq,y,x)` binary occupancy averaged over valid beats, then equally over valid phases |
| `Movies/FrequencyMeanOverBins` | `(y,x)` equal-weight mean of finite frequency-movie maps over all signed bins |
| `Movies/ValidFrequencyCounts`, `Movies/ValidPhaseCounts` | counts contributing to the two map reductions |

Masks and movie maps retain HD export coordinates. Numeric outputs use the
existing writer's float32 serialization convention;
interpolation, reference medians, thresholds and fractions use float64 internally.
Dataset attributes record timing sources/boundaries/durations, calibration,
geometry, reduction order, invalid policy, Otsu settings/library versions, tile
size, and the movie aggregation selections.

## Configuration and figures

Defaults and argument validation live inside `low_freq_spectra`; no global
application settings or front-end controls are added. Pipeline API callers can
override calculation parameters with `spectral_graph(cvi_parameters={...})` or pass keyword
arguments directly to `run_cvi_analysis`. These are the defaults:

```python
spectral_graph(cvi_parameters=dict(
    histogram_bins=256,
    working_memory_mb=64,
))
```

`doppler_cvi_heatmap.png` under `png/spectrograms/cvi/` has
its article EPS counterpart under `eps/spectrograms/cvi/`. It uses the shared
spectrogram article dimensions, typography and borders, has no title, and
displays a linear black-to-white grayscale fixed at 0--100% CVI. The horizontal
axis remains cardiac phase. The vertical axis uses ascending calibrated signed
kHz centers, with the display rows sorted accordingly; stored scientific arrays
keep original HD bin order. `spectrograms_gaps_monochrome` uses the same signed
frequency display coordinates while retaining its elapsed-time axis and power
scale. Waveform plots, spatial previews, preview buffers and preview HDF5
datasets are no longer produced.

Two grayscale AVI movies are exported at 10 fps under
`avi/spectrograms/cvi/`: `doppler_cvi_binary_phase.avi` traverses one normalized
cardiac cycle and averages bins `index_max1` through `index_max2`, inclusive;
`doppler_cvi_binary_frequency.avi` traverses all calibrated signed frequencies
in ascending order, averaging over cardiac phase. Stored map arrays retain
original source-bin order, and AVI metadata records the display mapping.

For these movies, each spatially intact linear-power sample is phase-aligned
and compared strictly against the persisted per-frequency CVI threshold. The
finite `PerBeat` cells define valid whole beat images. First average these
binary images over valid beats at each phase and frequency, then average
equally over available bins (phase movie) or phases (frequency movie). This
beat mean is specific to the movies; scientific CVI fractions retain their
original median-over-beats definition. No threshold is applied after averaging.
Missing bins/phases are omitted from these means; no valid contributions give
NaN. Pixels outside the disc remain NaN in HDF5 and black in the videos;
invalid included pixels appear magenta. Every frame uses a fixed 0--100%
grayscale scale with disc boundary and phase/bin/calibrated-frequency labels.
The AVI compression is lossy; the HDF5 fraction maps are authoritative.

`doppler_cvi_binary_frequency_mean` is an additional article PNG/EPS spatial
map under `spectrograms/cvi/`, obtained by averaging all finite frequency-movie
maps with equal bin weights. It retains native HD image coordinates and disc
support and uses a fixed linear black-to-white 0--100% scale, without another
threshold. The PNG/EPS contains only the image: no axes, labels, colorbar, disc
outline or surrounding margins; native pixel aspect ratio is preserved.
Frequencies with no finite value at a pixel are omitted there; no
valid frequency gives NaN. The saved HDF5 map is `Movies/FrequencyMeanOverBins`.
The exporter can also derive this map directly from `Movies/FrequencyMean` in
older results, reading one frame at a time rather than decoding the AVI.

## Resource behavior and validation

Four passes per frequency determine whole-image validity, finite reference
range, the reference histogram, and individual white counts. Constant/empty
references skip the histogram pass. Spatial tiles bound temporary arrays using
`working_memory_mb`, with a minimum of one pixel. The memory target excludes
retained fraction/QC outputs. The input service
provides independent immutable source reads; the spatial median reducer never
partitions CVI pixels. Repeated reads trade processing time for bounded working
memory and require no additional CVI scratch video.

Synthetic numerical tests compare with a fully materialized direct reference,
including NaNs/infinities, signed bins, negative power, rectangular discs,
strict equality, beat aggregation order and tile-size independence. Pipeline
tests cover detected-timing preference/fallback, missing inputs/timing,
serialization, artifacts, and existing spectral products. Real acquisition
validation remains necessary before scientific interpretation.
