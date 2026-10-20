# Pipeline context

## Boundary

Pipelines adapt typed inputs and shared scientific calculations to the runtime.
They declare dependencies/options, orchestrate calculations, own run-local state
keys, pack HDF5 datasets, and export product-specific artifacts. Generic DAG and
execution behavior belongs in `pipeline_engine`; reusable formulas belong in
`calculations`; external path aliases belong in `input_output/schema`.

## Discovery and package shape

`pipelines/__init__.py` scans top-level modules/packages and imports them.
Registration occurs through `@pipeline` or `@registerPipeline`; no registry list
needs editing. A helper-only top-level module must not register a pipeline.
Tutorial pipelines live under `pipelines/tutorials/` and are hidden from normal
selection.

Small pipelines may be one module. Larger pipelines use:

```text
pipeline_name/
  __init__.py       declaration only
  runner.py         orchestration and state ownership
  sources.py        typed input assembly, when needed
  outputs.py        schema packing/artifact export
  calculator.py     product-specific calculation, when not broadly reusable
```

Keep imports in `__init__.py` light enough that catalog discovery remains
predictable and reports missing optional dependencies as descriptor
availability rather than unrelated import crashes.

## Current product groups

| Group | Responsibility and important upstream state |
|---|---|
| `velocity` (hidden) | Physical retinal velocity and cardiac timing; produces `velocity`, `cardiac_cycles`, and typed `velocity_workflows` |
| `low_freq_spectra` (hidden) | Automatically follows velocity once per run when HD spectral video is available; owns shared spatial/gap spectra, beat alignment, filters, and spectral metrics |
| `topology_core` (hidden) | Canonical artery/vein prepared topology and combined legacy/intensity gap mask; produces `prepared_topology` and `capillary_analysis_zones` |
| `velocity_analysis` | Shared per-beat/segment/profile state plus selectable persisted continuous, quadrant, profile-fit, FFT, map, and movie products; requires velocity and topology |
| `spatial_gradient_moment0` | Independently selected moment0 gradient/lumen profiles; requires cardiac cycles and topology, not velocity analysis |
| `blood_volume_rate` | Two option families with different DAG dependencies: gradient edges and mask-derived geometry |
| `waveform_shape_metrics`, `absolute_waveform_metrics`, `lowrank_waveform_decomposition` | Downstream waveform metrics requiring `velocity_analysis` |
| `pdf_report` | Report assembly after waveform and shape outputs |
| `displacement_map` | Separate image-registration path producing magnitude videos and run-local dense-field artifacts |

The declarations in each package are authoritative if this table drifts.

## Shared scientific flow

`vessel_inputs.py` builds canonical `RetinalSourceData`: active HD volumes,
aligned vessel masks, optic-disc geometry, timing/pitch, DV background settings,
and velocity method. The velocity and velocity-analysis pipelines use compatible
sources.

The velocity pipeline estimates the physical velocity field and detects
zero-based cardiac-cycle boundaries. Native geometry lives in
`calculations/topology`; the spatial `SegmentSamplingPlan` is prepared by
`calculations/vessel_segments/sampling` and cached by `topology_core` in run
state. The existing `prepared_topology` DAG key is unchanged.
Velocity analysis consumes both results and produces the common
per-beat, segment, and profile state used by visible products.

The velocity estimator is product-specific and therefore remains in
`velocity/estimation.py` rather than `calculations`. It processes `(frame, y,
x)` volumes in bounded temporal chunks, uses moments or calibrated frequency
bands as its active source, applies the shared local-background difference, and
returns physical velocity in `mm/s`. Moment normalization uses each frame's
spatial M0 mean, floored at float32 epsilon (`eps('single')`); NaNs propagate.
In band mode, finite negative samples in
both HF and LF are replaced with zero during chunk validation, with a warning
naming the dataset and frame interval. Zero LF denominators map to ratio zero
without division. Both bands must remain finite; validation and calibration provenance are part of the
pipeline contract.

Both velocity workflows use OpenCV Telea background inpainting with a fixed
5-pixel radius on independent float32 frames within the estimator's bounded
chunks. Method and radius remain recorded in velocity provenance.

`low_freq_spectra/spectrograms.py` reads HD frequency-bin spectra once per frame for
full-frame `Global/SpatialAverage`/`Global/SpatialMedian` reductions and optional gap median
power. `low_freq_spectra/preprocessing.py` interpolates Boolean gap support with
nearest neighbors to the spectrogram video's last two dimensions. The low-frequency spectral
runner publishes `/ProcessingSpectrograms/Gaps/SpatialMedian` with `(time, freq)`
axes. HD frames are flipped in spatial Y then transposed into EF's published
lower-left `(x, y)` mask frame at the canonical input boundary. Empty native or
interpolated support skips gap analysis with a warning; segmentation and other
pipeline work continue. Nonfinite included samples produce NaN per pair;
finite negative gap spatial medians are clamped to zero after reduction and before
beat subdivision. Gap provenance records `spatial_aggregation_method="median"`.
Individual raw pixels are not clipped. The hidden `low_freq_spectra`
pipeline runs once after velocity and passes detected moments timing to gap analysis,
falling back to detected band-ratio timing when moments timing is unavailable,
producing `/ProcessingSpectrograms/Gaps/PerBeat/Raw` with `(phase, freq, beat)` axes.
It uses 128 linear phase samples without bridging NaNs. Missing detected timing in both workflows
warns and skips per-beat analysis while other products remain available. Timing
workflow and HDF5 source paths are recorded in per-beat/filter/aggregate metadata.
`TimeAggregated/MedianPhaseBeat` pools phase and beats into a `(freq,)` spectrum;
`TimeAggregated/MedianBeat` reduces beats only, preserving `(phase, freq)` and
phase coordinates. Both ignore NaNs and keep all-missing output cells as NaN.
`TimeAggregated/MedianPhaseBeatEarly` and `MedianPhaseBeatLate` take phase
medians of `MedianBeat` over closed intervals `[0,1/3]` and `[1/3,1]`, yielding
`(freq,)` spectra. They share the phase-window helper with metric summaries,
include a one-third boundary sample in both windows, and retain NaN where a
frequency has no valid samples in that window. Preserve source timing, units,
frequency order, window endpoints, and selected/valid sample counts.
`TimeAggregated/StdBeat`, `StdBeatEarly`, and `StdBeatLate` take population STD
over beats after each beat's whole/early/late phase median, reusing `PerBeat/Raw`.
Ignore NaNs; one valid beat gives zero, none gives NaN. Persist raw power units,
`ddof=0`, reduction order, windows, and valid beat counts per frequency.
Positive-frequency PNG/EPS overlays use the actual positive/DC half of the
source: `P=floor(N/2)` bins at spacing `fs/N`. Normalize each curve by the
Hv, the same NaN-aware central-quarter median of the full spectrum used by A0
(64 source bins use bins 24 through 39), before cropping from stored `Metrics/index_max1`
through bin `P-1` inclusive. Scale STD by the same per-curve divisor, retain
linear power. The runner publishes the full signed frequency mapping once at
`/ProcessingSpectrograms/FrequencyBinTable`, as a `(freq,)` array in `kHz` where
each zero-based array index identifies the source bin. Use the
shared conversions in `calculations/math/frequency.py` for tables, metrics, and
plots. Table publication requires valid sampling calibration, independently of
gap support, cardiac timing, and figure export. Missing/invalid stored peaks,
all-NaN neighborhoods, or nonpositive divisors omit the affected figures.
Median/early/late use solid/dashed/dotted signals and envelope patterns.
The full-domain heatmap/median frequency axes use full-spectrum spacing `fs/N` with DC zero,
positive-then-negative order, and kHz labels. Resolve `fs` from shared
`Meta/SamplingFrequencyHz/value`; publish it from HD timing parameters/config
before export if absent, since velocity analysis may run later. Missing/invalid
sampling metadata skips frequency figures while keeping other products.
The spectrum figure and pooled scalar metrics consume `MedianPhaseBeat`.
All gap metrics live directly under `/ProcessingSpectrograms/Gaps/Metrics`.
Phase curves `A0_t`, `A1_t`, `A2_t`, `R1_t`, `R2_t`, `Z1_t`, and `Z2_t` use
`MedianBeat` with the pooled peak indices fixed at every phase. Pooled scalar
metrics use the same seven names without `_t`; output heights `h0`/`h1`/`hv`
and `rho_h` are absent.
Each curve also has `_early` and `_late` scalar medians over closed phase
intervals (default `[0,1/3]` and `[1/3,1]`), and a signed `delta_<metric>` scalar.
All seven deltas are late-minus-early. Generic helpers in
`calculations/math/phase_windows.py` own closed interval masks, NaN-aware medians,
and subtraction. `run_gap_analysis` and `pack_spectrum_metrics` accept
`early_interval`/`late_interval` cycle-fraction bounds in `[0,1]`; use the same
windows for early/late spectra, STD envelopes, and metric summaries/deltas.
Keep subtraction order explicit in definitions, source paths, and calculation
metadata. Ignore NaNs within each window; missing windows/differences remain NaN.
Keep the source's excluded endpoint and include a one-third boundary sample in
both default windows when present. Curves retain phase/timing/peak provenance and also
produce shared `gaps_<metric>` PNG/EPS plots.
The valley height `hv` uses a NaN-aware median over a neighborhood centered on
the two central bins. Its default length is one quarter of the actual domain
(32 bins -> 8 samples, 64 bins -> 16 samples), rounded to the nearest even
length with ties upward and at least two samples. Any one-dimensional domain
with at least two bins is supported; an explicit even `valley_neighborhood_length`
overrides the default. Peak halves and central valley bins scale with the domain,
splitting at `floor(bin_count/2)` for odd lengths. All-NaN peak halves retain NaN
indexes, log a warning, and are identified by per-half valid-bin counts in metadata.
`A0` is valley-median height times peak-to-peak width in kHz; `A2` is the
peak-height mean times the same width. These swap the former A0/A2 meanings.
Use `df_khz=fs/(N*1000)` and width `(index_max2-index_max1)*df_khz`, with no
extra bin. `A1` is the trapezoidal integral between those peak bins with positive
`df_khz` spacing in stored HD order (unwrapped kHz, no signed-axis sorting).
Areas and their phase summaries/deltas have source-power-times-kHz units;
peak/valley heights retain power units. Preserve calibration, spacing, integration
method, and units in metadata. Artifact axis labels abbreviate `a.u.*kHz` to
`a.u.`. Resolve calibration from shared
`Meta/SamplingFrequencyHz/value`; missing/invalid calibration warns and skips
the Metrics family while retaining raw/aggregate spectra and filters.
`R1=A1/A0`, `R2=A2/A0`,
`Z1=(A1-A0)/A0`, and `Z2=(A2-A0)/A0` are dimensionless; zero/missing A0 gives
NaN. All seven metrics have pooled scalars, phase curves, early/late medians,
signed late-minus-early deltas, and matching PNG/EPS phase plots.
Metric metadata records the neighborhood.
Independent `PerBeat/BandLimited` keeps DC plus the configured non-DC harmonics
(default 12) via the shared bandlimiting API in signed mode, using
periodic phase boundaries. Original durations define each effective Hz cutoff.
Incomplete beat/bin series remain unchanged with
a warning. `PerBeat/SVDContinuous` concatenates the raw beats as `(beat * phase, freq)`,
subtracts the temporal mean per bin, removes the first two modes, restores the
mean, and retains the 2-D `(time, freq)` continuous output. Sample spacing is
uniform within each beat but can differ between beats; retain source durations
and within-beat phase coordinates. All-invalid bins/incomplete rows
are excluded from the SVD fit and retained unchanged; insufficient dimensions
warn and leave the output unchanged. Preserve filtered undershoots to retain
DC power. Shared spectra are written after workflow cleanup and have no
`ProcessingBandRatio` duplicate. Full-frame spatial products remain available when
gap support is empty or missing. Settings/support/order travel in attributes.
`topology_core` owns native gap support and its segmentation outputs:
`Segmentation/Gaps/Legacy` preserves the vessel/disc/cutoff mask,
`GapsFromM0ff` thresholds bright vessels in the temporally median HD flat-field
moment after percentile stretch, CLAHE, and Gaussian smoothing, dilates the
vessels with the legacy 7x7 square, then inverts to white gaps. It has no
anatomical exclusions. `Mask` ORs their vessel/exclusion complements before
inverting, equivalent to AND of the saved white-gap masks, restricted to the image-centered
80%-width disc and excluding the optic disc (DV mask or centered 1.91 mm
fallback), used in shared state, spectral reductions, and the mask figure.
Missing moment0ff falls back to Legacy; missing Legacy permits the intensity
mask with the same disc constraints. Preserve native HD
axes, finite-sample temporal medians, and lower-left output serialization.
`low_freq_spectra` declares that dependency alongside `velocity_workflows`.
The pipeline automatically enters plans that schedule velocity, stays hidden from
selection, and skips absent/legacy-group spectra. Runtime source assembly keeps
spectra validation independent of velocity/topology; malformed spectra fail only
the spectral descriptor and preserve safe velocity analyses.
Physical capillary velocity remains deferred.

`spatial_gradient_moment0` shares cardiac-cycle/topology state but owns its
gradient profiles and lumen edges. Do not put gradient work into velocity
analysis.
It requires the optional HD flat-field moment at `/moment0ff` or `/M0FF`, and
publishes under `/ProcessingMoments/SpatialGradientProfiles` and
`/ProcessingMoments/SpatialGradientMetrics` only when it is a direct target. When it
runs solely as a blood-volume-rate dependency, keep the products in `ctx.state`
without persisting those intermediate HDF5 families. Its branch lumen plots
come from the computed `(time, branch)` metric, not a separate plot-only
calculation.
`blood_volume_rate` retrieves the products required by each enabled option and
must keep segment labels, centers, branch/radius dimensions, and pixel scale
aligned.

## Velocity semantics

All selected downstream targets run for both automatic velocity workflows. `doppler_moments`
and `frequency_bands` represent physical velocity (`mm/s`); band mode first
converts `HF / LF` to frequency using the recorded
`band_ratio_frequency_scale_hz`. `velocity/semantics.py`
centralizes display and dataset interpretation. Moments publish under
`/ProcessingMoments`, bands under `/ProcessingBandRatio`, with isolated downstream state and
`moments`/`bandratio` artifact folders. The single velocity pipeline declares
itself as an execution-variant producer and publishes one typed
`ExecutionVariant` per successful method; the engine fans out its declared DAG
dependents. Each workflow detects its own cycles from signed whole-vessel
physical velocity after estimation and uses them for its downstream analyses;
timing provenance records physical velocity in `mm/s`. Segmentation remains
under `/Segmentation`. New velocity-derived outputs or plots must resolve
semantics from payload/provenance instead of hard-coding `mm/s`.

Velocity is never dimensionless. Reject legacy unit-`1` or relative-index
velocity metadata instead of adapting labels or output units for it.

Historically physical downstream pipeline names are not scheduling guards.
Retain clear method, calibration, and unit provenance for their results.

## State and output rules

- Give shared state keys one owning runner and a declared producer DAG key.
- Use typed objects for multi-field state and validate type/shape in consumers.
- Put durable products in the HDF5 result, not only in `ctx.state` or figures.
- Use `EyeFlowOutputPaths.active()` for centralized output families.
- Pack data with units and `dimDesc`; method-dependent velocity datasets also
  need velocity semantics.
- Options should control only their described product families. Express their
  intra-pipeline prerequisites with `PipelineOption.requires` and upstream data
  with `PipelineOption.dag_requires`.
- Do not infer whether a sibling pipeline ran from target selection; use
  `ctx.pipeline_scheduled` or retrieve declared state.

## Change navigation

| Change | Read/inspect | Usually skip |
|---|---|---|
| New pipeline or option | `CONTRIBUTING.md`, package `__init__.py`, runner, DAG tests | GUI views; discovery is automatic |
| Velocity source/semantics | `vessel_inputs.py`, velocity sources/runner/estimator, data contracts | displacement internals |
| Per-beat/segment output | velocity-analysis builder and packers, output schema, profile/segment tests | settings UI |
| Spatial-gradient lumen metric | spatial-gradient package, vessel-segment sampling/measurement, image filters, spatial-gradient/topology tests | CLI and report code unless output is exposed there |
| Blood-volume-rate formula | option declaration, runner, outputs, `calculations/blood_volume_rate.py`, BVR tests | unrelated waveform metric calculators |
| Metric family | that pipeline's runner/calculator/outputs plus waveform input contracts and matching test file | pipeline engine unless dependencies/options change |
| Displacement | `displacement_map/` and `test_displacement_map_pipeline.py` | topology and waveform metrics unless the dormant segment helper is activated |
| Profile fit | `calculations/vessel_segments/profiles/fits/quadratic.py`, then `velocity_analysis/outputs/profile_analysis.py`, dedicated doc/test | GUI and settings |
