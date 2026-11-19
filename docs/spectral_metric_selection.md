# Spectral shape metrics

The hidden `low_freq_spectra` pipeline publishes exploratory vertical power
fractions (5%, 25%, 50%), Shannon entropy, absolute-frequency power quantiles
(including F50/F90), and early-to-late Wasserstein-1 distance. These are
descriptive metrics, with no diagnostic threshold or automatic sample rejection.
The output schema and exact names live in `src/input_output/schema/eyeflow_output.py`
and `src/pipelines/low_freq_spectra/metrics.py`.

## Selection evidence

A local cohort analysis on 2026-10-28 used 10 AMD, 21 OAG, 9 DR, and 28
control acquisitions. Repeated acquisitions were reduced to eye medians, then
participant medians. Participant and eye identities came from filenames and
need study-metadata confirmation. The apparent participant-level AUC against
controls was 1.0 for several selected AMD descriptors, but the same sample
selected and evaluated the metrics. Existing SEII_50 also separated AMD in
this sample. No new candidate passed the stated 0.05 threshold after
Benjamini-Hochberg correction of 45 candidates per disease comparison (or
135 jointly); the best adjusted AMD values were about 0.064 and 0.105,
respectively. External participant-level validation is required before any
diagnostic accuracy claim. The cohort archive and analysis script are not
included in this repository; this was a local exploratory analysis.

## Numerical definitions

Selected additions are stored directly under
`/ProcessingSpectrograms/Gaps/Metrics` as float32; calculations use float64.
They consume raw `MedianBeat` and `MedianPhaseBeat`, fixed pooled peak indexes,
and calibrated signed-bin frequencies. Missing calibration or usable peaks
leaves dependent metrics unavailable or NaN under the pipeline's optional
product rules.

For a selected peak interval with midpoint `m=floor(N/2)-0.5` and pooled peak
indexes `i,j`, `Av_p` integrates the two intervals
`[i, i+(m-i)*p/100]` and `[j-(j-m)*p/100, j]` for `p=5,25,50`.
Use trapezoids with interpolated fractional boundaries and positive calibrated
kHz spacing in stored bin order. `SEv_p=Av_p_t/A_t` requires finite positive
`A_t`. Whole-cycle `SEIv_p` is the median of phasewise ratios, not a ratio of
pooled areas. No valley baseline is subtracted. Areas have power·kHz units;
fractions are dimensionless.

Distribution descriptors use full discrete mass at every bin in the inclusive
fixed peak interval. `ValleyNormalizedSpectrum_t=P/min(P)` and
`SpectralProbability_t=P/sum(P)` are distinct arrays. Entropy is
`-sum(p*log2(p))` in bits. Quantile curves combine bins with equal absolute
frequency and use the discrete inverse CDF, so no sub-bin resolution is
invented. F50 and F90 select their columns. Whole/early/late summaries take
medians of phasewise values. `WassersteinEarlyLate` compares separately
normalized early and late phase-median spectra on absolute frequency; it is an
unsigned distance in kHz.

Early and late phase windows are closed `[0,1/3]` and `[1/3,1]` and may share
a boundary sample. Deltas are late minus early. Nonfinite or negative included
power invalidates that phase's probability descriptors. A zero valley makes
valley normalization undefined, while positive-total probabilities can remain
valid. Missing windows and missing pooled peaks produce NaN rather than
silently dropping bins. See `src/pipelines/low_freq_spectra/metrics.py` and
focused tests for exact supporting datasets, attributes, and figure selection.
