# Exploratory selection of spectral shape metrics

The `low_freq_spectra` pipeline now calculates vertical power fractions at
5%, 25% and 50%, Shannon entropy, absolute-frequency power quantiles (including
F50/F90), and early-to-late Wasserstein-1 distance. These are exploratory
descriptors, with no diagnostic thresholds or sample rejection.

## Evidence used for selection

The local analysis used `doppler_manager_ef_cohorts(2).zip` on 2026-10-28.
Repeated acquisitions were first reduced to an eye median, then eyes to a
participant median, giving equal weight to each eye and each participant.
Participant and eye identities were inferred from filenames by removing the
acquisition date, acquisition suffix and eye suffix. These inferred identities
need confirmation against study metadata before inferential use.

| Cohort | Acquisitions | Eyes | Inferred participants |
|---|---:|---:|---:|
| AMD | 10 | 3 | 3 |
| OAG | 21 | 7 | 5 |
| DR | 9 | 4 | 3 |
| Controls | 28 | 9 | 5 |

The table reports descriptive participant-level AUC against controls, with
direction chosen from the observed data. Thus these are apparent separation
scores, not cross-validated diagnostic performance. AUC includes half credit
for ties. The AMD acquisition AUC column shows sensitivity to the aggregation
level; acquisitions are not independent observations.

| Metric | AMD AUC | AMD acquisition AUC | OAG AUC | DR AUC | AMD direction |
|---|---:|---:|---:|---:|---|
| SEIv_5 | 1.000 | 0.961 | 0.800 | 0.667 | Higher |
| SEIv_25 | 1.000 | 1.000 | 0.600 | 0.933 | Higher |
| SEIv_50 | 1.000 | 1.000 | 0.520 | 1.000 | Higher |
| SpectralEntropy | 1.000 | 1.000 | 0.600 | 0.933 | Lower |
| delta_SpectralEntropy | 1.000 | 1.000 | 0.680 | 0.933 | More negative |
| F50 | 1.000 | 0.991 | 0.620 | 0.967 | Lower |
| F90 | 1.000 | 0.982 | 0.560 | 0.967 | Lower |
| WassersteinEarlyLate | 1.000 | 1.000 | 0.640 | 0.867 | Higher |
| Existing SEII_50 | 1.000 | 0.996 | 0.560 | 1.000 | Higher |
| Existing SEI_50 | 0.867 | 0.864 | 0.720 | 0.667 | Higher |

AMD shows the clearest new signal: more power close to the two low-frequency
peaks, a less uniform spectrum, lower characteristic frequency magnitudes,
and greater early/late redistribution. This is a description of these data;
it does not establish a physiological cause. DR often shifts in the opposite
direction. OAG separation is weaker for the gain-invariant descriptors.

All of the selected AMD whole-cycle descriptors remain completely separated
when any single participant is omitted. This checks sensitivity to individual
participants; it is not independent validation. Existing SEII_50 already
separates AMD completely, so incremental predictive value has not been shown.

Exact two-sided permutation tests were calculated by enumerating participant
label assignments. Benjamini-Hochberg corrections covered 45 new scalar
candidates per disease comparison, and all 135 comparisons together.
The best AMD adjusted p-values were approximately 0.064 within AMD and 0.105
across comparisons. No new candidate met a 0.05 adjusted threshold.
Selection used these same observations, so external participant-level
validation is required before estimating diagnostic accuracy.

All acquisitions used 64 bins, FFT size 256 and sampling frequency 37037 Hz.
The included peak interval had 58 bins in AMD, OAG and controls. DR included
51- and 58-bin intervals, which can affect entropy's maximum and distribution
support. This does not explain AMD/control separation by support length, but
peak selection, preprocessing and acquisition/clinical covariates remain
possible influences.

SEIv_10 was omitted from regular runs because it was strongly redundant with
SEIv_5 and SEIv_25 (AMD/control participant Spearman correlation approximately
0.976). SEIv_5 retains a narrow peak measure and the better exploratory OAG
separation; SEIv_25 and SEIv_50 retain the stable AMD and broader DR signals.
Entropy, frequency quantiles and Wasserstein distance provide different
interpretations, though independence from SEIv has not been demonstrated.

Raw delta_Av_50 also separated OAG in this sample (participant AUC 1.000,
acquisition AUC 0.932). Its dependence on power gain makes it unsuitable as
a preferred cross-acquisition discriminator without additional calibration.
Av numerator outputs remain available as supporting measurements.

The local reproducibility script is `dev_scripts/analyze_spectral_cohorts.py`.
Its `cohort_metric_analysis/` outputs include filename grouping, participant
values, all rankings, correlations and a participant scatter plot. No archive
contents are changed by that analysis. Those development artifacts are local;
this document preserves the selection evidence in the repository.

## Regular execution and numerical definitions

The existing hidden `low_freq_spectra` pipeline produces these outputs after
the shared spectral preprocessing and beat alignment. CLI and GUI both use
the existing execution path. Missing calibration or beat timing follows the
existing optional-output policy. Calculations use raw `MedianBeat` and
`MedianPhaseBeat`, with fixed pooled peak indices, not filtered spectra.

All datasets live directly under `/ProcessingSpectrograms/Gaps/Metrics`.
The original 64 metrics retain their definitions. The selected additions give
63 more datasets, hence 127 when calibrated phase spectra are available.
Stored values are float32; integration and descriptor evaluation use float64.

For `p` in `{5,25,50}`, midpoint `m=floor(N/2)-0.5`, and pooled peak indices
`i,j`, integrate the intervals `[i,i+(m-i)*p/100]` and
`[j-(j-m)*p/100,j]`. Use trapezoids with linearly interpolated fractional
boundaries and positive calibrated kHz spacing in stored bin order.
`SEv_p=Av_p_t/A_t` requires finite, positive `A_t`. Whole-cycle `SEIv_p` is
the temporal median of that phasewise ratio, not a ratio of pooled areas.
No valley baseline is subtracted. At 50%, each interval occupies half its own
peak-to-midpoint span; their union does not fill the peak interval.

Each percentage adds `Av_p`, `Av_p_t`, `Av_p_early`, `Av_p_late`, `delta_Av_p`,
`SEv_p`, `SEIv_p`, `SEIv_p_early`, `SEIv_p_late`, and `delta_SEIv_p`.
Areas have source-power*kHz units; fractions are dimensionless.

Distribution descriptors use the inclusive fixed peak interval. Each HD band
contributes full discrete mass, including both endpoint bins; these are not
the trapezoidal endpoint weights used by Av and A.

- `ValleyMinimum_t` is the true minimum inside that interval, distinct from
  the original central-neighborhood median Hv.
- `ValleyNormalizedSpectrum_t=P/min(P)` has finite values at least 1.
  `SpectralProbability_t=P/sum(P)=R/sum(R)` is a separate probability array.
  Both have `(phase,freq)` axes and NaN outside the selected interval.
- `SpectralEntropy_t=-sum(p*log2(p))` is in bits on the original signed bins.
  Its maximum is log2(included bin count); the valley scale cancels.
- `PowerQuantileFrequency_t` has `(phase,power_percentile)` axes. Quantiles
  use absolute calibrated kHz coordinates, combining equal magnitudes and
  taking the discrete inverse CDF. `PowerQuantilePercent` is 0..100 in 1%
  increments. There is no invented sub-bin frequency resolution; 0% and 100%
  are the smallest and largest occupied magnitudes. `F50_t` and `F90_t`
  select the corresponding columns. `DistributionFrequencyMagnitude`
  records absolute bin coordinates in original order.
- Entropy, the quantile curve, F50 and F90 each also have a whole-cycle median,
  `_early`, `_late`, and `delta_` output. Median quantile curves summarize
  phasewise quantiles, rather than quantiles of a pooled spectrum.
- `WassersteinEarlyLate` compares separately normalized early/late
  phase-median spectra on absolute frequency. It is an unsigned distance in
  kHz: moving 20% of power by 5 kHz gives 1 kHz. Supporting
  `WassersteinEarlyProbability`, `WassersteinLateProbability`,
  `WassersteinEarlyCDF`, `WassersteinLateCDF`,
  `WassersteinEarlyQuantileFrequency`, `WassersteinLateQuantileFrequency`,
  and `WassersteinFrequencyMagnitude` describe the exact comparison.

Default early/late windows are closed `[0,1/3]` and `[1/3,1]`; shared boundary
samples belong to both. Summaries ignore invalid phases and use equal weight
per valid phase. Deltas are late minus early; missing windows give NaN.
Any nonfinite or negative included power invalidates that phase's probability
descriptors. A zero valley makes R undefined, while probabilities and entropy
remain defined if total power is positive. Missing pooled peaks leave
dependent measurements NaN. No bins are silently removed or shifted.

The standard artifact writer exports phase curves for Av/SEv at the selected
percentages, entropy, F50 and F90. It also exports
`gaps_power_quantiles` (whole/early/late median curves) and
`gaps_wasserstein_early_late` (early/late CDFs with their difference shaded).
Both PNG and EPS use the existing shared spectrogram artifact folders.
HDF5 metadata records source arrays, calibration, support, definitions,
normalization, phase windows and validity policies.
