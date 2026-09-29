# Artery velocity-profile analysis

Selecting `velocity_profile_analysis` ensures the upstream waveform pipeline
computes segments and per-beat products and publishes velocity profiles, even
when those selectable outputs were disabled. Analysis reads the work HDF5 dataset:

`/Processing/VelocityProfiles/Artery/TransverseVelocityProfileMasked/value`

Input axes are `(x, time, beat, branch, radius)`. Every profile is fit independently;
there is no averaging across time, beats, branches or radii.

## Fit definition

Use `x = 0, ..., Nx-1` and minimize `sum(w * (v - (a*x*x+b*x+c))**2)` over
finite observed samples. For `u=x/(Nx-1)` and `d=abs(2*u-1)`, use `w=1-d^p`.
The default `p=2` gives a quadratic decrease from weight one at the domain
center to weight zero at either border. Direct API callers may supply any finite
`p>0`. With an even number of samples, the two samples nearest the center have
equal maximum sampled weights below one. Weights remain tied to the complete
input domain when some samples are missing. A one-sample domain has weight one.

The solver scales design rows and observations by `sqrt(w)` and works in
float64 on centered/scaled coordinates. Stored coefficients refer to the original
zero-based index coordinate. Time slabs contain at most 256 profiles per
beat/branch/radius, and profiles sharing finite masks reuse the least-squares
solve within each slab. The entire input dataset is not materialized.

## Output schema

Each name below is written to
`/Processing/VelocityProfileAnalysis/Artery/<name>/value` with axes
`(time, beat, branch, radius)`. Float arrays are float32 and counts are int32.

| Names | Definition |
|---|---|
| `a`, `b`, `c` | Quadratic coefficients; coefficient `b` is distinct from the input beat axis |
| `fit_rss` | Sum of squared residuals over finite samples |
| `fit_rmse` | `sqrt(fit_rss/n_fit_samples)` |
| `fit_r_squared` | `1-fit_rss/SST`, using the ordinary observed mean |
| `fit_weighted_rss` | Sum of `w*residual**2` |
| `fit_weighted_rmse` | `sqrt(fit_weighted_rss/sum(w))` |
| `fit_weighted_r_squared` | `1-fit_weighted_rss/weighted_SST`, using the weighted observed mean |
| `n_fit_samples` | Number of finite input samples |
| `index_center` | Fractional vertex index for a downward-opening, nondegenerate fit |
| `index_left_zero`, `index_right_zero` | Sorted fractional real roots, not clipped to the input range |
| `Qv_fit`, `Qv` | Unweighted signed sums of fitted and observed values on identical observed support |
| `n_area_samples` | Size of that observed support |

Observed support consists of finite input samples whose integer indexes are
inclusively between the roots. Indexes outside the sampled domain contribute
nothing. There is no extrapolation, quadrature, absolute value, or pixel-spacing
factor in either area sum. In particular, negative observed samples retain their
sign. Coefficients and area sums use index coordinates, not physical distances.

All datasets include axis, source path, zero-based-index, model, weighting, and
integration metadata. The `weight_power` attribute records the production value
`2.0`. Only artery analysis is produced.

## Invalid or degenerate profiles

- Fewer than three finite samples, rank deficiency or nonfinite coefficients:
  numerical outputs are NaN; the observed sample count is retained.
- Upward-opening or effectively linear fits retain coefficients and diagnostics,
  but geometry and areas are NaN. Negligible normalized curvature is at most
  `64*float64_epsilon*max(1,max(abs(v)))` in magnitude.
- Downward fits retain the vertex even without two distinct finite real roots.
  A discriminant within float64 roundoff of zero is treated as a repeated root;
  the tolerance is `64*epsilon*max(1,alpha**2,beta**2,abs(4*alpha*gamma))` for
  normalized coefficients `(alpha,beta,gamma)`.
- No usable roots or no observed indexes between roots: NaN areas and zero area count.
- Constant observations have undefined R-squared, returned as NaN. Unweighted
  R-squared may be negative because the fit optimizes weighted residuals.
- Empty time/beat/branch/radius axes are preserved. Missing input or a wrong
  dimensionality raises an explicit error before outputs are returned.
