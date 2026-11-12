# Velocity-profile quadratic analysis

The `velocity_analysis.velocity_profile_analysis` option requires velocity
profiles and segments. It fits each profile independently from the artery and
vein `VelocityProfiles/.../Transversal/Masked/VelocityProfile/value` families
under `/ProcessingMoments` or `/ProcessingBandRatio`. Input axes are
`(x, time, beat, branch, radius)`. The reusable solver is
`src/calculations/vessel_segments/profiles/fits/quadratic.py`; the pipeline
adapter selects source paths and writes the schema.

## Fit and outputs

For zero-based sample index `x`, minimize
`sum(w * (v - (a*x*x + b*x + c))**2)` over finite observations. With
`u=x/(Nx-1)` and `d=abs(2*u-1)`, use `w=1-d**p`; production uses `p=2`.
Weights are tied to the complete input domain even when samples are missing.
The one-sample domain has weight one. The solver uses float64 internally and
stores coefficients for the original zero-based index coordinate.

Each vessel has output
`/Processing{Moments|BandRatio}/VelocityProfileAnalysis/{Artery|Vein}/<name>/value`
with `(time, beat, branch, radius)` axes. Float outputs are float32; counts are
int32.

| Names | Meaning |
|---|---|
| `a`, `b`, `c` | Quadratic coefficients |
| `fit_rss`, `fit_rmse`, `fit_r_squared` | Ordinary residual diagnostics |
| `fit_weighted_rss`, `fit_weighted_rmse`, `fit_weighted_r_squared` | Weighted diagnostics |
| `n_fit_samples` | Number of finite observations |
| `index_center`, `index_left_zero`, `index_right_zero` | Vertex and sorted real roots of a downward-opening fit |
| `Qv_fit`, `Qv`, `n_area_samples` | Signed fitted and observed sums and support count |

`Qv_fit` and `Qv` sum over the same observed finite indexes inclusively between
the roots and inside the sampled domain. They use index coordinates, with no
extrapolation, absolute value, quadrature, or pixel-spacing factor. Negative
observations keep their sign. The model is a weighted polynomial fit, without
a Poiseuille-flow constraint.

## Invalid profiles

Fewer than three finite samples, rank deficiency, or nonfinite coefficients
give NaN numerical outputs while preserving the sample count. Upward-opening
or effectively linear fits retain coefficients and diagnostics but have NaN
geometry and areas. Downward fits may retain a vertex without usable roots;
then areas are NaN and the area count is zero. Constant observations have
undefined R-squared. Missing input or wrong dimensionality raises an error;
empty non-profile axes are preserved. For curvature/discriminant tolerances and
solver chunking, use the solver and its numerical-reference tests.
