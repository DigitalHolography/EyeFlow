# Spatial-gradient preprocessing experiments

All experimental settings live in `preprocessing_config.json` beside this
file. The hidden `spatial_gradient_moment0` pipeline reads that file on every
EyeFlow run. To use a configuration stored elsewhere, set
`EYEFLOW_SPATIAL_GRADIENT_CONFIG` to its path.

The active `legacy_existing` preset preserves the repository's reference
behavior: existing geometry crop and interpolation/rotation, temporal median,
Sobel magnitude, temporal Gaussian, then the existing transverse profile and
peak detector. It is intentionally retained as a control. Change only
`active_preset` to select one of the A–J starting points.

## Controlled overrides

`overrides` is merged after the selected preset. It is the quickest place to
change one parameter without copying a preset:

```json
"overrides": {
  "gaussian3d": {
    "sigma_x": 1.0,
    "sigma_y": 3.0,
    "sigma_t": 2.0
  },
  "gradient": {
    "direction": "x",
    "signed": true
  }
}
```

Stage order comes only from `pipeline`. For the central interpolation test,
compare these orders while keeping every parameter equal:

```json
["crop", "gaussian3d", "gradient", "rotate"]
["crop", "gaussian3d", "rotate", "gradient"]
```

`rotate` contains the repository's existing resize, center padding, and
rotation operations. They are kept together so no extra resampling is added.
The shared implementation currently supports its baseline bilinear
interpolation only; unsupported interpolation names fail explicitly.

## Sweeps

Set `sweep.enabled` to `true`. `presets` may select several named starting
points, and `parameters` defines a Cartesian product of dotted configuration
paths. The checked-in disabled sweep already lists the complete initial A–J
set, including median/mean windows 5, 9, and 15 and derivative-Gaussian scales
1, 1.5, and 2. Replace that list when defining a focused Cartesian sweep:

```json
"sweep": {
  "enabled": true,
  "presets": ["experiment_j_gradient_before_rotation"],
  "parameters": {
    "gaussian3d.sigma_x": [0.5, 1.0, 1.5, 2.0],
    "gaussian3d.sigma_y": [1.0, 2.0, 3.0, 4.0],
    "gaussian3d.sigma_t": [1.0, 2.0, 3.0, 4.0],
    "gradient.method": ["sobel", "scharr", "derivative_gaussian"]
  }
}
```

For pipeline-order sweeps, a list is one candidate value, so `pipeline` takes
a list of lists. Named J presets are usually clearer for this comparison.

Configured sweeps process all valid artery and vein segments because the
verification target is the existing masked lumen-size result for every branch.
This is substantially more expensive than the active preset alone. The
selected active preset still feeds the normal EyeFlow HDF5 outputs; sweep
variants remain sidecar-only and do not create one HDF5 per configuration.

## Outputs and interpretation

Set `output.save_intermediates` or `output.save_final_float_tiff` to enable
active-run artifacts. Each sweep experiment gets a non-overwriting directory under
`spatial_gradient_experiments` containing:

- the exact resolved `configuration.json`;
- `artery_masked_lumen_size_time_branch.npz` and the equivalent vein file;
- one `masked_lumen_size_all_branches.png` graph containing every artery and
  vein branch;
- one appended row in `lumen_size_results.csv` with flattened settings and
  per-vessel validity, mean, standard deviation, and branch count.

The scientific source is the unchanged
`Masked/tbkr/lumen/size` result with dimensions `(time, beat, branch, radius)`.
Verification applies a NaN-aware arithmetic mean over beat and radius only,
leaving `(time, branch)`. Time is the normalized within-beat phase created by
the existing per-beat interpolation. All branches are drawn on the same axes;
arteries use solid lines and veins use dashed lines.

Scientific arrays remain float32 and missing results remain NaN. Lumen size is
still produced by the existing masked transverse profile, unchanged two-peak
detector, and right-minus-left peak separation. The sweep only changes how
that existing result is reduced and visualized for verification.

NLM and bilateral filtering are optional and disabled by default because they
are comparatively expensive. Anisotropic diffusion was not added: the current
dependency set has no dedicated, reliable implementation, and adding a fragile
dependency solely for that experiment would work against reproducibility.
