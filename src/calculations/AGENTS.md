# Calculation map

`calculations` owns reusable numerical operations. Keep it independent of
Tkinter, settings persistence, pipeline selection, `.holo` layout, and result
paths. Product-specific orchestration and output packing belong in the
pipeline package.

| Task | Start with | Typical consumer |
|---|---|---|
| Cardiac timing, per-beat signals, waveform metrics | `blood_flow_velocity/signal_analysis/` | velocity and waveform pipelines |
| Native vessel geometry, branch identity, optic disc, gap zones | `topology/` | `pipelines/topology_core/` |
| Segment sampling and transforms | `vessel_segments/sampling/` | velocity and spatial-gradient analyses |
| Segment measurement and profiles | `vessel_segments/measurement/`, `vessel_segments/profiles/` | velocity analysis |
| Weighted quadratic profile fit | `vessel_segments/profiles/fits/quadratic.py` | `pipelines/velocity_analysis/outputs/profile_analysis.py` |
| Shared filtering, frequency, phase windows | `math/` | velocity and spectral pipelines |
| Blood-volume-rate formulas | `blood_volume_rate.py` | `pipelines/blood_volume_rate/` |
| CPU/GPU selection | `compute_backend.py` | sampling and image calculations |

Native topology defines spatial geometry; segment sampling and measurement
consume it. Retinal velocity estimation itself lives in
`pipelines/velocity/estimation.py`. For the profile fit's source and output
contracts, read [velocity-profile analysis](../../docs/velocity_profile_analysis.md).

Treat axes, NaN handling, sign, units, index base, and boundary rules as part of
each numerical API. Search callers and focused numerical-reference tests before
changing one. Functions should not mutate caller arrays unless documented.
