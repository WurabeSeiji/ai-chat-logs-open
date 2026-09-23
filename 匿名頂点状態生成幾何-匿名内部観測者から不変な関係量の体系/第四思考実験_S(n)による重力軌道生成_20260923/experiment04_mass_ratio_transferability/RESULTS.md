| q | nu | radial RMS/a | PN precession | S(n) precession | precession residual | PN peri shrink | S(n) peri shrink |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0.25000000 | 0.5084% | 0.331147 | 0.340531 | 0.009385 | 1.017150e-03 | 9.335690e-04 |
| 2 | 0.22222222 | 0.4967% | 0.332513 | 0.340496 | 0.007983 | 9.040754e-04 | 8.296700e-04 |
| 4 | 0.16000000 | 0.4704% | 0.331694 | 0.340419 | 0.008725 | 6.500540e-04 | 5.970898e-04 |
| 10 | 0.08264463 | 0.4379% | 0.333864 | 0.340322 | 0.006458 | 3.353587e-04 | 3.082393e-04 |

The radial orbit mismatch decreases from about 0.508% at q=1 to about 0.438% at q=10. The same S(n) rule therefore remains stable under the three new mass ratios without case-specific fitting. This is consistent with the conservative S(n) angular rule being Schwarzschild/test-particle-like: as q increases, nu decreases and the finite-mass-ratio PN reference moves toward the test-particle regime. This is an interpretation of the trend, not a derivation.

All cases preserve `det S = 1` to floating-point precision (max deviation <= ~5.6e-16).

## Raw data

For every q there are three pre-plot row files in `raw/`: `pn_integrator_*.csv` (raw PN RK4 rows), `pn_plot_*.csv` (uniform-time PN rows used for plotting), and `sn_*.csv` (every S(n) step including p,e and local radiation-reaction diagnostics). Therefore additional metrics can be computed without rerunning the dynamics.

## Reproduction

Run `python code/run_all_mass_ratio.py`. The plotting program reads the saved CSV files and does not recompute the dynamics.
