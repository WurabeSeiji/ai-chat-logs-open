# ABC Multigauge Interference Readout Result v1

## Stage Structure

1. Reproduce the one-collision ABC elastic reflection map.
2. Reconstruct p-like and E-like quantities from complex interference ratios over multiple gauges.
3. Reconstruct R-like stable residuals from calibrated multigauge intensity readouts.
4. Check conservation, label preservation, compensated square closure, and t/R separation.

## Verdict

- baseline_collision_valid: `True`
- label_modes_preserved: `True`
- closure_preserved: `True`
- p_reconstructed_all_gauges: `True`
- E_reconstructed_all_gauges: `True`
- R_reconstructed_all_gauges: `True`
- p_reflection_valid: `True`
- E_preserved: `True`
- R_preserved: `True`
- R_gauge_stable: `True`
- t_R_separation_valid: `True`
- single_gauge_only_used: `False`
- multigauge_measurement_valid: `True`

## Key Numerical Values

- p_max_abs_error: `2.5202062658991053e-14`
- E_max_abs_error: `2.2315482794965646e-14`
- R_max_abs_error: `4.4408920985006262e-16`
- R_max_gauge_std: `2.2204460492503131e-16`
- closure_residual_abs: `0.0000000000000000e+00`
- separation_ratio_time: `3.6838616474030686e-30`

## Stage Readout Summary

| stage | particle | p_mean | p_std | E_mean | E_std | R_mean | R_std | t_mean |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| initial | A | 9.9999999999999800e-01 | 7.6855384341228018e-15 | 9.9999999999999800e-01 | 7.6855384341228018e-15 | 9.9999999999999989e-01 | 2.2204460492503131e-16 | -2.0000000000000001e-01 |
| pre_collision | A | 9.9999999999999800e-01 | 7.6855384341228018e-15 | 9.9999999999999800e-01 | 7.6855384341228018e-15 | 9.9999999999999989e-01 | 2.2204460492503131e-16 | -1.5999999999999848e-02 |
| collision_cell | A | 9.9999999999999800e-01 | 7.6855384341228018e-15 | 9.9999999999999800e-01 | 7.6855384341228018e-15 | 9.9999999999999989e-01 | 2.2204460492503131e-16 | -1.4999999999999847e-02 |
| collision_map | A | -9.9999999999999767e-01 | 8.6400570246037599e-15 | 9.9999999999999800e-01 | 7.6855384341228018e-15 | 9.9999999999999989e-01 | 2.2204460492503131e-16 | -1.4999999999999847e-02 |
| post_collision | A | -9.9999999999999767e-01 | 8.6400570246037599e-15 | 9.9999999999999800e-01 | 7.6855384341228018e-15 | 9.9999999999999989e-01 | 2.2204460492503131e-16 | -1.3999999999999846e-02 |
| final | A | -9.9999999999999767e-01 | 8.6400570246037599e-15 | 9.9999999999999800e-01 | 7.6855384341228018e-15 | 9.9999999999999989e-01 | 2.2204460492503131e-16 | 2.0000000000000032e-01 |
| initial | B | -9.9999999999999767e-01 | 8.6400570246037599e-15 | 9.9999999999999800e-01 | 7.6855384341228018e-15 | 9.9999999999999989e-01 | 2.2204460492503131e-16 | -2.0000000000000001e-01 |
| pre_collision | B | -9.9999999999999767e-01 | 8.6400570246037599e-15 | 9.9999999999999800e-01 | 7.6855384341228018e-15 | 9.9999999999999989e-01 | 2.2204460492503131e-16 | -1.5999999999999848e-02 |
| collision_cell | B | -9.9999999999999767e-01 | 8.6400570246037599e-15 | 9.9999999999999800e-01 | 7.6855384341228018e-15 | 9.9999999999999989e-01 | 2.2204460492503131e-16 | -1.4999999999999847e-02 |
| collision_map | B | 9.9999999999999800e-01 | 7.6855384341228018e-15 | 9.9999999999999800e-01 | 7.6855384341228018e-15 | 9.9999999999999989e-01 | 2.2204460492503131e-16 | -1.4999999999999847e-02 |
| post_collision | B | 9.9999999999999800e-01 | 7.6855384341228018e-15 | 9.9999999999999800e-01 | 7.6855384341228018e-15 | 9.9999999999999989e-01 | 2.2204460492503131e-16 | -1.3999999999999846e-02 |
| final | B | 9.9999999999999800e-01 | 7.6855384341228018e-15 | 9.9999999999999800e-01 | 7.6855384341228018e-15 | 9.9999999999999989e-01 | 2.2204460492503131e-16 | 2.0000000000000032e-01 |

## R Gain Readout Probe

| R_gain | R_mean | R_std | expected | abs error |
|---:|---:|---:|---:|---:|
| 1.0 | 1.0000000000000000e+00 | 2.5639502485114184e-16 | 1.0000000000000000e+00 | 0.0000000000000000e+00 |
| 10.0 | 1.0000000000000000e+01 | 2.0511601988091347e-15 | 1.0000000000000000e+01 | 0.0000000000000000e+00 |
| 100.0 | 1.0000000000000000e+02 | 2.4613922385709617e-14 | 1.0000000000000000e+02 | 0.0000000000000000e+00 |

## Files

| kind | file |
|---|---|
| JSON | `abc_multigauge_interference_readout_result_v2.json` |
| stage CSV | `abc_multigauge_interference_readout_timeline_v2.csv` |
| gauge CSV | `abc_multigauge_interference_readout_gauge_sweep_v2.csv` |
| summary CSV | `abc_multigauge_interference_readout_stage_summary_v2.csv` |
| event CSV | `abc_multigauge_interference_readout_events_v2.csv` |
| R gain CSV | `abc_multigauge_interference_readout_r_gain_sweep_v2.csv` |
| invariants plot | `abc_multigauge_interference_readout_invariants_v2.png` |
| t/R plot | `abc_multigauge_interference_readout_tr_separation_v2.png` |
