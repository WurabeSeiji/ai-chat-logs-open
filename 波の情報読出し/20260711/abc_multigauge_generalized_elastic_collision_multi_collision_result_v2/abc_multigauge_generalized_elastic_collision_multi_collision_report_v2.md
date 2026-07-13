# ABC Multigauge Generalized Elastic Collision Multi Collision v1

## Purpose

This experiment repeats the generalized R-weighted elastic collision map across multiple AB collisions with wall returns.
It checks whether R*p, R*p^2, relative phase-gradient flip, R*E_tau, and R stability survive repeated collision cycles.

## Aggregate Verdict

- case_count: `4`
- all_cases_valid: `True`
- completed_target_collisions_all_cases: `True`
- individual_readout_valid_all_cases: `True`
- closure_preserved_all_cases: `True`
- P_R_preserved_each_collision_all_cases: `True`
- K_R_preserved_each_collision_all_cases: `True`
- relative_flip_each_collision_all_cases: `True`
- E_tau_R_preserved_each_collision_all_cases: `True`
- R_preserved_each_collision_all_cases: `True`
- max_ab_collision_count: `6`
- max_wall_reflection_count: `8`
- max_p_abs_error: `1.0336176359260207e-13`
- max_E_abs_error: `3.341771304121721e-13`
- max_R_abs_error: `2.255973186038318e-13`
- max_R_gauge_std: `9.768953270904304e-14`
- max_within_particle_separation_ratio_time: `1.4212404920023786e-27`
- max_P_R_error: `3.552713678800501e-14`
- max_K_R_error: `1.056932319443149e-13`
- max_relative_flip_error: `2.220446049250313e-14`
- max_E_tau_R_error: `0.0`
- max_R_A_error: `4.440892098500626e-16`
- max_R_B_error: `8.881784197001252e-16`
- single_gauge_only_used: `False`
- generalized_multi_collision_valid: `True`

## Case Summary

| case | AB collisions | wall reflections | max R*p err | max R*p^2 err | max relative err | valid |
|---|---:|---:|---:|---:|---:|---|
| c01_A1.00_B1.00_u1.00_v-1.00 | 6 | 5 | 0.0000000000000000e+00 | 0.0000000000000000e+00 | 0.0000000000000000e+00 | `True` |
| c02_A1.00_B2.00_u1.40_v-0.60 | 6 | 6 | 3.5527136788005009e-14 | 7.9936057773011271e-14 | 1.6875389974302379e-14 | `True` |
| c03_A1.00_B2.00_u1.20_v0.20 | 6 | 8 | 2.4646951146678475e-14 | 3.5527136788005009e-14 | 1.7097434579227411e-14 | `True` |
| c04_A1.50_B1.00_u1.80_v-0.20 | 6 | 6 | 3.3750779948604759e-14 | 1.0569323194431490e-13 | 2.2204460492503131e-14 | `True` |

## Files

| kind | file |
|---|---|
| JSON | `abc_multigauge_generalized_elastic_collision_multi_collision_result_v2.json` |
| case CSV | `abc_multigauge_generalized_elastic_collision_multi_collision_cases_v2.csv` |
| collision CSV | `abc_multigauge_generalized_elastic_collision_multi_collision_readouts_v2.csv` |
| gauge CSV | `abc_multigauge_generalized_elastic_collision_multi_collision_gauge_rows_v2.csv` |
| P error plot | `abc_multigauge_generalized_elastic_collision_multi_collision_p_errors_v2.png` |
| summary plot | `abc_multigauge_generalized_elastic_collision_multi_collision_summary_v2.png` |
