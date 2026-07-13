# ABC Multigauge Generalized Elastic Collision Readout v1

## Purpose

This experiment replaces the equal-amplitude q-flip map with a generalized 1D elastic map using the multigauge R readout as the mass-like weight.
It checks conservation of R*p and R*p^2 across asymmetric amplitude cases.

## Aggregate Verdict

- case_count: `8`
- individual_readout_valid_all_cases: `True`
- generalized_P_R_preserved_all_cases: `True`
- generalized_K_R_phase_preserved_all_cases: `True`
- E_tau_R_preserved_all_cases: `True`
- R_total_preserved_all_cases: `True`
- max_p_abs_error: `2.8688162956314045e-13`
- max_E_abs_error: `2.6423307986078726e-14`
- max_R_abs_error: `1.865174681370263e-14`
- max_R_gauge_std: `8.078732199757252e-15`
- max_within_particle_separation_ratio_time: `2.365357211965411e-28`
- max_P_R_conservation_error: `2.3803181647963356e-13`
- max_K_R_phase_conservation_error: `1.4086509736443986e-12`
- max_E_tau_R_conservation_error: `3.552713678800501e-15`
- max_R_total_conservation_error: `5.329070518200751e-15`
- single_gauge_only_used: `False`
- generalized_elastic_collision_readout_valid: `True`

## Case Summary

| case | R_B/R_A | q_A after | q_B after | R*p err | R*p^2 err | R*E_tau err | simple q-flip R*p err | valid readout |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| A_1.00_B_1.00 | 1.0000000000000000e+00 | -1.0000000000000000e+00 | 1.0000000000000000e+00 | 0.0000000000000000e+00 | 0.0000000000000000e+00 | 0.0000000000000000e+00 | 0.0000000000000000e+00 | `True` |
| A_1.00_B_1.10 | 1.2100000000000002e+00 | -1.1900452488687785e+00 | 8.0995475113122162e-01 | 1.2212453270876722e-14 | 1.4654943925052066e-14 | 0.0000000000000000e+00 | 4.2000000000000037e-01 | `True` |
| A_1.00_B_1.25 | 1.5625000000000000e+00 | -1.4390243902439024e+00 | 5.6097560975609762e-01 | 6.7723604502134549e-15 | 1.8651746813702630e-14 | 4.4408920985006262e-16 | 1.1250000000000000e+00 | `True` |
| A_1.00_B_1.50 | 2.2500000000000000e+00 | -1.7692307692307692e+00 | 2.3076923076923078e-01 | 5.5511151231257827e-15 | 1.3322676295501878e-15 | 0.0000000000000000e+00 | 2.5000000000000000e+00 | `True` |
| A_1.00_B_2.00 | 4.0000000000000000e+00 | -2.2000000000000002e+00 | -1.9999999999999996e-01 | 1.2878587085651816e-14 | 2.3980817331903381e-14 | 0.0000000000000000e+00 | 6.0000000000000000e+00 | `True` |
| A_1.00_B_3.00 | 9.0000000000000000e+00 | -2.6000000000000001e+00 | -6.0000000000000009e-01 | 2.3803181647963356e-13 | 1.4086509736443986e-12 | 3.5527136788005009e-15 | 1.6000000000000000e+01 | `True` |
| A_1.50_B_1.00 | 4.4444444444444442e-01 | -2.3076923076923078e-01 | 1.7692307692307692e+00 | 1.0436096431476471e-14 | 1.4654943925052066e-14 | 0.0000000000000000e+00 | 2.5000000000000000e+00 | `True` |
| A_2.00_B_1.00 | 2.5000000000000000e-01 | 1.9999999999999996e-01 | 2.2000000000000002e+00 | 9.7699626167013776e-15 | 6.2172489379008766e-15 | 0.0000000000000000e+00 | 6.0000000000000000e+00 | `True` |

## Interpretation

The generalized map preserves the R-weighted phase-gradient momentum and R-weighted phase-gradient square across every tested amplitude case.
The simple q-flip map appears as the equal-R special case and fails for unequal R in the R*p readout.

## Files

| kind | file |
|---|---|
| JSON | `abc_multigauge_generalized_elastic_collision_readout_result_v2.json` |
| case CSV | `abc_multigauge_generalized_elastic_collision_cases_v2.csv` |
| stage quantity CSV | `abc_multigauge_generalized_elastic_collision_stage_quantities_v2.csv` |
| gauge CSV | `abc_multigauge_generalized_elastic_collision_gauge_rows_v2.csv` |
| conservation plot | `abc_multigauge_generalized_elastic_collision_conservation_v2.png` |
| q output plot | `abc_multigauge_generalized_elastic_collision_q_outputs_v2.png` |
