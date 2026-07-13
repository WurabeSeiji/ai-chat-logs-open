# ABC Multigauge Generalized Elastic Collision Velocity Sweep v1

## Purpose

This experiment tests the generalized R-weighted elastic collision map under asymmetric initial phase gradients.
It checks that the construction is not limited to the initial +1/-1 counter-propagating condition.

## Aggregate Verdict

- case_count: `9`
- collision_reached_all_cases: `True`
- individual_readout_valid_all_cases: `True`
- P_R_preserved_all_cases: `True`
- K_R_phase_preserved_all_cases: `True`
- relative_gradient_flipped_all_cases: `True`
- E_tau_R_preserved_all_cases: `True`
- R_total_preserved_all_cases: `True`
- max_p_abs_error: `4.1033842990145786e-13`
- max_E_abs_error: `2.6423307986078726e-14`
- max_R_abs_error: `1.865174681370263e-14`
- max_R_gauge_std: `8.078732199757252e-15`
- max_within_particle_separation_ratio_time: `2.365357211965411e-28`
- max_P_R_conservation_error: `2.8910207561239076e-13`
- max_K_R_phase_conservation_error: `1.5258905250448151e-12`
- max_relative_flip_error: `1.2434497875801753e-14`
- max_E_tau_R_conservation_error: `4.884981308350689e-15`
- max_R_total_conservation_error: `4.884981308350689e-15`
- single_gauge_only_used: `False`
- velocity_sweep_generalized_collision_valid: `True`

## Case Summary

| case | R_B/R_A | q before A/B | q after A/B | R*p err | R*p^2 err | relative flip err | valid |
|---|---:|---|---|---:|---:|---:|---|
| c01_A1.00_B1.00_u1.00_v-1.00 | 1.0000000000000000e+00 | 1 / -1 | -1 / 1 | 0.0000000000000000e+00 | 0.0000000000000000e+00 | 0.0000000000000000e+00 | `True` |
| c02_A1.00_B1.00_u1.40_v-0.60 | 1.0000000000000000e+00 | 1.4 / -0.6 | -0.6 / 1.4 | 0.0000000000000000e+00 | 0.0000000000000000e+00 | 0.0000000000000000e+00 | `True` |
| c03_A1.00_B1.00_u0.50_v-1.70 | 1.0000000000000000e+00 | 0.5 / -1.7 | -1.7 / 0.5 | 1.8052226380405045e-13 | 6.2705396430828841e-13 | 0.0000000000000000e+00 | `True` |
| c04_A1.00_B2.00_u1.00_v-1.00 | 4.0000000000000000e+00 | 1 / -1 | -2.2 / -0.2 | 1.2878587085651816e-14 | 2.3980817331903381e-14 | 1.2434497875801753e-14 | `True` |
| c05_A1.00_B2.00_u1.40_v-0.60 | 4.0000000000000000e+00 | 1.4 / -0.6 | -1.8 / 0.2 | 1.8651746813702630e-14 | 1.4210854715202004e-14 | 2.2204460492503131e-15 | `True` |
| c06_A1.00_B2.00_u1.20_v0.20 | 4.0000000000000000e+00 | 1.2 / 0.2 | -0.4 / 0.6 | 2.4424906541753444e-14 | 1.5543122344752192e-15 | 1.1324274851176597e-14 | `True` |
| c07_A2.00_B1.00_u0.80_v-1.50 | 2.5000000000000000e-01 | 0.8 / -1.5 | -0.12 / 2.18 | 5.2846615972157451e-14 | 8.8817841970012523e-16 | 1.0214051826551440e-14 | `True` |
| c08_A1.50_B1.00_u1.80_v-0.20 | 4.4444444444444442e-01 | 1.8 / -0.2 | 0.569231 / 2.56923 | 2.8910207561239076e-13 | 1.5258905250448151e-12 | 5.3290705182007514e-15 | `True` |
| c09_A1.00_B3.00_u0.80_v-0.40 | 9.0000000000000000e+00 | 0.8 / -0.4 | -1.36 / -0.16 | 7.9936057773011271e-14 | 3.8191672047105385e-14 | 4.4408920985006262e-16 | `True` |

## Files

| kind | file |
|---|---|
| JSON | `abc_multigauge_generalized_elastic_collision_velocity_sweep_result_v2.json` |
| case CSV | `abc_multigauge_generalized_elastic_collision_velocity_sweep_cases_v2.csv` |
| stage quantity CSV | `abc_multigauge_generalized_elastic_collision_velocity_sweep_stage_quantities_v2.csv` |
| gauge CSV | `abc_multigauge_generalized_elastic_collision_velocity_sweep_gauge_rows_v2.csv` |
| error plot | `abc_multigauge_generalized_elastic_collision_velocity_sweep_errors_v2.png` |
| q output plot | `abc_multigauge_generalized_elastic_collision_velocity_sweep_q_outputs_v2.png` |
