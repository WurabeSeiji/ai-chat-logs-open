# ABC Multigauge Generalized Elastic Collision Extreme R Sweep v1

## Purpose

This experiment sweeps the readout R ratio over extreme asymmetric amplitude conditions.
It checks whether the generalized R-weighted elastic map and multigauge readout survive large R contrast.

## Aggregate Verdict

- case_count: `12`
- all_cases_valid: `True`
- collision_reached_all_cases: `True`
- individual_readout_valid_all_cases: `True`
- P_R_preserved_all_cases: `True`
- K_R_phase_preserved_all_cases: `True`
- relative_gradient_flipped_all_cases: `True`
- E_tau_R_preserved_all_cases: `True`
- R_total_preserved_all_cases: `True`
- max_R_dynamic_range: `64.0`
- min_R_ratio_B_over_A: `0.015625`
- max_R_ratio_B_over_A: `64.0`
- max_p_abs_error: `1.000310945187266e-13`
- max_E_abs_error: `2.2315482794965646e-14`
- max_R_abs_error: `2.842170943040401e-14`
- max_R_gauge_std: `1.7404671430534633e-14`
- max_within_particle_separation_ratio_time: `7.143466159887607e-30`
- max_P_R_conservation_error: `6.465938895416912e-13`
- max_K_R_phase_conservation_error: `1.2789769243681803e-12`
- max_relative_flip_error: `2.6867397195928788e-14`
- max_E_tau_R_conservation_error: `8.881784197001252e-16`
- max_R_total_conservation_error: `3.552713678800501e-15`
- single_gauge_only_used: `False`
- extreme_R_sweep_valid: `True`

## Case Summary

| case | R_B/R_A | dynamic range | q after A/B | R*p err | R*p^2 err | valid |
|---|---:|---:|---|---:|---:|---|
| c01_A1.000_B0.125_u1.00_v-0.50 | 1.5625000000000000e-02 | 6.4000000000000000e+01 | 0.95384615 / 2.4538462 | 1.0103029524088925e-14 | 1.9984014443252818e-14 | `True` |
| c02_A1.000_B0.250_u1.00_v-0.50 | 6.2500000000000000e-02 | 1.6000000000000000e+01 | 0.82352941 / 2.3235294 | 1.0436096431476471e-14 | 1.4432899320127035e-14 | `True` |
| c03_A1.000_B0.500_u1.00_v-0.50 | 2.5000000000000000e-01 | 4.0000000000000000e+00 | 0.4 / 1.9 | 5.4400928206632670e-15 | 8.4376949871511897e-15 | `True` |
| c04_A1.000_B1.000_u1.00_v-0.50 | 1.0000000000000000e+00 | 1.0000000000000000e+00 | -0.5 / 1 | 0.0000000000000000e+00 | 0.0000000000000000e+00 | `True` |
| c05_A1.000_B2.000_u1.00_v-0.50 | 4.0000000000000000e+00 | 4.0000000000000000e+00 | -1.4 / 0.1 | 4.5519144009631418e-15 | 2.9531932455029164e-14 | `True` |
| c06_A1.000_B4.000_u1.00_v-0.50 | 1.6000000000000000e+01 | 1.6000000000000000e+01 | -1.8235294 / -0.32352941 | 2.3714363805993344e-13 | 2.0161650127192843e-13 | `True` |
| c07_A1.000_B8.000_u1.00_v-0.50 | 6.4000000000000000e+01 | 6.4000000000000000e+01 | -1.9538462 / -0.45384615 | 1.3145040611561853e-13 | 1.9895196601282805e-13 | `True` |
| c08_A0.250_B1.000_u1.00_v-0.50 | 1.6000000000000000e+01 | 1.6000000000000000e+01 | -1.8235294 / -0.32352941 | 1.4821477378745840e-14 | 1.2601031329495527e-14 | `True` |
| c09_A0.500_B1.000_u1.00_v-0.50 | 4.0000000000000000e+00 | 4.0000000000000000e+00 | -1.4 / 0.1 | 1.1379786002407855e-15 | 7.3829831137572910e-15 | `True` |
| c10_A2.000_B1.000_u1.00_v-0.50 | 2.5000000000000000e-01 | 4.0000000000000000e+00 | 0.4 / 1.9 | 2.1760371282653068e-14 | 3.3750779948604759e-14 | `True` |
| c11_A4.000_B1.000_u1.00_v-0.50 | 6.2500000000000000e-02 | 1.6000000000000000e+01 | 0.82352941 / 2.3235294 | 1.6697754290362354e-13 | 2.3092638912203256e-13 | `True` |
| c12_A8.000_B1.000_u1.00_v-0.50 | 1.5625000000000000e-02 | 6.4000000000000000e+01 | 0.95384615 / 2.4538462 | 6.4659388954169117e-13 | 1.2789769243681803e-12 | `True` |

## Files

| kind | file |
|---|---|
| JSON | `abc_multigauge_generalized_elastic_collision_extreme_R_sweep_result_v1.json` |
| case CSV | `abc_multigauge_generalized_elastic_collision_extreme_R_sweep_cases_v1.csv` |
| stage quantity CSV | `abc_multigauge_generalized_elastic_collision_extreme_R_sweep_stage_quantities_v1.csv` |
| gauge CSV | `abc_multigauge_generalized_elastic_collision_extreme_R_sweep_gauge_rows_v1.csv` |
| error plot | `abc_multigauge_generalized_elastic_collision_extreme_R_sweep_errors_v1.png` |
| output plot | `abc_multigauge_generalized_elastic_collision_extreme_R_sweep_outputs_v1.png` |
