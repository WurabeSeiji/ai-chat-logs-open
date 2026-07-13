# ABC Multigauge Generalized Elastic Collision Noise Robustness v1

## Purpose

This experiment injects deterministic readout-side noise into p/E/R gauge rows after the physical state has been simulated.
Zero-mean gauge noise is expected to cancel by multigauge averaging; common bias is expected to remain detectable.

## Aggregate Verdict

- case_count: `4`
- noise_mode_count: `2`
- noise_level_count: `6`
- total_summary_rows: `48`
- max_gauge_count: `116`
- zero_mean_multigauge_valid_all: `True`
- common_bias_detection_floor: `1e-10`
- common_bias_detected_all_above_floor: `True`
- zero_mean_max_p_mean_abs_error: `1.4477308241112041e-13`
- zero_mean_max_E_mean_abs_error: `6.439293542825908e-15`
- zero_mean_max_R_mean_abs_error: `3.552713678800501e-14`
- zero_mean_max_P_R_error: `1.674216321134736e-13`
- zero_mean_max_K_R_error: `8.322231792590173e-13`
- zero_mean_max_relative_flip_error: `1.687538997430238e-14`
- biased_control_max_p_mean_abs_error: `0.00021624903162420495`
- biased_control_max_E_mean_abs_error: `9.971756487914263e-05`
- biased_control_max_R_mean_abs_error: `0.00039675791940574356`
- single_gauge_only_used: `False`
- noise_robustness_valid: `True`

## Zero-Mean Summary

| noise level | max p mean err | max R mean err | max R*p^2 err | valid all cases |
|---:|---:|---:|---:|---|
| 0.0e+00 | 1.4432899320127035e-13 | 3.5527136788005009e-14 | 8.2867046558021684e-13 | `True` |
| 1.0e-12 | 1.4477308241112041e-13 | 3.5527136788005009e-14 | 8.2867046558021684e-13 | `True` |
| 1.0e-10 | 1.4477308241112041e-13 | 3.5527136788005009e-14 | 8.3044682241961709e-13 | `True` |
| 1.0e-08 | 1.4477308241112041e-13 | 3.5527136788005009e-14 | 8.3133500083931722e-13 | `True` |
| 1.0e-06 | 1.4477308241112041e-13 | 3.5527136788005009e-14 | 8.3222317925901734e-13 | `True` |
| 1.0e-04 | 1.4432899320127035e-13 | 3.5527136788005009e-14 | 8.2778228716051672e-13 | `True` |

## Common Bias Control

| noise level | max p mean err | max R mean err | detected all cases |
|---:|---:|---:|---|
| 1.0e-12 | 2.1973534103381098e-12 | 3.7871927816013340e-12 | `False` |
| 1.0e-10 | 1.8466961293484019e-10 | 3.8134029267666847e-10 | `True` |
| 1.0e-08 | 2.3795240000623608e-08 | 3.9659822093085495e-08 | `True` |
| 1.0e-06 | 1.7035717845281795e-06 | 3.9197052124073650e-06 | `True` |
| 1.0e-04 | 2.1624903162420495e-04 | 3.9675791940574356e-04 | `True` |

## Files

| kind | file |
|---|---|
| JSON | `abc_multigauge_generalized_elastic_collision_noise_robustness_result_v2.json` |
| summary CSV | `abc_multigauge_generalized_elastic_collision_noise_robustness_summary_v2.csv` |
| quantity CSV | `abc_multigauge_generalized_elastic_collision_noise_robustness_stage_quantities_v2.csv` |
| gauge CSV | `abc_multigauge_generalized_elastic_collision_noise_robustness_gauge_rows_v2.csv` |
| zero mean plot | `abc_multigauge_generalized_elastic_collision_noise_robustness_zero_mean_v2.png` |
| bias plot | `abc_multigauge_generalized_elastic_collision_noise_robustness_bias_control_v2.png` |
