# ABC Multigauge Interference Readout Asymmetric Amplitude Sweep v1

## Purpose

This diagnostic keeps the simple q-flip collision map and changes only A/B representative amplitudes.
It checks whether individual p/E/R readouts remain reconstructable, and whether R-weighted total momentum remains compatible with the simple equal-mass reflection map.

## Aggregate Verdict

- case_count: `8`
- asymmetric_case_count: `7`
- individual_multigauge_valid_all_cases: `True`
- weighted_energy_preserved_all_cases: `True`
- R_total_preserved_all_cases: `True`
- equal_case_weighted_momentum_preserved: `True`
- asymmetric_cases_detect_weighted_momentum_failure: `True`
- max_p_abs_error: `3.6193270602780103e-14`
- max_E_abs_error: `2.6423307986078726e-14`
- max_R_abs_error: `3.552713678800501e-15`
- max_R_gauge_std: `2.082963028648268e-15`
- max_global_R_contrast_ratio_time: `1195.4814536929925`
- max_within_particle_separation_ratio_time: `3.683878850562613e-30`
- max_weighted_p_collision_error: `16.000000000000036`
- asymmetric_amplitude_diagnostic_valid: `True`

## Case Summary

| case | R_B/R_A | individual valid | within-particle Var(R)/Var(t) | weighted P error | weighted E error | R total error | simple q flip P preserved |
|---|---:|---|---:|---:|---:|---:|---|
| A_1.00_B_1.00 | 1.0000000000000000e+00 | `True` | 0.0000000000000000e+00 | 0.0000000000000000e+00 | 0.0000000000000000e+00 | 0.0000000000000000e+00 | `True` |
| A_1.00_B_1.10 | 1.2100000000000002e+00 | `True` | 0.0000000000000000e+00 | 4.2000000000000137e-01 | 0.0000000000000000e+00 | 0.0000000000000000e+00 | `False` |
| A_1.00_B_1.25 | 1.5625000000000000e+00 | `True` | 3.6838788505626130e-30 | 1.1249999999999964e+00 | 0.0000000000000000e+00 | 0.0000000000000000e+00 | `False` |
| A_1.00_B_1.50 | 2.2500000000000000e+00 | `True` | 0.0000000000000000e+00 | 2.5000000000000124e+00 | 0.0000000000000000e+00 | 0.0000000000000000e+00 | `False` |
| A_1.00_B_2.00 | 4.0000000000000000e+00 | `True` | 0.0000000000000000e+00 | 5.9999999999999858e+00 | 0.0000000000000000e+00 | 0.0000000000000000e+00 | `False` |
| A_1.00_B_3.00 | 9.0000000000000000e+00 | `True` | 0.0000000000000000e+00 | 1.6000000000000036e+01 | 0.0000000000000000e+00 | 0.0000000000000000e+00 | `False` |
| A_1.50_B_1.00 | 4.4444444444444442e-01 | `True` | 0.0000000000000000e+00 | 2.5000000000000124e+00 | 0.0000000000000000e+00 | 0.0000000000000000e+00 | `False` |
| A_2.00_B_1.00 | 2.5000000000000000e-01 | `True` | 0.0000000000000000e+00 | 5.9999999999999858e+00 | 0.0000000000000000e+00 | 0.0000000000000000e+00 | `False` |

## Interpretation

The individual multigauge p/E/R readout remains valid in every amplitude case.
However, the R-weighted total momentum is conserved by the simple q-flip map only in the equal-amplitude case.
The global R variance is not used as a failure condition here, because unequal A/B amplitudes are the diagnostic signal itself.
Thus the experiment separates readout validity from collision-law compatibility.

## Files

| kind | file |
|---|---|
| JSON | `abc_multigauge_interference_readout_asymmetric_amplitude_sweep_result_v2.json` |
| case CSV | `abc_multigauge_interference_readout_asymmetric_amplitude_cases_v2.csv` |
| stage quantity CSV | `abc_multigauge_interference_readout_asymmetric_amplitude_stage_quantities_v2.csv` |
| gauge CSV | `abc_multigauge_interference_readout_asymmetric_amplitude_gauge_rows_v2.csv` |
| conservation plot | `abc_multigauge_interference_readout_asymmetric_amplitude_conservation_v2.png` |
| cases plot | `abc_multigauge_interference_readout_asymmetric_amplitude_cases_v2.png` |
