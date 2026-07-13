# ABC Multigauge Interference Readout Robustness Sweep v1

## Purpose

This sweep changes only the readout-gauge family while keeping the same one-collision ABC trajectory.
It checks whether p-like, E-like, and R-like readouts are stable under reference phase, readout center, width, and gain changes.

## Aggregate Verdict

- case_count: `5`
- total_gauge_count: `130`
- all_cases_valid: `True`
- max_p_abs_error_all_cases: `3.0331293032759277e-13`
- max_E_abs_error_all_cases: `3.0331293032759277e-13`
- max_R_abs_error_all_cases: `1.5765166949677223e-14`
- max_R_gauge_std_all_cases: `5.288392122597181e-15`
- max_separation_ratio_time_all_cases: `1.3338999651354898e-27`
- single_gauge_only_used: `False`
- robustness_sweep_valid: `True`

## Case Summary

| case | gauges | p max err | E max err | R max err | R std | Var(R)/Var(t) | valid |
|---|---:|---:|---:|---:|---:|---:|---|
| baseline_default | 8 | 2.5202062658991053e-14 | 2.2315482794965646e-14 | 4.4408920985006262e-16 | 2.2204460492503131e-16 | 3.6838616474030686e-30 | `True` |
| phase_center_grid | 45 | 2.3137047833188262e-13 | 2.1471713296250527e-13 | 5.7731597280508140e-15 | 1.7775899925091300e-15 | 2.8276279738980642e-28 | `True` |
| width_gain_grid | 36 | 1.0469403122215226e-13 | 6.5503158452884236e-14 | 4.4408920985006262e-16 | 2.5639502485114184e-16 | 4.9118384674168162e-30 | `True` |
| near_lobe_offset | 25 | 3.0331293032759277e-13 | 3.0331293032759277e-13 | 1.5765166949677223e-14 | 5.2883921225971813e-15 | 1.3338999651354898e-27 | `True` |
| mixed_readout_grid | 16 | 1.6031620475587260e-13 | 1.3566925360919413e-13 | 1.1102230246251565e-14 | 2.5727487310015434e-15 | 1.0176441594089155e-27 | `True` |

## Files

| kind | file |
|---|---|
| JSON | `abc_multigauge_interference_readout_robustness_sweep_result_v2.json` |
| case CSV | `abc_multigauge_interference_readout_robustness_sweep_cases_v2.csv` |
| gauge CSV | `abc_multigauge_interference_readout_robustness_sweep_gauge_rows_v2.csv` |
| stage summary CSV | `abc_multigauge_interference_readout_robustness_sweep_stage_summary_v2.csv` |
| error plot | `abc_multigauge_interference_readout_robustness_sweep_errors_v2.png` |
| stability plot | `abc_multigauge_interference_readout_robustness_sweep_stability_v2.png` |
