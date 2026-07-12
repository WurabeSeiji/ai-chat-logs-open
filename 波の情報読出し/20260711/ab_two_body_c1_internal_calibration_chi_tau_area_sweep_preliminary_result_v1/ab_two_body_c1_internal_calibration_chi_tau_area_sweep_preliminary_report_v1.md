# AB二体 c=1 内部較正 chi-tau 面積スイープ予備実験レポート v1

## Summary

- valid: `True`
- case_summary_count: `288`
- max_Q_closed_abs: `0.0000000000000000e+00`
- disabled_max_area: `0.0000000000000000e+00`
- locked_max_area: `7.1054273576010019e-15`
- independent_min_area: `2.4367633602631385e-03`
- c1_max_epsilon_c_abs: `3.0629256479024480e-03`
- tau_disabled_rank1_all_cases: `True`
- tau_locked_rank1_all_cases: `True`
- tau_independent_rank2_all_cases: `True`
- tau_is_step_used_any: `False`
- external_c_used_any: `False`
- f_A_or_f_B_used_any: `False`

## Main reading

`tau_disabled_control` and `tau_locked_*` are controls. They keep the readout effectively one-dimensional and do not generate a `chi-tau` area.

`tau_independent_*` modes generate a two-dimensional readout surface. Among them, `tau_independent_c1` is the internally calibrated case: the one-period RMS exchange ratio between `chi_read` and `tau_read` is approximately one.

The power-candidate fit is reported only after a nonzero area sweep is present. In the present preliminary construction, the c1 surface gives alpha values:

```text
-2, -2
```

Negative alpha means that the measured area-like readout grows with the initial phase deviation rather than decays with it. Therefore this preliminary experiment establishes the `chi-tau` surface control, but it does not yet produce an inverse-power decay law.

## Output files

| kind | file |
|---|---|
| JSON | `ab_two_body_c1_internal_calibration_chi_tau_area_sweep_preliminary_result_v1.json` |
| series CSV | `ab_two_body_c1_internal_calibration_chi_tau_area_sweep_series_v1.csv` |
| case summary CSV | `ab_two_body_c1_internal_calibration_chi_tau_area_sweep_case_summary_v1.csv` |
| power candidate CSV | `ab_two_body_c1_internal_calibration_chi_tau_area_sweep_power_candidates_v1.csv` |
| chi-tau surface | `ab_two_body_c1_internal_calibration_chi_tau_surface_v1.png` |
| area sweep | `ab_two_body_c1_internal_calibration_area_sweep_v1.png` |
| c calibration | `ab_two_body_c1_internal_calibration_error_v1.png` |
| power candidate | `ab_two_body_c1_internal_calibration_power_candidate_v1.png` |
| readout leak | `ab_two_body_c1_internal_calibration_readout_leak_v1.png` |
