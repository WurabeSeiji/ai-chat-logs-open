# AB二体 chi-tau 面積逆数補償診断予備実験レポート v1

## Summary

- valid: `True`
- diagnostic_case_count: `288`
- area_valid_case_count: `144`
- native_positive2_count: `0`
- constructed_reciprocal_positive2_count: `2`
- c1_readout_off_area_min: `3.8252030346360663e-03`
- c1_readout_off_area_max: `3.4426827311724595e+00`

## Best native alpha candidate

```json
{
  "candidate": "native_max_epsilon_c",
  "fit_valid": true,
  "loglog_slope": -0.011291564598436962,
  "power_candidate_alpha": 0.011291564598436962,
  "loglog_intercept": -33.74994690512192,
  "log_rmse": 0.012683962617493227,
  "case_count": 6,
  "min_value": 2.220446049250313e-15,
  "max_value": 2.3314683517128287e-15,
  "protocol": "Protocol_B",
  "tau_mode": "tau_independent_c1",
  "readout_mode": "readout_off",
  "candidate_kind": "native",
  "alpha_near_positive_2": false
}
```

## Reading

This diagnostic deliberately separates native readouts from constructed reciprocal controls.

`1/A_chi_tau` gives an alpha near `+2` by construction. This is not counted as a native discovery.

In the present dataset, no native readout candidate naturally gives alpha near `+2`.

Therefore the strict reading is:

```text
chi-tau area exists.
1 / chi-tau area has inverse-square scaling by construction.
native inverse-area compensation has not yet been detected.
```
