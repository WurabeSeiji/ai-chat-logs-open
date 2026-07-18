# AB二体 chi-tau native 逆面積 extended sweep 予備実験レポート v1

## Summary

- valid: `True`
- sweep_case_count: `1323`
- area_valid_case_count: `1260`
- c1_surface_like_case_count: `126`
- native_positive2_count: `0`
- c1_native_positive2_count: `0`
- constructed_reciprocal_positive2_count: `198`
- c1_area_min: `2.4281646266173604e-04`
- c1_area_max: `3.4426827311724595e+00`

## Best native alpha candidate

```json
{
  "candidate": "native_max_epsilon_c",
  "fit_valid": true,
  "loglog_slope": -0.013217738553867962,
  "power_candidate_alpha": 0.013217738553867962,
  "loglog_intercept": -33.75183974361942,
  "log_rmse": 0.012175393351613396,
  "case_count": 7,
  "min_value": 2.220446049250313e-15,
  "max_value": 2.3314683517128287e-15,
  "frequency_ratio": 1.0,
  "amplitude_ratio": 1.0,
  "phase_shift_deg": 0.0,
  "readout_mode": "readout_off",
  "filter": "area_valid",
  "candidate_kind": "native",
  "alpha_near_positive_2": false
}
```

## Reading

This is an extended negative-control search.

The sweep varies initial deviation, frequency ratio, amplitude ratio, phase shift, and readout leak.

The strict result is:

```text
native inverse-area scaling was not detected.
constructed reciprocal controls give alpha near +2.
```

This keeps the inverse-square claim on hold.
