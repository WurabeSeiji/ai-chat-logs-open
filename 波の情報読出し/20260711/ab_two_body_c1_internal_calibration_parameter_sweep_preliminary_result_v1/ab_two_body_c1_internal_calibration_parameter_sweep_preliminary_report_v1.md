# AB二体 c=1 内部較正パラメータスイープ予備実験レポート v1

## Summary

- valid: `True`
- sweep_case_count: `756`
- readout_off_case_count: `252`
- rank2_readout_off_count: `246`
- c1_surface_like_readout_off_count: `13`
- c1_locked_like_readout_off_count: `8`
- min_c_error_readout_off: `1.1102230246251565e-15`
- max_area_readout_off: `1.9685536479742288e-01`

## Best readout_off config

```json
{
  "frequency_ratio": 1.0,
  "amplitude_ratio": 1.0,
  "phase_shift_deg": 75.0,
  "amp_freq_product": 1.0,
  "readout_mode": "readout_off",
  "per_step_leak": 0.0,
  "max_abs_A_chi_tau": 0.024750884918694482,
  "rank_chi_tau": 2,
  "max_epsilon_c_abs": 1.1102230246251565e-15,
  "locked_like": false,
  "c1_like": true,
  "c1_surface_like": true,
  "decay_rate_envelope": 0.0,
  "chi_min": -0.17453292519943295,
  "chi_max": 0.17453292519943295,
  "tau_min": -0.17453292519943295,
  "tau_max": 0.17453292519943295,
  "tau_is_step_used": false,
  "external_c_used": false,
  "absolute_background_axis_used": false,
  "f_A_or_f_B_used": false
}
```

## Reading

This sweep shows that internal `c=1` calibration alone is not sufficient.

Some locked-like configurations can satisfy the RMS exchange ratio while producing no independent `chi-tau` area.

Therefore, the necessary readout conditions are:

```text
c=1 calibration
rank_chi_tau = 2
A_chi_tau != 0
```

The result strengthens the control discipline for the next experiment.
