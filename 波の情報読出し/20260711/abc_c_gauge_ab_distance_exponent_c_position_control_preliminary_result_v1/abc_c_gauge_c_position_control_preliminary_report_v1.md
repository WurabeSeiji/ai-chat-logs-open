# ABC C-gauge C-position control preliminary report v1

## Summary

- valid: `True`
- config_count: `160`
- gauge_valid_count: `70`
- gauge_valid_position_pair_count: `23`
- max_pair_abs_alpha_difference: `0.11776992863447344`

## Gauge-valid counts by position

```json
{
  "a_side_large": 6,
  "a_side_large_pi_flip": 6,
  "a_side_small": 11,
  "b_side_large": 6,
  "b_side_large_pi_flip": 6,
  "b_side_small": 11,
  "symmetric": 12,
  "symmetric_pi_flip": 12
}
```

## tau_ABC alpha classes

```json
{
  "proportional_like_alpha_minus1": 70
}
```

## Reading

This is a C-position control test.

No inverse or inverse-square term is injected.

The strict reading is:

```text
C position controls preserve proportional-like tau_ABC readout in gauge-valid cases. Large asymmetric C placement mostly fails by contamination rather than producing a stable inverse law.
```
