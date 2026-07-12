# ABC C-gauge AB distance exponent preliminary report v1

## Summary

- valid: `True`
- config_count: `160`
- gauge_valid_count: `41`
- resolution_valid_count: `60`
- clock_valid_count: `100`
- disturbance_valid_count: `101`
- inverse_like_alpha1_count: `0`
- inverse_square_like_alpha2_count: `0`
- proportional_like_alpha_minus1_count: `41`

## Alpha classes in gauge-valid cases

```json
{
  "proportional_like_alpha_minus1": 41
}
```

## Alpha classes by relation time

```json
{
  "tau_ABC": {
    "proportional_like_alpha_minus1": 41
  },
  "tau_AB": {
    "other": 41
  },
  "tau_AC": {
    "proportional_like_alpha_minus1": 41
  },
  "tau_BC": {
    "proportional_like_alpha_minus1": 41
  }
}
```

## Max alpha delta from tau_ABC

```json
{
  "tau_AB": 0.6588265529372437,
  "tau_AC": 2.220446049250313e-16,
  "tau_BC": 4.440892098500626e-16
}
```

## Reading

This preliminary test does not inject `1/L` or `1/L^2`.

It checks whether a third wave `C` can serve as an independent space-time gauge for reading the AB relation compensation.

The main time readout is `tau_ABC`.

Relation-time diagnostics for `tau_AB`, `tau_AC`, and `tau_BC` are also recorded, but they are not used as the main success criterion.

The strict reading is:

```text
The C gauge has a finite eligibility window. Within this minimal AB-harmonic preliminary model, gauge-valid cases recover proportional-like alpha=-1, not inverse or inverse-square. This is a gauge-eligibility test, not a final physical-law derivation.
```
