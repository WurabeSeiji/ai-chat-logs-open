# ABC C-gauge relation decomposition preliminary report v1

## Summary

- valid: `True`
- config_count: `240`
- decomposition_valid_count: `180`
- AB_dominant_valid_count: `48`
- non_AB_dominant_count: `132`
- inverse_or_inverse_square_in_AB_dominant_count: `0`

## Projection mode counts

```json
{
  "a_assists_b_opposes": 45,
  "a_opposes_b_assists": 45,
  "c_assists_ab": 45,
  "c_opposes_ab": 45
}
```

## Pair tau_ABC alpha classes in AB-dominant cases

```json
{
  "proportional_like_alpha_minus1": 48
}
```

## Pair tau_ABC alpha classes in all decomposition-valid cases

```json
{
  "constant_like_alpha0": 3,
  "other": 46,
  "proportional_like_alpha_minus1": 131
}
```

## Native AB alpha classes

```json
{
  "proportional_like_alpha_minus1": 180
}
```

## C-bias alpha classes

```json
{
  "constant_like_alpha0": 150,
  "unfit": 30
}
```

## C-asymmetry alpha classes

```json
{
  "constant_like_alpha0": 150,
  "unfit": 30
}
```

## Common-mode direct-control alpha classes

```json
{
  "other": 63,
  "proportional_like_alpha_minus1": 117
}
```

## Relation-time alpha deltas in AB-dominant cases

```text
max |tau_AB - tau_ABC| alpha delta: 0.6588265529372439
max |tau_AC - tau_ABC| alpha delta: 6.661338147750939e-16
max |tau_BC - tau_ABC| alpha delta: 6.661338147750939e-16
```

## Reading

No inverse or inverse-square term is injected.

The strict reading is:

```text
Separating f_AB, f_AC, f_BC, and f_ABC confirms that tau_ABC can be used as a representative time gauge, while f_ABC should not be added as a direct circular term. In AB-dominant projected cases, the pair readout remains proportional-like rather than inverse-like.
```
