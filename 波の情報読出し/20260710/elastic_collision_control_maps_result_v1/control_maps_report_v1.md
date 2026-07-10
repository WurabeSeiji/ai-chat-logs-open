# Control Maps Result v1

## Verdict

- tested_maps: `4`
- reflection_valid_maps: `['reflection']`

## Cases

| map | q_reversed | identity_labels_preserved | identity_labels_swapped | left_slot_mode | right_slot_mode | reflection_spatial_pattern | transmission_spatial_pattern | reflection_valid |
|---|---|---|---|---:|---:|---|---|---|
| reflection | True | True | False | 1 | 2 | True | False | True |
| transmission | False | True | False | 2 | 1 | False | True | False |
| label_exchange_reflection | True | False | True | 2 | 1 | False | True | False |
| transmission_with_label_exchange | False | False | True | 1 | 2 | True | False | False |
