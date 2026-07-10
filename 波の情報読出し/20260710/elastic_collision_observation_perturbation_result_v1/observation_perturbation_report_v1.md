# Observation Perturbation Result v1

## Verdict

- total_cases: `8`
- model_valid_cases: `5`
- collision_map_valid_cases: `7`
- invalid_cases: `3`

## Cases

| case | max_chi/eps_C | max_tau/eps_C | tau_diff/eps_AB | observation_bound_respected | collision_cell_reached | collision_map_valid | model_valid |
|---|---:|---:|---:|---|---|---|---|
| none | 0 | 0 | 0 | True | True | True | True |
| common_spatial_within_C | 0.8 | 0 | 0 | True | True | True | True |
| differential_spatial_within_C | 1 | 0 | 0 | True | True | True | True |
| differential_time_within_C | 0 | 1 | 0.2 | True | True | True | True |
| mixed_corner_within_C | 1 | 1 | 0.2 | True | True | True | True |
| over_C_time_still_in_AB | 0 | 3 | 0.6 | False | True | True | False |
| over_C_time_breaks_AB | 0 | 8 | 1.6 | False | False | False | False |
| over_C_spatial_only | 10 | 0 | 0 | False | True | True | False |
