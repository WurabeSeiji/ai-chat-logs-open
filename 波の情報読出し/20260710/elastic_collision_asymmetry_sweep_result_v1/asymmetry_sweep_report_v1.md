# Asymmetry Sweep Result v1

## Verdict

- total_cases: `10`
- valid_cases: `7`
- invalid_cases: `3`

## Cases

| case | A_A | A_B | Nh_chi_A | Nh_chi_B | tau_gap_initial | omega_A | omega_B | min_tau_gap | collision_cell_reached | spatial_crossed_without_time_cell | case_valid |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|
| symmetric_baseline | 1 | 1 | 99 | 99 | 0 | 1 | 1 | 0 | True | False | True |
| amplitude_asymmetry | 1 | 2 | 99 | 99 | 0 | 1 | 1 | 0 | True | False | True |
| spatial_harmonic_asymmetry | 1 | 1 | 199 | 99 | 0 | 1 | 1 | 0 | True | False | True |
| temporal_harmonic_asymmetry | 1 | 1 | 99 | 99 | 0 | 1 | 1 | 0 | True | False | True |
| small_time_phase_offset | 1 | 1 | 99 | 99 | 0.02 | 1 | 1 | 0.02 | True | False | True |
| large_time_phase_offset | 1 | 1 | 99 | 99 | 0.05 | 1 | 1 | 0.05 | False | True | False |
| small_omega_mismatch | 1 | 1 | 99 | 99 | 0 | 1 | 1.05 | 0 | True | False | True |
| large_omega_mismatch | 1 | 1 | 99 | 99 | 0 | 1 | 1.5 | 0 | False | True | False |
| combined_asymmetry_pass | 1.5 | 0.8 | 199 | 99 | 0.02 | 1 | 1.1 | 0.0015 | True | False | True |
| combined_asymmetry_fail | 1.5 | 0.8 | 199 | 99 | 0.05 | 1 | 1.5 | 0.05 | False | True | False |
