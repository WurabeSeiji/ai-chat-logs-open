# Eta Resolution Sweep Result v1

## Verdict

- total_cases: `88`
- valid_cases: `59`
- invalid_cases: `29`
- alias_collision_cases: `29`
- non_alias_failures: `0`
- first_eta_samples_all_pairs_valid: `64`

## Interpretation

The internal identification vibration is readable only up to the eta-readout sampling resolution.
When the mode difference is a multiple of the eta sample count, the two label modes alias and the readout becomes ambiguous.
This is a readout-resolution failure, not a failure of the collision map itself.

## Cases

| eta samples | m_A | m_B | diff | alias | min purity | detected A0 | detected B0 | valid |
|---:|---:|---:|---:|---|---:|---|---|---|
| 4 | 1 | 2 | 1 | False | 1 | 1 | 2 | True |
| 4 | 1 | 5 | 4 | True | 0.5 | ambiguous | ambiguous | False |
| 4 | 1 | 9 | 8 | True | 0.5 | ambiguous | ambiguous | False |
| 4 | 1 | 17 | 16 | True | 0.5 | ambiguous | ambiguous | False |
| 4 | 1 | 33 | 32 | True | 0.5 | ambiguous | ambiguous | False |
| 4 | 8 | 9 | 1 | False | 1 | 8 | 9 | True |
| 4 | 15 | 31 | 16 | True | 0.5 | ambiguous | ambiguous | False |
| 4 | 16 | 32 | 16 | True | 0.5 | ambiguous | ambiguous | False |
| 4 | 17 | 33 | 16 | True | 0.5 | ambiguous | ambiguous | False |
| 4 | 5 | 29 | 24 | True | 0.5 | ambiguous | ambiguous | False |
| 4 | 7 | 31 | 24 | True | 0.5 | ambiguous | ambiguous | False |
| 6 | 1 | 2 | 1 | False | 1 | 1 | 2 | True |
| 6 | 1 | 5 | 4 | False | 1 | 1 | 5 | True |
| 6 | 1 | 9 | 8 | False | 1 | 1 | 9 | True |
| 6 | 1 | 17 | 16 | False | 1 | 1 | 17 | True |
| 6 | 1 | 33 | 32 | False | 1 | 1 | 33 | True |
| 6 | 8 | 9 | 1 | False | 1 | 8 | 9 | True |
| 6 | 15 | 31 | 16 | False | 1 | 15 | 31 | True |
| 6 | 16 | 32 | 16 | False | 1 | 16 | 32 | True |
| 6 | 17 | 33 | 16 | False | 1 | 17 | 33 | True |
| 6 | 5 | 29 | 24 | True | 0.5 | ambiguous | ambiguous | False |
| 6 | 7 | 31 | 24 | True | 0.5 | ambiguous | ambiguous | False |
| 8 | 1 | 2 | 1 | False | 1 | 1 | 2 | True |
| 8 | 1 | 5 | 4 | False | 1 | 1 | 5 | True |
| 8 | 1 | 9 | 8 | True | 0.5 | ambiguous | ambiguous | False |
| 8 | 1 | 17 | 16 | True | 0.5 | ambiguous | ambiguous | False |
| 8 | 1 | 33 | 32 | True | 0.5 | ambiguous | ambiguous | False |
| 8 | 8 | 9 | 1 | False | 1 | 8 | 9 | True |
| 8 | 15 | 31 | 16 | True | 0.5 | ambiguous | ambiguous | False |
| 8 | 16 | 32 | 16 | True | 0.5 | ambiguous | ambiguous | False |
| 8 | 17 | 33 | 16 | True | 0.5 | ambiguous | ambiguous | False |
| 8 | 5 | 29 | 24 | True | 0.5 | ambiguous | ambiguous | False |
| 8 | 7 | 31 | 24 | True | 0.5 | ambiguous | ambiguous | False |
| 12 | 1 | 2 | 1 | False | 1 | 1 | 2 | True |
| 12 | 1 | 5 | 4 | False | 1 | 1 | 5 | True |
| 12 | 1 | 9 | 8 | False | 1 | 1 | 9 | True |
| 12 | 1 | 17 | 16 | False | 1 | 1 | 17 | True |
| 12 | 1 | 33 | 32 | False | 1 | 1 | 33 | True |
| 12 | 8 | 9 | 1 | False | 1 | 8 | 9 | True |
| 12 | 15 | 31 | 16 | False | 1 | 15 | 31 | True |
| 12 | 16 | 32 | 16 | False | 1 | 16 | 32 | True |
| 12 | 17 | 33 | 16 | False | 1 | 17 | 33 | True |
| 12 | 5 | 29 | 24 | True | 0.5 | ambiguous | ambiguous | False |
| 12 | 7 | 31 | 24 | True | 0.5 | ambiguous | ambiguous | False |
| 16 | 1 | 2 | 1 | False | 1 | 1 | 2 | True |
| 16 | 1 | 5 | 4 | False | 1 | 1 | 5 | True |
| 16 | 1 | 9 | 8 | False | 1 | 1 | 9 | True |
| 16 | 1 | 17 | 16 | True | 0.5 | ambiguous | ambiguous | False |
| 16 | 1 | 33 | 32 | True | 0.5 | ambiguous | ambiguous | False |
| 16 | 8 | 9 | 1 | False | 1 | 8 | 9 | True |
| 16 | 15 | 31 | 16 | True | 0.5 | ambiguous | ambiguous | False |
| 16 | 16 | 32 | 16 | True | 0.5 | ambiguous | ambiguous | False |
| 16 | 17 | 33 | 16 | True | 0.5 | ambiguous | ambiguous | False |
| 16 | 5 | 29 | 24 | False | 1 | 5 | 29 | True |
| 16 | 7 | 31 | 24 | False | 1 | 7 | 31 | True |
| 24 | 1 | 2 | 1 | False | 1 | 1 | 2 | True |
| 24 | 1 | 5 | 4 | False | 1 | 1 | 5 | True |
| 24 | 1 | 9 | 8 | False | 1 | 1 | 9 | True |
| 24 | 1 | 17 | 16 | False | 1 | 1 | 17 | True |
| 24 | 1 | 33 | 32 | False | 1 | 1 | 33 | True |
| 24 | 8 | 9 | 1 | False | 1 | 8 | 9 | True |
| 24 | 15 | 31 | 16 | False | 1 | 15 | 31 | True |
| 24 | 16 | 32 | 16 | False | 1 | 16 | 32 | True |
| 24 | 17 | 33 | 16 | False | 1 | 17 | 33 | True |
| 24 | 5 | 29 | 24 | True | 0.5 | ambiguous | ambiguous | False |
| 24 | 7 | 31 | 24 | True | 0.5 | ambiguous | ambiguous | False |
| 32 | 1 | 2 | 1 | False | 1 | 1 | 2 | True |
| 32 | 1 | 5 | 4 | False | 1 | 1 | 5 | True |
| 32 | 1 | 9 | 8 | False | 1 | 1 | 9 | True |
| 32 | 1 | 17 | 16 | False | 1 | 1 | 17 | True |
| 32 | 1 | 33 | 32 | True | 0.5 | ambiguous | ambiguous | False |
| 32 | 8 | 9 | 1 | False | 1 | 8 | 9 | True |
| 32 | 15 | 31 | 16 | False | 1 | 15 | 31 | True |
| 32 | 16 | 32 | 16 | False | 1 | 16 | 32 | True |
| 32 | 17 | 33 | 16 | False | 1 | 17 | 33 | True |
| 32 | 5 | 29 | 24 | False | 1 | 5 | 29 | True |
| 32 | 7 | 31 | 24 | False | 1 | 7 | 31 | True |
| 64 | 1 | 2 | 1 | False | 1 | 1 | 2 | True |
| 64 | 1 | 5 | 4 | False | 1 | 1 | 5 | True |
| 64 | 1 | 9 | 8 | False | 1 | 1 | 9 | True |
| 64 | 1 | 17 | 16 | False | 1 | 1 | 17 | True |
| 64 | 1 | 33 | 32 | False | 1 | 1 | 33 | True |
| 64 | 8 | 9 | 1 | False | 1 | 8 | 9 | True |
| 64 | 15 | 31 | 16 | False | 1 | 15 | 31 | True |
| 64 | 16 | 32 | 16 | False | 1 | 16 | 32 | True |
| 64 | 17 | 33 | 16 | False | 1 | 17 | 33 | True |
| 64 | 5 | 29 | 24 | False | 1 | 5 | 29 | True |
| 64 | 7 | 31 | 24 | False | 1 | 7 | 31 | True |
