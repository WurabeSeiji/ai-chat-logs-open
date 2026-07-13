# Identification Vibration Robustness Result v1

## Verdict

- total_cases: `36`
- valid_cases: `20`
- invalid_cases: `16`
- first_failure_leakage: `0.35`

## Interpretation

The identification channel is preserved while the target mode remains the dominant internal vibration.
The transition point is the expected mode-mixture boundary near leakage = 0.5.

## Cases

| m_A | m_B | leakage | min_purity | label_preserved | label_swapped | case_valid |
|---:|---:|---:|---:|---|---|---|
| 1 | 2 | 0 | 1 | True | False | True |
| 1 | 2 | 0.01 | 0.99 | True | False | True |
| 1 | 2 | 0.05 | 0.95 | True | False | True |
| 1 | 2 | 0.1 | 0.9 | True | False | True |
| 1 | 2 | 0.2 | 0.8 | True | False | True |
| 1 | 2 | 0.35 | 0.65 | True | False | False |
| 1 | 2 | 0.49 | 0.51 | True | False | False |
| 1 | 2 | 0.5 | 0.5 | False | True | False |
| 1 | 2 | 0.51 | 0.49 | False | True | False |
| 1 | 3 | 0 | 1 | True | False | True |
| 1 | 3 | 0.01 | 0.99 | True | False | True |
| 1 | 3 | 0.05 | 0.95 | True | False | True |
| 1 | 3 | 0.1 | 0.9 | True | False | True |
| 1 | 3 | 0.2 | 0.8 | True | False | True |
| 1 | 3 | 0.35 | 0.65 | True | False | False |
| 1 | 3 | 0.49 | 0.51 | True | False | False |
| 1 | 3 | 0.5 | 0.5 | False | True | False |
| 1 | 3 | 0.51 | 0.49 | False | True | False |
| 2 | 5 | 0 | 1 | True | False | True |
| 2 | 5 | 0.01 | 0.99 | True | False | True |
| 2 | 5 | 0.05 | 0.95 | True | False | True |
| 2 | 5 | 0.1 | 0.9 | True | False | True |
| 2 | 5 | 0.2 | 0.8 | True | False | True |
| 2 | 5 | 0.35 | 0.65 | True | False | False |
| 2 | 5 | 0.49 | 0.51 | True | False | False |
| 2 | 5 | 0.5 | 0.5 | False | True | False |
| 2 | 5 | 0.51 | 0.49 | False | True | False |
| 8 | 9 | 0 | 1 | True | False | True |
| 8 | 9 | 0.01 | 0.99 | True | False | True |
| 8 | 9 | 0.05 | 0.95 | True | False | True |
| 8 | 9 | 0.1 | 0.9 | True | False | True |
| 8 | 9 | 0.2 | 0.8 | True | False | True |
| 8 | 9 | 0.35 | 0.65 | True | False | False |
| 8 | 9 | 0.49 | 0.51 | True | False | False |
| 8 | 9 | 0.5 | 0.5 | False | True | False |
| 8 | 9 | 0.51 | 0.49 | False | True | False |
