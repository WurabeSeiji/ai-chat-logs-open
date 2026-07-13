# Cell Resolution Sweep Result v1

## Verdict

- total_cases: `60`
- valid_cases: `57`
- invalid_cases: `3`
- offgrid_valid_cases: `27`
- offgrid_invalid_cases: `3`

## Interpretation

The finite collision cell is reliably detected when the calculation step is not larger than the cell width.
Aligned cases can pass accidentally by landing exactly on the center; off-grid cases expose skipped-cell failures.

## Cases

| d0 | Nh_chi_AB | delta_s | epsilon_chi_AB | min_abs_chi_gap | collision_cell_reached | crossed_without_detection | case_valid |
|---:|---:|---:|---:|---:|---|---|---|
| 0.2 | 19 | 0.02 | 0.157079632679 | 0.12 | True | False | True |
| 0.2 | 19 | 0.01 | 0.157079632679 | 0.14 | True | False | True |
| 0.2 | 19 | 0.005 | 0.157079632679 | 0.15 | True | False | True |
| 0.2 | 19 | 0.002 | 0.157079632679 | 0.156 | True | False | True |
| 0.2 | 19 | 0.001 | 0.157079632679 | 0.156 | True | False | True |
| 0.2 | 39 | 0.02 | 0.0785398163397 | 0.04 | True | False | True |
| 0.2 | 39 | 0.01 | 0.0785398163397 | 0.06 | True | False | True |
| 0.2 | 39 | 0.005 | 0.0785398163397 | 0.07 | True | False | True |
| 0.2 | 39 | 0.002 | 0.0785398163397 | 0.076 | True | False | True |
| 0.2 | 39 | 0.001 | 0.0785398163397 | 0.078 | True | False | True |
| 0.2 | 99 | 0.02 | 0.0314159265359 | 4.16333634234e-17 | True | False | True |
| 0.2 | 99 | 0.01 | 0.0314159265359 | 0.02 | True | False | True |
| 0.2 | 99 | 0.005 | 0.0314159265359 | 0.03 | True | False | True |
| 0.2 | 99 | 0.002 | 0.0314159265359 | 0.028 | True | False | True |
| 0.2 | 99 | 0.001 | 0.0314159265359 | 0.03 | True | False | True |
| 0.2 | 199 | 0.02 | 0.0157079632679 | 4.16333634234e-17 | True | False | True |
| 0.2 | 199 | 0.01 | 0.0157079632679 | 6.24500451352e-17 | True | False | True |
| 0.2 | 199 | 0.005 | 0.0157079632679 | 0.01 | True | False | True |
| 0.2 | 199 | 0.002 | 0.0157079632679 | 0.012 | True | False | True |
| 0.2 | 199 | 0.001 | 0.0157079632679 | 0.014 | True | False | True |
| 0.2 | 399 | 0.02 | 0.00785398163397 | 4.16333634234e-17 | True | False | True |
| 0.2 | 399 | 0.01 | 0.00785398163397 | 6.24500451352e-17 | True | False | True |
| 0.2 | 399 | 0.005 | 0.00785398163397 | 1.97758476261e-16 | True | False | True |
| 0.2 | 399 | 0.002 | 0.00785398163397 | 0.004 | True | False | True |
| 0.2 | 399 | 0.001 | 0.00785398163397 | 0.006 | True | False | True |
| 0.2 | 999 | 0.02 | 0.00314159265359 | 4.16333634234e-17 | True | False | True |
| 0.2 | 999 | 0.01 | 0.00314159265359 | 6.24500451352e-17 | True | False | True |
| 0.2 | 999 | 0.005 | 0.00314159265359 | 1.97758476261e-16 | True | False | True |
| 0.2 | 999 | 0.002 | 0.00314159265359 | 3.05311331772e-16 | True | False | True |
| 0.2 | 999 | 0.001 | 0.00314159265359 | 0.002 | True | False | True |
| 0.203 | 19 | 0.02 | 0.157079632679 | 0.126 | True | False | True |
| 0.203 | 19 | 0.01 | 0.157079632679 | 0.146 | True | False | True |
| 0.203 | 19 | 0.005 | 0.157079632679 | 0.156 | True | False | True |
| 0.203 | 19 | 0.002 | 0.157079632679 | 0.154 | True | False | True |
| 0.203 | 19 | 0.001 | 0.157079632679 | 0.156 | True | False | True |
| 0.203 | 39 | 0.02 | 0.0785398163397 | 0.046 | True | False | True |
| 0.203 | 39 | 0.01 | 0.0785398163397 | 0.066 | True | False | True |
| 0.203 | 39 | 0.005 | 0.0785398163397 | 0.076 | True | False | True |
| 0.203 | 39 | 0.002 | 0.0785398163397 | 0.078 | True | False | True |
| 0.203 | 39 | 0.001 | 0.0785398163397 | 0.078 | True | False | True |
| 0.203 | 99 | 0.02 | 0.0314159265359 | 0.006 | True | False | True |
| 0.203 | 99 | 0.01 | 0.0314159265359 | 0.026 | True | False | True |
| 0.203 | 99 | 0.005 | 0.0314159265359 | 0.026 | True | False | True |
| 0.203 | 99 | 0.002 | 0.0314159265359 | 0.03 | True | False | True |
| 0.203 | 99 | 0.001 | 0.0314159265359 | 0.03 | True | False | True |
| 0.203 | 199 | 0.02 | 0.0157079632679 | 0.006 | True | False | True |
| 0.203 | 199 | 0.01 | 0.0157079632679 | 0.006 | True | False | True |
| 0.203 | 199 | 0.005 | 0.0157079632679 | 0.006 | True | False | True |
| 0.203 | 199 | 0.002 | 0.0157079632679 | 0.014 | True | False | True |
| 0.203 | 199 | 0.001 | 0.0157079632679 | 0.014 | True | False | True |
| 0.203 | 399 | 0.02 | 0.00785398163397 | 0.006 | True | False | True |
| 0.203 | 399 | 0.01 | 0.00785398163397 | 0.006 | True | False | True |
| 0.203 | 399 | 0.005 | 0.00785398163397 | 0.006 | True | False | True |
| 0.203 | 399 | 0.002 | 0.00785398163397 | 0.006 | True | False | True |
| 0.203 | 399 | 0.001 | 0.00785398163397 | 0.006 | True | False | True |
| 0.203 | 999 | 0.02 | 0.00314159265359 | 0.006 | False | True | False |
| 0.203 | 999 | 0.01 | 0.00314159265359 | 0.006 | False | True | False |
| 0.203 | 999 | 0.005 | 0.00314159265359 | 0.004 | False | True | False |
| 0.203 | 999 | 0.002 | 0.00314159265359 | 0.002 | True | False | True |
| 0.203 | 999 | 0.001 | 0.00314159265359 | 0.002 | True | False | True |
