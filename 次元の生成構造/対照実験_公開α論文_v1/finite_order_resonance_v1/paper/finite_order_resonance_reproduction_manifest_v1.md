# Finite-Order Resonance Reproduction Manifest v1

## Record

- Concept DOI: 10.5281/zenodo.21421366
- Version DOI: 10.5281/zenodo.21421367
- Release date: 2026-07-18

## Paper

- Japanese Markdown, TeX, and PDF
- English Markdown, TeX, and PDF

## Core Source Code

- `run_minimal_system_B_gray_direct_check_v5.py`: uniform $R$ sweep, scattering update, and gray error
- `run_system_A_localization_exchange_R_sweep_preliminary_v1.py`: vector wave-packet localization-transfer calculation
- `phase5_eigenphase_resonance_v2.py`: analytic eigenvalues and finite-order roots
- `run_root_centered_delta1e-12_v1.py`: $R_{124,23}$ root-centered sweep
- `run_R122_23_root_centered_resolution_comparison_v1.py`: $R_{122,23}$ resolution comparison
- `run_high_order_roots_centered_delta1e-12_v1.py`: $R_{567,107}$ and $R_{620,117}$ sweeps
- `run_two_physical_roots_multiprecision_v1.py`: 50- and 80-digit recurrence calculation

## Data

- Full-range sweep statistics at $\Delta R=10^{-7}$
- $R_{124,23}$ root-centered data and summary
- $R_{122,23}$ root-centered data and summary
- $R_{567,107}$ and $R_{620,117}$ root-centered data and summary
- 50- and 80-digit data and summary for $R_{124,23}$ and $R_{620,117}$

## Figures

- Full-range candidate-depth figure
- $R_{124,23}$ root-centered figure
- $R_{122,23}$ resolution-comparison figure
- $R_{620,117}$ root-centered figure
- 80-digit two-root comparison figure

## Minimal Reproduction Sequence

1. Run the System B v5 uniform sweep over the range and spacing listed in the paper.
2. Run the three root-centered scripts without altering their analytic center values.
3. Verify that the deepest grid point coincides with the analytic root.
4. Run the multiprecision script at 50 and 80 digits.
5. Verify that the central depth grows with arithmetic precision and that the best prefix ends at $2n-1$ for the even roots.

