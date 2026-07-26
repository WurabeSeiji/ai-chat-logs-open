# README — Paper 7: Emergence of a Three-Direction Space in an N-Body Relational-Wave Closed System

- Concept DOI: 10.5281/zenodo.21578401
- Version DOI (v1.0): 10.5281/zenodo.21578402
- Author: Noriaki Kihara (WF System Co., Ltd.), ORCID 0009-0004-6753-4020
- License: CC-BY-4.0

This record contains the full paper text, figures, the complete reproduction programs, and the numerical results.

## Contents

### Paper (JP / EN)
- `paper7_ja.md`, `paper7_ja.tex`, `paper7_ja.pdf` — Japanese
- `paper7_en.md`, `paper7_en.tex`, `paper7_en.pdf` — English

### Figures
- `figure1_compare_N5_N40_N300.png` — splitting fraction f (N=5,40,300)
- `figure2_compare_N5_N40_N300.png` — five-component occupation (stack)
- `figure3_compare_N5_N40_N300.png` — five-component occupation (log)
- `transverse_growth_compare_N5_N40_N300.png` — transverse-perturbation response

### Reproduction programs (Python 3)
Dependency order (upper depends on lower):
1. `run_n_scaling_lowrank_v1.py` — original engine (K=WJW^T, Cayley update). SHA-256 `ba0fc19b03caf06d16e97e2cd5da499ed2ed8e95288f6dbf2062bedcaf11176d`. Unmodified.
2. `run_plane_flow_exact_v1.py` — fixed parent basis, dense eig (N=5,40)
3. `run_plane_flow_approx_v1.py` — fixed parent basis, low-rank JG (N=300)
4. `run_n300_dimension_saturation_v2.py` — dominant plane B_dom via Gram reduction G=W^T W
5. `run_paper7_5color_timeseries.py` — five-component occupation time series (Figures 1/2/3)
6. `run_paper7_transverse.py` — transverse stability (perturbation in S4(t)^perp, Benettin)
7. `run_paper7_transverse_cached.py` — transverse stability, baseline-cached (N=300)
8. `run_paper7_exact_vs_approx_N40.py` — Section 6.2 exact-vs-low-rank check
9. `make_paper7_figures.py` — figure generation from CSVs

### Numerical results
- `paper7_long_timeseries_N00005.csv`, `_N00040.csv`, `_N00300.csv` — five-component occupation
- `transverse_stability_timeseries_N00005.csv`, `_N00040.csv`, `_N00300.csv` — transverse response
- `N_comparison_table.csv`, `transverse_stability_summary.csv` — summaries
- `paper7_longtime_and_transverse_stability_report.md` — numerical report (values-only)

## Execution conditions
- Python 3.x, NumPy, SciPy (1.13), Matplotlib.
- Deterministic: seed = 40260722 + 1000*N; initial seed delta = 1e-15; Cayley gamma = tan(pi/144).
- Common absolute time axis 0..55000 (crossing not shifted to 0).
- Run (from a directory where the imported modules resolve; the programs use relative sys.path inserts matching the source repository layout):
  - `python3 run_paper7_5color_timeseries.py 5` (also 40, 300)
  - `python3 run_paper7_transverse.py 5` / `run_paper7_transverse_cached.py 300`
  - `python3 make_paper7_figures.py`
- For N=40, the dense-matrix and low-rank methods agree within double precision (max deviation 1.78e-15).

## Notes
- Numerical diagnostics: conservation error ≤ 2.0e-15; five-component projection-closure error ≤ 2.2e-16.
- Section 6.4 documents an open implementation issue in the transverse Benettin warm-start synchronization (affects the quantitative lambda_perp only, not rank Q:2→4, occupations, or the closure).
