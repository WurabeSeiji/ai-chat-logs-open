# Independent rank-3 geometry audit (2026-08-31)

This folder is an independent reproducibility package for the hm-series geometry question.
It deliberately does **not** use `hm_series_k.csv` or pass13-pass17 result files as numerical input.
The only numerical inputs are the raw `data/hm_N*/states_treatment.npz` files.

## Layout

- `program/reproduce_hm_rank3_audit.py`: independent readout + Gram-rank audit.
- `program/run_all.sh`: canonical N=3..16 runner.
- `results/validation_hm_N6_edge_wavelengths.csv`: N=6 readout from the raw NPZ for windows 4096/8192/16384.
- `results/validation_hm_N6_rank3_group_baseline.csv`: N=6 baseline geometry result.
- `docs/analysis_hm_rank3_audit_20260831.md`: interpretation, limitations, and next test.
- `SHA256SUMS.txt`: hashes of this package plus the N=6 raw source used for validation.

## Canonical execution

Run from the project checkout where this audit folder and `data/` are siblings:

```bash
bash independent_rank3_geometry_audit_20260831/program/run_all.sh
```

Requirements: Python 3, NumPy, SciPy.

## Source-data contract

Expected source path for each N:

```text
data/hm_N{N}/states_treatment.npz
```

Each NPZ must contain array `Z` with shape `(T, N(N-1)/2)`. The N=6 validation source had shape `(40001, 15)`.

N=6 raw-source SHA-256 used for the validation files:

```text
0d523192ddd2a2f2f1ec6a67f19016099a212de75290774bc876b4d5e4c2fe02  data/hm_N6/states_treatment.npz
```

## Important scope boundary

The included integer-k audit reproduces the old **2% group-uniform ansatz only as a baseline**. A PSD solution of arbitrary rank is not called success. The physical target for N>4 is rank(B)=3. A full edge-wise independent integer search is not claimed here and must be implemented separately.
