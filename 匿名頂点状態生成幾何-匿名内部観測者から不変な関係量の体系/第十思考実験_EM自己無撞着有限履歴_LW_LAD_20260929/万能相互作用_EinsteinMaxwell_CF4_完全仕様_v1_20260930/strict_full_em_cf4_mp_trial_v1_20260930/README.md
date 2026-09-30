# Strict full Einstein-Maxwell CF4 MP trial v1

Reproduce with:

```bash
python3 run_strict_em_cf4_mp_v1.py
```

The selected case is the exact two-center extremal Majumdar-Papapetrou equilibrium with M_a=M_b=Q_a=Q_b=0.5 and coordinate separation 6. The dynamics use the full first-order GH + Maxwell state and the specified CF4 exp(B)exp(A) update.

This run is a validation trial, not a physical result: the Maxwell Gauss residual at the present N_cut=275 is too large, so the run is classified NUMERICALLY_UNRESOLVED. All raw outputs are retained.
