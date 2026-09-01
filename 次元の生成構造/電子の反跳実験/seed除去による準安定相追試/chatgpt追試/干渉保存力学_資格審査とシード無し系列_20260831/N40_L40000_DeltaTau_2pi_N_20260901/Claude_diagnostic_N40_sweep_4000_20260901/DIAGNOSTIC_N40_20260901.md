# N=40, Delta tau=2pi/N diagnostic execution

- N: 40
- Requested L: 40000
- Exact update: original `np.linalg.eigh(H)` frozen-H step, no clipping/renormalization.
- Initial state: hm_N40 `parent_v.npz`, r^2=1/15.
- Completed diagnostic range in this environment: 0..150.
- First Hperp/H > 0.05: **step 109**.
- H_total: 51.999999999999986 -> 52.00000000000113.
- Hperp/H at step 100: 0.0085609679875533352.
- Hperp/H at step 150: 0.95674755978805104.
- global closure at step 100: 0.99966319878140686.

## Interpretation
The inflation threshold is already crossed at step 109; therefore the new high-symmetry N=40 initializer does not suppress the numerical-seed inflation in this binary64/eigh implementation over the early-time window. The norm remains conserved to roundoff through step 150.

## Important status
This is not presented as the requested completed L=40000 trajectory. A full 40000-step run with a 780x780 Hermitian eigendecomposition at every step exceeds the single-turn execution budget here. No dynamics-changing acceleration has been substituted.
