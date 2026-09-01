# 64-bit denominator control with fixed 124 baseline

## Conditions

- N = 3..16
- 500 steps
- state dtype: complex128
- real dtype: float64
- initial state: original saved `states_treatment.npz["Z"][0]` for each hm_N
- update: `H=A*(conj(z)[:,None]*z[None,:])`, Hermitian eigendecomposition with `numpy.linalg.eigh`, then `exp(-i*(2*pi/denominator)*w)`
- denominators: N-2, N-1, N, N+1, N+2, and fixed 124

## Purpose

This rerun regenerates all dynamics in 64-bit arithmetic and adds the original fixed denominator 124 as a sixth curve in every N panel. No initial state is regenerated.

## Selected onset (>0.05) steps

| N | N-2 | N-1 | N | N+1 | N+2 | 124 |
|---:|---:|---:|---:|---:|---:|---:|
| 3 | 92 | 259 | -1 | -1 | -1 | -1 |
| 4 | 51 | 79 | 121 | 163 | 233 | -1 |
| 5 | 134 | 172 | 207 | 257 | 299 | -1 |
| 6 | 74 | 94 | 120 | 143 | 180 | -1 |
| 7 | 87 | 114 | 138 | 170 | 203 | -1 |
| 8 | 106 | 134 | 159 | 187 | 224 | -1 |
| 9 | 129 | 160 | 190 | 218 | 259 | -1 |
| 10 | 135 | 162 | 188 | 220 | 251 | -1 |
| 11 | 160 | 185 | 212 | 247 | 284 | -1 |
| 12 | 185 | 215 | 249 | 285 | 325 | -1 |
| 13 | 206 | 243 | 274 | 311 | 354 | -1 |
| 14 | 236 | 273 | 311 | 348 | 385 | -1 |
| 15 | 268 | 311 | 341 | 386 | 422 | -1 |
| 16 | 300 | 333 | 376 | 435 | 461 | -1 |

The fixed 124 curve remains near the numerical floor over 500 steps for all N, while the N±2 family reaches macroscopic Hperp/H within the 500-step window for most N. The detailed interpretation is intentionally left separate from this reproduction record.
