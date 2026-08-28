# Complex-simplex decompactification test — N=5 and N=16

The dynamics were not altered. Both runs use the no-K/sigma-normalization engine, seedless parent state Z0=v, the same GAMMA, and 5000 steps.

## Geometric readout

For each state z_ij, set d_ij^2=z_ij^2 and form the centered complex symmetric Gram matrix B=-1/2 J D^2 J. Its Takagi values are the singular values s_k(B); canonical simplex axis scales are r_k=sqrt(s_k).

The state is also split relative to the initial parent plane as Z=Z_parallel+Z_perp. The same complex-simplex reconstruction is applied separately to Z_perp. This is a readout only; it does not feed back into the dynamics.

## N=5
- full simplex rank: {'4': 5001}
- max |Z^T Z|: 5.385e-14
- H_total range: 0.8596491228067914 .. 0.8596491228071432
- A_perp: 1.331e-16 -> 0.000000 (max 0.000000)
- R_perp(Takagi): 8.154e-17 -> 0.000000
- R_perp early log growth rate: 0.000119/step
- growing perpendicular canonical axes measured: 4 of 4
- axis growth-rate range: 0.000042 .. 0.000114/step

## N=16
- full simplex rank: {'15': 5001}
- max |Z^T Z|: 2.094e-13
- H_total range: 3.6350618310580534 .. 3.6350618310696481
- A_perp: 2.462e-16 -> 0.000000 (max 0.000000)
- R_perp(Takagi): 1.198e-16 -> 0.000000
- R_perp early log growth rate: 0.000374/step
- growing perpendicular canonical axes measured: 15 of 15
- axis growth-rate range: 0.000029 .. 0.000341/step

## Interpretation constrained by the data

The perpendicular component is not merely a normalized plotting ratio: its raw amplitude A_perp grows exponentially from numerical-noise scale to O(1), while H_total remains conserved.

When that perpendicular component is read as a complex simplex, its canonical distance scales also grow exponentially. Therefore the decompactification-like reading survives a direct complex-distance reconstruction.

However, this particular Takagi-axis readout does not select only three expanding axes: N=5 shows all 4 non-null centered-simplex axes growing, and N=16 shows all 15. Thus the earlier 'three readable directions' phenomenon is not identical to the number of Takagi axes of the full complex simplex; it must be a separate readout/rank-selection structure if the two are to be connected.