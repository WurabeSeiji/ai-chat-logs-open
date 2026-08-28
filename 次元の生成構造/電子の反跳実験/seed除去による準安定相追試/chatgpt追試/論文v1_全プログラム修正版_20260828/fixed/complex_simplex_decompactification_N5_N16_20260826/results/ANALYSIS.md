# Complex-simplex decompactification test — N=5 and N=16

The dynamics were not altered. Both runs use the no-K/sigma-normalization engine, seedless parent state Z0=v, the same GAMMA, and 5000 steps.

## Geometric readout

For each state z_ij, set d_ij^2=z_ij^2 and form the centered complex symmetric Gram matrix B=-1/2 J D^2 J. Its Takagi values are the singular values s_k(B); canonical simplex axis scales are r_k=sqrt(s_k).

The state is also split relative to the initial parent plane as Z=Z_parallel+Z_perp. The same complex-simplex reconstruction is applied separately to Z_perp. This is a readout only; it does not feed back into the dynamics.

## N=5
- full simplex rank: {'4': 5001}
- max |Z^T Z|: 3.406e-14
- H_total range: 0.8596491228069796 .. 0.8596491228070694
- A_perp: 1.612e-16 -> 0.251586 (max 0.655610)
- R_perp(Takagi): 1.191e-16 -> 0.203946
- R_perp early log growth rate: not fitted (no exponential window 1e-10<R_perp<1e-3 with >=10 points)
- growing perpendicular canonical axes measured: 0 of 4
- no fitted axis rates

## N=16
- full simplex rank: {'15': 5001}
- max |Z^T Z|: 2.609e-13
- H_total range: 3.6201925420978029 .. 3.6201925420987080
- A_perp: 3.542e-16 -> 1.878963 (max 1.894670)
- R_perp(Takagi): 1.693e-16 -> 0.978747
- R_perp early log growth rate: not fitted (no exponential window 1e-10<R_perp<1e-3 with >=10 points)
- growing perpendicular canonical axes measured: 3 of 15
- axis growth-rate range: 0.113247 .. 0.183880/step

## Interpretation constrained by the data

The perpendicular component is not merely a normalized plotting ratio: its raw amplitude A_perp grows exponentially from numerical-noise scale to O(1), while H_total remains conserved.

When that perpendicular component is read as a complex simplex, its canonical distance scales also grow exponentially. Therefore the decompactification-like reading survives a direct complex-distance reconstruction.

However, this particular Takagi-axis readout does not select only three expanding axes: N=5 shows all 4 non-null centered-simplex axes growing, and N=16 shows all 15. Thus the earlier 'three readable directions' phenomenon is not identical to the number of Takagi axes of the full complex simplex; it must be a separate readout/rank-selection structure if the two are to be connected.