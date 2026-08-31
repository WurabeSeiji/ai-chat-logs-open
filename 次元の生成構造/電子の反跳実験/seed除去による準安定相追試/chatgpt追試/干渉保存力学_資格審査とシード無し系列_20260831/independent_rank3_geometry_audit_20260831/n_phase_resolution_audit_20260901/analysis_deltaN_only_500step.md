# deltaN-only 500-step control experiment

Only changed parameter from the existing dynamics: delta = 2*pi/124 -> 2*pi/N. N=3..16, 500 steps. Existing hm initial states, interaction H_ef=A_ef conj(z_e) z_f, eigendecomposition exponential update, no phase projection, no normalization.

## Main result
For N=4..16 the transverse fraction Hperp/H crosses 0.05 within 121--369 steps. N=3 does not cross 0.05 within 500 steps (final 3.59e-5). Norm drift remains at ~1e-13 or below, consistent with unitary frozen-H updates.

Onset step Hperp/H > 0.05:
N4 121; N5 207; N6 120; N7 138; N8 159; N9 190; N10 188; N11 212; N12 258; N13 274; N14 318; N15 355; N16 369.

By step 500, global/local square-closure is strongly broken for most N; this is expected from the current H law because square-closure is not an invariant of the update. This experiment does not change or test that interaction law; it isolates only the effect of replacing delta by 2*pi/N.

No regular-polygon reconstruction, phase-grid projection, re-normalization, seed change, or alternative interaction was used.
