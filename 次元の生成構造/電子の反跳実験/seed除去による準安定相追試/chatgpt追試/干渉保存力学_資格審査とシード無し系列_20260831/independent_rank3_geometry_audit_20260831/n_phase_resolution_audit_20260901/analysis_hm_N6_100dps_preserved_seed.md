# hm_N6 100-digit evolution audit (corrected)

- N = 6
- steps = 500
- delta = 2*pi/6
- initial condition = original `data/hm_N6/parent_v.npz` complex128 values, promoted exactly from the stored binary64 values; no regeneration and no reseeding
- p,q observation plane = computed in float64 exactly as original `pass2_run.py`, then promoted exactly
- interaction = H_ef = A_ef conj(z_e) z_f, diagonal 0
- update = z_next = exp(-i delta H(z)) z, H frozen within each step
- arithmetic = mpmath 100 decimal digits; Hermitian eigendecomposition (`eighe`) and exponential phase action

Key values:

- step 0: Hperp/H = 3.6386025402876e-32
- step 1: Hperp/H = 3.64056958269e-32
- step 50: Hperp/H = 1.15681829138e-21
- step 100: Hperp/H = 9.04394493504e-8
- step 200: Hperp/H = 0.162228927626
- step 250: Hperp/H = 0.162229123447
- step 500: Hperp/H = 0.162229123600 (approximately)

The corrected high-precision evolution preserves the original float64 seed floor near 1e-32. Unlike the float64 500-step delta=2*pi/N run, the second rise from the ~0.16 plateau to the ~O(1) plateau does not occur by step 500. This isolates subsequent floating-point noise injection as a possible trigger for the second transition; this is an observation, not yet a proof.
