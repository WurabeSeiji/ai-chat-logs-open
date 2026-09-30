# Strict full Einstein-Maxwell CF4 trial v1

Case: two-center extremal Majumdar-Papapetrou, Q/M=1 for both centers.

This is not a reduced force model. The evolved state contains g_ab, Pi_ab, Phi_iab, H_a, Theta_a, E^i, B^i on the finite cardinal basis. The update is the specified state-dependent K(psi) and fourth-order commutator-free product exp(B)exp(A).

## Numerical audit
- nodes: 275
- K(psi)psi - F_direct L2: 0.000000e+00
- K(psi)psi - F_direct Linf: 0.000000e+00
- initial orbit radius readout: 3
- final orbit radius readout: 3
- relative orbit-radius drift: 0.000000e+00
- max |GW out-abs|: 1.324879e-12
- max |EM out-abs|: 1.982163e-06
- max raw Psi4 RMS: 1.582593e-01

The exact continuum MP solution is static. Therefore all nonzero drift and flux in this short run measure numerical truncation/discretization error of the trial implementation. They must not be interpreted as physical radiation.
