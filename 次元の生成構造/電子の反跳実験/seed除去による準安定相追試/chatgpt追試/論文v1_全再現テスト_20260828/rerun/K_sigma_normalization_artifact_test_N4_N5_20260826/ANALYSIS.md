# K/sigma normalization artifact test — N=4 and N=5

## Purpose
Test whether the exponential growth of the parent-plane-orthogonal component is caused by the physical normalization K -> K/sigma_max.

## Controlled comparison
Both branches use exactly the same:
- make_parent implementation
- initial parent v
- zero-closure seed g
- DELTA = 1e-12
- initial Z0
- GAMMA = tan(pi/144)
- initial power-iteration probe vector wp0
- 5000 steps
- observables and fitting windows

The only dynamics difference is:

Normalized:
    gn = GAMMA / sigma
    A2 = (sigma / GAMMA) * J + G

Raw K:
    gn = GAMMA
    A2 = (1 / GAMMA) * J + G

## Result

N=4:
- H_perp exponential rate normalized = 0.054884390/step
- H_perp exponential rate raw K      = 0.148644695/step
- A_perp exponential rate normalized = 0.027442195/step
- A_perp exponential rate raw K      = 0.074322347/step
- H_perp reaches 1e-8 at step 708 normalized vs 262 raw K.

N=5:
- H_perp exponential rate normalized = 0.049359169/step
- H_perp exponential rate raw K      = 0.172506395/step
- A_perp exponential rate normalized = 0.024679585/step
- A_perp exponential rate raw K      = 0.086253197/step
- H_perp reaches 1e-8 at step 757 normalized vs 217 raw K.

For both N=4 and N=5:
- H_total stays at 1 to floating-point precision in both branches.
- |Z^T Z| stays near machine precision.
- Removing K/sigma does NOT remove exponential growth.
- Removing K/sigma makes the early exponential growth substantially faster.

## Conclusion
The observed exponential amplification of the orthogonal amplitude is not created by K/sigma normalization.
Within this controlled low-N test, K/sigma normalization suppresses the growth rate rather than producing it.

The physically meaningful statement supported here is:
A_perp = ||Z_perp|| grows exponentially from a microscopic seed while the total Hermitian norm is conserved.
Therefore H_perp = A_perp^2 also grows exponentially, with twice the logarithmic rate of A_perp.

This test does not by itself establish a cosmological interpretation. It establishes that the exponential instability survives removal of the sigma normalization and is therefore not an artifact of that normalization.
