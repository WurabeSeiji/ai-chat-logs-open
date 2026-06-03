# Realigning Uncertainty, Sp(2,ℝ), Wick Rotation, and Stone's Theorem Through the Single Lens of "Phase-Space Area" — Paper 11 (Observation Paper)

Classical signal theory (Fourier analysis, sampling) and quantum theory, when read as structures on the phase plane (q, p), share a single **common currency: area**. Both the action ∮ p dq (a conserved quantity) and the area spanned by the uncertainty (a fluctuation) are measured as the **symplectic area** spanned by a conjugate pair.

Using only this single lens — "area as the common currency" — I have published on Zenodo an **observation paper** that realigns structures already present in standard mathematics: the Heisenberg–Gabor uncertainty relation, the Robertson inequality, the Sp(2,ℝ) ≅ SU(1,1) symmetry, Stone's theorem, and Wick rotation.

This paper does NOT modify standard quantum theory, does NOT deny Wick rotation, and changes NO observable prediction. It does not derive the metric signature, the imaginary unit, or the complex structure. It is merely a record of re-reading known structures through one lens.

Publication information:

- DOI (Concept, always resolves to latest): https://doi.org/10.5281/zenodo.20521566
- DOI (v1.0, latest): https://doi.org/10.5281/zenodo.20521567
- Zenodo page: https://zenodo.org/records/20521567
- License: CC BY 4.0
- Format: md / tex / pdf × JA/EN + 2 figures = 8 files

---

## Starting Point: Area as a Common Currency

On the phase plane (q, p), two seemingly unrelated quantities are both measured in units of the same "area":

・**Conserved quantity**: the action ∮ p dq along a closed orbit is the area enclosed by that orbit. Bohr–Sommerfeld quantization ∮ p dq = (n + ½) h is the discretization of this area.
・**Fluctuation**: the uncertainty Δq Δp ≥ ℏ/2 is the minimal area a state occupies in phase space. The √detΣ of the covariance matrix Σ is that area.

The two are measured in the same **unit** of area, but their numerical values do not coincide (ℏ/2, πℏ, 2πℏ appear with coefficients differing by 2π / 4π). Keeping this distinction intact, the only operation of this paper is to realign known structures using area as a common language.

---

## Three Things Observed

### (1) The ½ that appears in three settings, and the metaplectic representation

The number ½ appears in several settings that look mutually unrelated:

・the ½ in semiclassical quantization ∮ p dq = (n + ½) h
・the equality condition of the Robertson inequality (the Gaussian ground state)
・the Maslov index (a phase counting orbit caustics)
・the SU(1,1) lowest weight (the Bargmann index)

These are **not numerically the same object**. The Maslov index is system-dependent, the floor ½ carries the dimension of action, and the weight is dimensionless — they must not be conflated. Yet they are not mere coincidence either: each can be understood as one of three cross-sections of a structure connected to the **half-integer weight** of the double cover of Sp(2,ℝ) — the metaplectic (Weil) representation.

This paper does not claim "therefore ½ originates from X." It stays at the observation that they all touch the same covering structure.

### (2) Boost = non-compact subgroup = squeeze

The full set of area-preserving linear transformations is the group Sp(2,ℝ) ≅ SU(1,1). The symplectic area √detΣ is kept invariant under this group action, and squeezing, chirping, and rotation are all actions of this group.

Of these, the **non-compact subgroup (squeezing)** is isomorphic to a Lorentz boost on phase space:

　　v/c = tanh η　　(rapidity parameter)
　　k = e^(2η)　　(shape parameter)

This holds rigorously as the action of the standard squeeze operator S(η) = exp( (η/2)(a² − a†²) ) on the canonical variables. It is the correspondence by which the squeezed states used routinely in quantum optics can be read as "boosts" in phase space.

### (3) A generator-level re-reading of the imaginary unit i and the hyperbolic structure

Whether a two-dimensional surface is circular (metric signature +) or hyperbolic (signature −) is determined by whether the generator acting on that surface is of type

　　J² = −I (rotation, compact)
　　K² = +I (boost, non-compact)

and the two are carried into each other by

　　θ = i η　　(the standard form of Wick rotation).

Stone's theorem gives the generator of a continuous unitary evolution the form G = iH (anti-self-adjoint).

In this paper's view, the i appearing here is **part of a presupposed structure** — the almost-complex structure (the Kähler / compatible triple) fixed by the compatibility of the area ω and the metric g. This paper does NOT derive this i. It merely **re-reads** the very same i identified by Stone's theorem, in a form observationally equivalent to Wick rotation and standard quantum theory.

Note also that because of the uncertainty floor ½, there is no point-like state of zero area in phase space.

---

## On the Figures

Two figures accompany this paper:

・Figure 1: how a squeeze transformation deforms a circle in phase space into an ellipse while preserving the area (√detΣ). It visualizes the squeeze = boost correspondence.
・Figure 2: how the sign of the same generator (J² = −I vs. K² = +I) yields different projected images — compact (circle) vs. hyperbolic (open hyperbola).

Both are illustrations of standard content and contain no new claims.

---

## What This Paper Does NOT Claim (Explicit)

- any change to the predictions of standard quantum theory or special relativity
- any claim that Wick rotation or the Minkowski metric is wrong
- any derivation of the metric signature, the imaginary unit, or the complex structure
- any prediction of a new physical constant, cross section, or decay rate, or proof of a new mathematical theorem

What this paper offers is one way of reading known structures (a juxtaposition under the common currency of area). Evaluation and interpretation are left to the reader. It is published as an observation record that, having passed multiple rounds of review by four AIs (Claude.ai / ChatGPT / Gemini / Grok), removes overclaims and states its limits explicitly.

---

## Related Resources

- Zenn article (more technical): https://zenn.dev/noriaki_kihara/articles/phase-space-area-symplectic-uncertainty
- Structural correspondences between signal/control theory and quantum theory (Paper 10, companion): https://doi.org/10.5281/zenodo.20521598
- note index article (all papers): https://note.com/kiharanoriaki/n/nc1619291b690

---

Author: Noriaki Kihara
WF System Co., Ltd. / ORCID: 0009-0004-6753-4020

Paper DOI (Concept): https://doi.org/10.5281/zenodo.20521566
Paper DOI (v1.0, latest): https://doi.org/10.5281/zenodo.20521567
Zenodo page: https://zenodo.org/records/20521567

---

#PhaseSpace #SymplecticGeometry #UncertaintyPrinciple #Sp2R #SU11 #WickRotation #StonesTheorem #SqueezedStates #QuantumOptics #MetaplecticRepresentation #PhaseSpaceArea #ObservationPaper #Zenodo #TheoreticalPhysics #MathematicalPhysics #Preprint
