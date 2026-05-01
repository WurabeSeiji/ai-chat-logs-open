# Geometric Observations on the Second-Order Correction of the α Identity — Supplement to Paper 7

Following the companion papers (Paper 7 and Paper 8):

　　α⁻¹ = 137 + (π²/2)·α (8.7 ppb precision)
　　α-identity ≅ Wilson Lattice Gauge Theory (structural correspondence on a 4D chain complex)

I have published a supplement on Zenodo concerning the geometric interpretation of the residual 8.7 ppb left in Paper 7.

This is a formal supplement to Paper 7, but it is **an observation paper, NOT a proof paper**. It is not intended to "strengthen" or "extend" the claims of Paper 7.

Publication information:

- DOI (Concept, always resolves to latest): https://doi.org/10.5281/zenodo.19933729
- DOI (v1, latest): https://doi.org/10.5281/zenodo.19933730
- Zenodo page: https://zenodo.org/records/19933730
- License: CC BY 4.0
- Format: md / tex / pdf × JA/EN = 6 files

---

## The Open Problem: Origin of the 8.7 ppb Residual

The α identity in Paper 7

　　α⁻¹ = 137 + (π²/2)·α

predicts α at 8.7 ppb precision. This is approximately 1/3000 of the deviation from Eddington-style integer fitting (2.6×10⁻⁴), strikingly precise. Yet a residual of 6.3×10⁻¹⁰ remains. Paper 7 left this residual as an open problem, identifying it as a candidate for an α³-order geometric correction with coefficient c₃ ≈ 1.6×10⁻³.

This supplement provides **one geometric interpretation** of this residual.

---

## Partial Response to Anticipated Concerns

The structure of Paper 7 may be subject to such concerns as:

- Concern A: "π²/2 as the unit ball volume may have been chosen post-hoc to match observation."
- Concern B: "Without a geometric origin for the residual, the coefficient π²/2 cannot be excluded as numerical fitting."

This supplement offers a **partial response** to these concerns. It is not a complete refutation but an observation that one geometric interpretation is consistent with the structure of Paper 7.

---

## Five Observations

A summary of what was recorded chronologically. For details, refer to the paper PDF on Zenodo.

### Observation 1: An Alternative Reading of π²/2

In Paper 7, π²/2 was introduced as the volume of the unit 4-ball V₄(1). However, π²/2 also appears naturally from a different route.

The surface measure of the 4D ball boundary S³:

　　A_{S³}(R) = 2π²R³

Normalized by 4R³:

　　A_{S³}(R) / (4R³) = π²/2

This means π²/2 admits a reading as "**the normalization core of the boundary measure of the 4D ball, normalized by R³**" — yielding a **scale-invariant value** for any radius R.

Numerically, V₄(1) = π²/2 is identical, but as interpretations:
- V₄(1): includes "the choice of unit"
- A_{S³}/(4R³): scale-invariant, independent of R

This is the presentation of one reading; we do not claim it is "more correct."

### Observation 2: Asymptotic Expansion of the Volume Deficit G(R)

On the 4D Euclidean integer lattice, place a unit hypercube cell of side 1 centered at each lattice point. Let N_full(R) denote the number of cells fully contained in the 4D ball of radius R. Define the volume deficit G(R) as:

　　G(R) := V₄(R) − N_full(R) = (π²/2)R⁴ − N_full(R)

This corresponds to the contribution of cells partially cut by the spherical surface.

Using the outward support function of a unit hypercube with respect to the unit normal u on the spherical surface:

　　τ(u) = (1/2)·(|u₁| + |u₂| + |u₃| + |u₄|)

G(R) admits the asymptotic expansion:

　　G(R) = (16π/3)R³ − (3π²/4 + 9π/2)R² + O(R)

where:
- First term (16π/3)R³: main contribution from the continuous ball boundary
- Second term −(3π²/4 + 9π/2)R²: second-order effect of "shaving"

Closed forms of the sphere integrals:

　　∫_{S³} τ(u) dΩ = 16π/3
　　∫_{S³} τ(u)² dΩ = π²/2 + 3π

These can be derived rigorously (see §3 of the supplement paper).

### Observation 3: The Sign of the Second Term Is Necessarily Negative

The coefficient of the second term is:

　　−(3π²/4 + 9π/2) ≈ −7.402 − 14.137 ≈ −21.539

**Under this approximation model, this is necessarily negative** (since ∫τ² dΩ > 0 with the −3/2 coefficient in front).

This means, under the approximation "the spherical boundary is shaved inward by the hypercube structure," the continuous ball approximation **overestimates** the volume deficit.

### Observation 4: Numerical Verification

For integer R = 50, 100, 200, 300, computing N_full(R) and G(R) yields:

- R = 50: G(R)/R³ ≈ 16.31 (theoretical 16π/3 ≈ 16.755), second coefficient ≈ −22.11
- R = 100: G(R)/R³ ≈ 16.55, second coefficient ≈ −20.22
- R = 200: G(R)/R³ ≈ 16.66, second coefficient ≈ −19.13
- R = 300: G(R)/R³ ≈ 16.68, second coefficient ≈ −21.42

The second coefficient fluctuates around the theoretical value −21.539 (the fluctuation is consistent with finite-size effects of the discrete lattice).

The numerical experiment is consistent with the asymptotic expansion. However, consistency does not show "this interpretation is correct"; it only shows "**this interpretation is not numerically inconsistent**."

### Observation 5: Sign and Magnitude Consistency with Observed Residual

Introduce a correction term δ in the equation of Paper 7:

　　α⁻¹ = 137 + (π²/2 − δ)·α

Back-calculating from CODATA 2022 value α⁻¹ = 137.035999177:

　　δ = π²/2 − α⁻¹·(α⁻¹ − 137) ≈ 0.001619

Observed facts:
- The coefficient of the second term, −21.539, is **negative**
- The δ back-calculated from observation is also **subtracted from π²/2**
- Their signs are consistent
- δ/(π²/2) ≈ 3.3×10⁻⁴ falls within the natural order-of-magnitude range for discrete corrections of hypercube lattices

However, **the complete derivation of δ is not achieved**. The normalization rule converting the R² deficit term to the coefficient term in the α equation is not formulated in this supplement.

---

## Important Caveats

This supplement does NOT claim:

- That π²/2 admits no interpretation other than the one presented here
- That the residual of 8.7 ppb is fully derived under our interpretation
- That other interpretive possibilities are excluded
- That this supplement "proves," "strengthens," or "extends" the claims of Paper 7

What this supplement presents is **one reading**. The possibility that the same equations and numbers admit other readings is not denied.

---

## Implications for Paper 7

### Positioning as Reinforcement

The contribution of this supplement is limited. This supplement does NOT:

- **Modify** the numerical fact of 8.7 ppb precision in Paper 7
- **Justify** the choice of π²/2 as a coefficient (it merely presents an alternative reading of the same numerical value)
- **Solve** the residual problem (no complete derivation of δ is provided)

However, this supplement does provide:

- A geometric viewpoint that π²/2 admits reading as a scale-invariant boundary measure core
- The observation that, under this viewpoint, the sign and magnitude of the residual are interpretable
- A partial response to the concerns about Paper 7

### Effect on Defense Lines

This supplement does **NOT depend on** the claims of Paper 7:

- The 8.7 ppb observation of Paper 7 does not change regardless of acceptance or rejection of this supplement's interpretation.
- Even if this supplement is partially rejected, Paper 7 remains unaffected.
- This supplement is a reinforcement attempt and does not modify the logical foundation of Paper 7.

---

## Open Questions Remaining

Gap G1: A normalization rule converting the R² deficit term into the coefficient term δ of the α equation is not formulated.

The second term in the asymptotic expansion is an R²-order volume correction, while δ is a dimensionless coefficient in the α equation. To connect the two, a mapping rule is needed between volume corrections at the specific scale R = 3 (where Paper 7 yields 137) and the self-consistency condition of the α equation. This supplement does not provide this rule.

Other open questions:
- Analytical computation of higher-order terms (O(R), O(1))
- Possibility of deriving an isomorphic correction term from other geometric structures (4D cross-polytope, root system lattices, etc.)
- Re-interpretation in the framework of Wilson structural correspondence

---

## Related Papers (Zenodo, CC BY 4.0)

- Paper 7 (α identity, 8.7 ppb precision): https://doi.org/10.5281/zenodo.19869266
- Paper 8 (α isomorphism, Wilson structural correspondence): https://doi.org/10.5281/zenodo.19880467
- **Paper 7 Supplement (this paper, observation paper)**: https://doi.org/10.5281/zenodo.19933729
- Working Paper (Thought Experiment Note on 6D extension): https://doi.org/10.5281/zenodo.19902677

## Related Zenn Articles

- α identity (Paper 7): https://zenn.dev/noriaki_kihara/articles/alpha-identity-4d-geometry
- α isomorphism (Paper 8, Wilson categorical correspondence): https://zenn.dev/noriaki_kihara/articles/alpha-isomorphism-lattice-gauge
- BH Thermodynamics Programme (six core papers): https://zenn.dev/noriaki_kihara/articles/bh-thermodynamics-projection
- Working Paper (6D extension thought experiment): https://zenn.dev/noriaki_kihara/articles/rethinking-xyztRQ-thought-experiment

## Related Note Articles

- α identity (Paper 7, Japanese): https://note.com/kiharanoriaki/n/n19cc13927c51
- α identity (Paper 7, English): https://note.com/kiharanoriaki/n/n01d98237ddae
- α isomorphism (Paper 8, Japanese): https://note.com/kiharanoriaki/n/n87df7a4977e7
- α isomorphism (Paper 8, English): https://note.com/kiharanoriaki/n/n9c37c6f99a0a
- Working Paper (Japanese): https://note.com/kiharanoriaki/n/nc75457ca96f7
- Working Paper (English): https://note.com/kiharanoriaki/n/n39fecf32607f

---

## Conclusion

When future researchers, perhaps a hundred years from now, provide a more complete geometric explanation of the structure of Paper 7, this supplement may be referenced as an intermediate observation. That is sufficient.

This research programme will continue to be carefully built up as an independent researcher's work.

---

Author: Noriaki Kihara
WF System Co., Ltd. / ORCID: 0009-0004-6753-4020

Paper DOI (Concept): https://doi.org/10.5281/zenodo.19933729
Paper DOI (v1, latest): https://doi.org/10.5281/zenodo.19933730
Zenodo page: https://zenodo.org/records/19933730

---

#FineStructureConstant #137 #Alpha #LatticeGaugeTheory #Wilson #FourDimensionalGeometry #AsymptoticExpansion #ObservationPaper #Supplement #Zenodo #TheoreticalPhysics #MathematicalPhysics #BlackHoleThermodynamics #Preprint
