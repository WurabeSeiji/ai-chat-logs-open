# Dual Geometry of Wavelength Space and Frequency Space: Reciprocal Duality and 4D Lattice Counting (Observational Model)

I have published two short observational papers as a new series on Zenodo. They place the simplest reciprocal duality λₙ = 1/νₙ between wavelength components λₙ and frequency components νₙ, impose a constant sum-of-squares (norm) condition on each of the wavelength space and the frequency space, and then observe—using only elementary arguments—what degrees of freedom, symmetries, and counting structures arise.

These papers do not derive any physical entity such as spacetime, mass, energy, or momentum. Here λₙ and νₙ are model variables for describing the dual geometry; νₙ is not identified with physical time-frequency, energy, or momentum. No new prediction or theorem is offered. Everything stays at the level of observing geometric and topological structure. Evaluation and interpretation are left to the reader.

Publication details:

- Paper 1: Dual Geometry of Wavelength Space and Frequency Space (v0.3)
- Zenodo: https://zenodo.org/records/20588037
- Concept DOI: https://doi.org/10.5281/zenodo.20588036
- Version DOI: https://doi.org/10.5281/zenodo.20588037

- Paper 2: Radius Sweep of Fully-Inscribed Unit-Cell Counts on a 4-Dimensional Lattice (v0.1)
- Zenodo: https://zenodo.org/records/20588039
- Concept DOI: https://doi.org/10.5281/zenodo.20588038
- Version DOI: https://doi.org/10.5281/zenodo.20588039

- Paper 3 (supplement to Paper 1): Closed Four-Degree-of-Freedom Structure and Its Correspondence with 4-Dimensional Lattice Counting (v0.3)
- Zenodo: https://zenodo.org/records/20589515
- Concept DOI: https://doi.org/10.5281/zenodo.20589261
- Version DOI: https://doi.org/10.5281/zenodo.20589515

- License: CC BY 4.0
- Format: md / tex / pdf in Japanese and English, plus figures / CSV

---

## 1. Setup — reciprocal duality and norm preservation

We place a componentwise reciprocal duality between wavelength and frequency components.

　λₙ = 1 / νₙ

We then impose a constant sum-of-squares condition on each space.

　Σ λₙ² = Λ²　(norm radius Λ of the wavelength space)
　Σ νₙ² = N²　(norm radius N of the frequency space)

When one component is large, the corresponding other component is small. This is not a symmetry of identical form but a dual symmetry that includes inversion.

The following figure is a schematic of this dual constraint.

![Figure 1: Schematic of the λ–ν dual constraint](./figure1_lambda_nu_dual_constraint_EN.png)

▲ Figure 1: Schematic of the λ–ν dual constraint. In one dimension the solution is fixed to essentially a single point; a 5-component / 1-constraint system retains four degrees of freedom. Introducing uncertainty, it can be represented as an observational region with thickness around the ideal constraint (the magnitude of that thickness is not derived here).

## 2. One dimension is almost a point; four degrees of freedom at 5 components / 1 constraint

In one dimension the positive solution is fixed to

　λ = Λ,　ν = 1/Λ　(with ΛN = 1)

leaving almost no room for redistribution among components.

By contrast, imposing one sum-of-squares constraint on five components leaves

　5 − 1 = 4

degrees of freedom. These four degrees of freedom are an important object of observation, but we do not immediately identify them with the 4-dimensionality of physical spacetime (they are degrees of freedom that appear purely on the dual geometry).

The existence condition, from a Cauchy–Schwarz type inequality, is

　ΛN ≥ 5　(for general d components, ΛN ≥ d)

with equality when all components are equal.

## 3. Logarithmic representation — the hyperbola becomes a sign-reversal symmetry

In the ordinary representation the reciprocal duality is the hyperbola λν = 1. But introducing the logarithmic variables

　qₙ = log λₙ,　pₙ = log νₙ

since log λₙ = −log νₙ, we obtain

　pₙ = −qₙ

a straight line of slope −1, that is, a sign-reversal symmetry.

![Figure 2: Logarithmic representation of the λ–ν duality](./figure2_log_representation_lambda_nu_duality_EN.png)

▲ Figure 2: Logarithmic representation of the λ–ν duality. The hyperbola λν = 1 of the ordinary representation becomes the straight line p = −q in the logarithmic variables q = log λ, p = log ν. In an observed state with uncertainty, it can be represented as a band with thickness around the ideal line (that thickness is not determined here).

## 4. Counting on a 4D lattice — N₀(3) = 137

Choosing the frequency space or the wavelength space to be the 4-dimensional integer lattice ℤ⁴ turns the sum-of-squares condition into a unit-cell counting problem.

The condition that a unit cell of side 1 (half-width ½ in each direction) is fully inscribed in the 4-dimensional hyperball of radius R is

　Σ (|kᵢ| + ½)² ≤ R²

Defining this fully-inscribed cell count as N₀(R), at integer radii

　N₀(1) = 1,　N₀(2) = 9,　N₀(3) = 137

Paper 2 sweeps R from 0.5 to 10.0 in steps of 0.5 and tabulates, together with the circumscribing diameter 2R, the stacked diagonal length 2ρ(R), and the radial margin R − ρ(R) (a reproducible pseudocode and a CSV are included).

### What the 137 at R = 3 is made of

Multiplying both sides by 4 gives Σ (2|kᵢ| + 1)² ≤ 36, and the counts per shell are:

- shell 4: 1
- shell 12: 8
- shell 20: 24
- shell 28: 40
- shell 36: 64

Total: 1 + 8 + 24 + 40 + 64 = 137

This 137 is not introduced under any assumed correspondence with a physical constant. It is a pure counting result from the 4-dimensional integer lattice, the unit cell, and the full-inscription condition.

## 5. Uncertainty is placed only as a "form"

Placing an integer lattice with a fractional fluctuation νₙ = mₙ + εₙ on the frequency side induces, on the wavelength side, an opposite-sign fluctuation

　Δλₙ ≈ − εₙ / mₙ²

Replacing the full-inscription indicator with a weight function Wδ having thickness near the boundary, the effective count Nδ(R) tends to N₀(R) in the limit δ → 0.

However, this work derives neither the value of the uncertainty width δ nor the explicit form of Wδ. It only defines the counting formalism.

## 6. Supplement (Paper 3): correspondence between the 5-component constraint and the 4D lattice counting

The 5-component sum-of-squares constraint of Paper 1, Σ xₙ² = R², defines a 4-dimensional hypersphere S⁴_R (a closed object with four degrees of freedom) inside a 5-dimensional space. The counting of Paper 2, on the other hand, is carried out inside the 4-dimensional ball B⁴_R (Σ uᵢ² ≤ R²). Paper 3 is a short supplement that organizes this "5-component constraint → four DOF → 4D lattice counting" correspondence, with no physical interpretation added.

There is a single point. For a point that already satisfies the constraint (‖λ‖ = Λ), applying the radial projection of the same radius, Π_R(y) = R·y/‖y‖, gives

　λ′ = Λ·λ/‖λ‖ = λ,　ν′ = ν

i.e. the identity map. So this radial projection is not an operation that changes values, but merely a geometric rephrasing that reads a constraint-satisfying point as a point on a closed four-DOF structure of constant radius. Geodesic cell partitions on S⁴_R and area-based partitions are not treated in this supplement and are left to future work.

## 7. Positioning

- A consistent counting that simultaneously discretizes the frequency and wavelength spaces and satisfies the reciprocal duality is left to future work (the papers define only the single-space counting).
- The four degrees of freedom, p = −q, and N₀(3) = 137 are observations; no identification with physical quantities or physical constants is claimed.
- Related work is referred to only as mathematical background, not as physical grounds: symmetric treatment of time and frequency (Gabor 1946), communication theory (Shannon 1948), geometry of numbers (Aliev–Henk 2023), sum of four squares (Hirschhorn 1987).

---

Author: Noriaki Kihara
WF System Co., Ltd. / ORCID: 0009-0004-6753-4020

Zenn article (more technical, with rendered formulas): https://zenn.dev/noriaki_kihara/articles/wavelength-frequency-dual-geometry

---

#WavelengthSpace #FrequencySpace #ReciprocalDuality #FourierAnalysis #LogarithmicRepresentation #SignReversalSymmetry #4DLattice #LatticePointProblem #GeometryOfNumbers #UnitCellCounting #SumOfFourSquares #RadialProjection #Geometry #MathematicalPhysics #TheoreticalPhysics #ObservationalPaper #ThoughtExperiment #Preprint #Zenodo
