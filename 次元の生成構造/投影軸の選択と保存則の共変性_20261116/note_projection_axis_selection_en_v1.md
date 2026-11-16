# Which Axes Are Spacetime and Which Are Internal? — Two Radii Carried by All Waves and the Covariance of Conservation Laws

In the previous article, I introduced a roadmap of thought experiments attempting to reconstruct spacetime and elementary particles solely from a bounded Hilbert space and the algebraic structures of quaternions and octonions, without assuming a background spacetime.

However, that previous paper was primarily an organizational map of open problems rather than a source of direct physical predictions. Following an immediate follow-up investigation into the "two R's" of quaternions, clear structural predictions emerged, which I have now published as a new preprint on Zenodo (Concept DOI: https://doi.org/10.5281/zenodo.24273723).

Here is an overview of the four core claims and predictions.

## 1. All Waves Simultaneously Reside on Both R₀ and R₁

The polar form of a quaternion

　o = A(cos θ + n sin θ)　(n ∈ S², n² = −1)

contains two radii of distinct character:

- R₀ (projection axis; radius of the 1D phase circle S¹): the radius of the circle on which the real part cos θ and the magnitude sin θ reside.
- R₁ (spatial side; curvature radius of the 2D directional sphere S² inside pure-imaginary 3D space Im(ℍ) ≅ ℝ³): the radius carrying direction and extension spanned by the unit imaginary vector n.

Initially, we imagined a division of labor where only a scalar cos wave (a Higgs-like uniform field) lives on R₀ while 3D space lives on R₁. Closer inspection showed that **there is no reason for only one type of wave to occupy R₀; rather, every wave in the Hilbert space simultaneously resides on both R₀ and R₁ from the outset**.

Because every wave carries both components, quaternion multiplication o₁ o₂ directly couples 1D and 3D: the 3D inner product n₁ · n₂ alters the scalar amplitude on R₀, while the phases on R₀ control the 3D twist (cross product n₁ × n₂).

## 2. Projection Asymmetry in the Hopf Fibration

When all waves reside on both R₀ and R₁, the mapping between them is asymmetric due to the Hopf fibration S¹ ↪ S³ → S²:

- Onto R₁ (the directional sphere S²), every point of S³ maps surjectively (full projection).
- Along R₀ (the 1D phase fiber S¹), rotating the global phase leaves the image on S² collapsed to a single point (zero projection as spatial displacement).
- The phase θ controls only the weight magnitude sin θ appearing on R₁ without fixing the direction n (partial projection).

Consequently, to an observer reading the world from R₁, the R₀ degrees of freedom appear not as spatial directions, but as a uniform background (scalar reference, phase, and temporal progression). Under the standard Riemannian submersion, geometry fixes the ratio to R₀ : R₁ = 2 : 1.

## 3. Selecting 1-, 2-, or 4-Dimensional Projection Axes

There is no reason to fix the dimension of R₀ to 1. If we choose a multiplicatively closed subalgebra as the projection axis R₀, Hurwitz's theorem (real normed division algebras exist only in dimensions 1, 2, 4, 8) restricts the choices inside the octonions 𝕆 (8D) to **1, 2, or 4 dimensions (ℝ, ℂ, ℍ)**:

- 1D axis (ℝ): opens a 3D space (ℍ) or 7D space (𝕆).
- 2D axis (ℂ): yields ℂP¹ ≅ S² in ℍ, and decomposes the remaining 6D of 𝕆 into ℂ³, where SU(3) (color symmetry) emerges as the subgroup of G₂ preserving that complex axis.
- 4D axis (ℍ): makes the remainder in 𝕆 the quaternionic projective line ℍP¹ ≅ S⁴, i.e., closed 4D spacetime.

Thus, the distinction between "spacetime" and "internal state space" is not a fixed stage, but swaps roles depending on which subalgebra is chosen as the projection axis.

## 4. Covariance of Conservation Laws and Reduction to Minimal Systems

Connecting this to the null condition ∑ xₙ² = 0 from our earliest April models, transposing the components chosen as R₀ to the right-hand side is itself what generates the relative minus sign (Lorentzian metric signature).

Because the original sum-of-squares equation is completely symmetric with respect to which axes are grouped onto the right-hand side, we predict:

**Covariance of conservation laws under axis selection**: Whenever a conserved quantity holds under one choice of projection axis, a corresponding quantity is likewise conserved when another axis (or another 1-, 2-, or 4-dimensional subalgebra) is re-selected as the projection axis.

This stands one level above the Coleman–Mandula theorem, operating at the pre-spacetime stage before a background S-matrix crystallizes.

Practically, this justifies studying the simple minimal system with 5 axes compressed into R'²:

　x² + y² + z² = R'²

Any conservation law for R'² established there is inherited in higher dimensions as a **composite conserved quantity over all five remaining compressed dimensions**:

　R'² = t² + R² + Q₁² + Q₂² + Q₃²

Even in the octonions 𝕆, norm multiplicativity |o₁ o₂| = |o₁||o₂| (the Degen–Graves eight-square identity) holds strictly, ensuring total conservation is never broken while non-associativity drives internal energy redistribution across axes.

## Paper Information

- Title: Projection Axis Selection and Covariance of Conservation Laws ― Two Radii Carried by All Waves, Projection Asymmetry of 1-, 2-, and 4-Dimensional Subalgebras, and Reduction to Minimal Systems ―
- Concept DOI (always latest): https://doi.org/10.5281/zenodo.24273723
- Version DOI (v1.0): https://doi.org/10.5281/zenodo.24273724
- English PDF: https://zenodo.org/records/24273724/files/projection_axis_selection_covariance_conservation_en.pdf
- Japanese PDF: https://zenodo.org/records/24273724/files/%E6%8A%95%E5%BD%B1%E8%BB%B8%E3%81%AE%E9%81%B8%E6%8A%9E%E3%81%A8%E4%BF%9D%E5%AD%98%E5%89%87%E3%81%AE%E5%85%B1%E5%A4%89%E6%80%A7.pdf
- GitHub folder: https://github.com/WurabeSeiji/ai-chat-logs-open/tree/main/%E6%AC%A1%E5%85%83%E3%81%AE%E7%94%9F%E6%88%90%E6%A7%8B%E9%80%A0/%E6%8A%95%E5%BD%B1%E8%BB%B8%E3%81%AE%E9%81%B8%E6%8A%9E%E3%81%A8%E4%BF%9D%E5%AD%98%E5%89%87%E3%81%AE%E5%85%B1%E5%A4%89%E6%80%A7_20261116
- Previous article (EN): https://note.com/kiharanoriaki/n/n368eb19c8f5d

#TheoreticalPhysics #QuantumMechanics #EmergentSpacetime #Octonions #Quaternions #HilbertSpace #HopfFibration #ConservationLaws #Covariance #IndependentResearch #OpenScience #Zenodo
