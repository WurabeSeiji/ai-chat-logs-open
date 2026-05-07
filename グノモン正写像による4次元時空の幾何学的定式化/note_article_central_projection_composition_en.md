# Dimensional Reduction via Central Projection — Algebraizing "Successive Central Projections" Correctly

## Introduction

Central projection (also known as gnomonic projection) is a classical geometric operation that radially projects points of Euclidean space onto a sphere.

In various physical and mathematical frameworks, dimensional reduction from high to low dimensions is sometimes described as "successive central projections along multiple axes." For example: "select axes x₁, x₂, ..., x_k and apply central projection successively to reduce from n to n−k dimensions."

However, this description has a **critical ambiguity**. The fact that the first operation and the second-and-subsequent operations are essentially different kinds of maps is easily overlooked.

This article explains the algebraic formulation that resolves this ambiguity. The full paper is available on Zenodo as a preprint (CC BY 4.0, bilingual md/tex/pdf in both Japanese and English, 6 files total).

- DOI (Concept, auto-redirects to the latest version): https://doi.org/10.5281/zenodo.20060728
- DOI (v1): https://doi.org/10.5281/zenodo.20060729
- GitHub: https://github.com/WurabeSeiji/ai-chat-logs-open/tree/main/グノモン正写像による4次元時空の幾何学的定式化

---

## Part 1: The Decisive Difference Between First and Subsequent Operations

A rigorous observation of dimensional reduction by central projection reveals a **two-stage structure**.

### First Stage: True Central Projection

The map sending a point P of n-dimensional Euclidean space ℝⁿ radially onto the (n−1)-dimensional sphere Sⁿ⁻¹(r₁) of radius r₁ centered at the origin:

　　π(P) = (r₁ / ‖P‖) × P

This is a genuine dimension-reducing map from ℝⁿ to Sⁿ⁻¹(r₁). **This operation occurs only once.**

### Second Stage: Cuts on the Sphere

Once a point lies on the sphere, applying "axis-wise operations" is no longer a central projection.

The correct operation is the **simultaneity-section cut** σ_i, which fixes the axis x_i component of a point X on the sphere Sⁿ⁻¹(r) to a value c:

　　σ_i : Sⁿ⁻¹(r) → Sⁿ⁻²(r')、　 r' = √(r² − c²)

This operation is intrinsic to the sphere and does not directly reference Euclidean background coordinates.

### Essential Difference

- Central projection π: Euclidean background ℝⁿ → sphere Sⁿ⁻¹(r₁)
- Cut σ_i: sphere Sⁿ⁻¹(r) → sphere Sⁿ⁻²(r')

Both are sometimes called "dimension reduction," but **the input space differs**. The first maps Euclidean points onto the sphere; subsequent operations are algebraic operations completed within the sphere.

---

## Part 2: Cuts Are Completely Commutative

### Main Theorem: Commutativity

The first main result of the paper is that **cuts along distinct axes yield the same result regardless of order**.

Mathematically, for two distinct axes x_i, x_j on the sphere Sⁿ⁻¹(r), the cuts σ_i, σ_j satisfy:

　　σ_j ∘ σ_i = σ_i ∘ σ_j

### Why Commutativity Holds

Intuitively: a cut along axis x_i "removes the x_i component" but leaves other axis values untouched. Therefore the value of the x_j component is unchanged before and after the cut σ_i.

Order A: cut x_i then x_j
- Radius: r → √(r² − c_i²) → √(r² − c_i² − c_j²)

Order B: cut x_j then x_i
- Radius: r → √(r² − c_j²) → √(r² − c_j² − c_i²) = √(r² − c_i² − c_j²)

The final radii are the same. The remaining coordinate components are the same (the residue after both x_i and x_j are removed). Hence the maps coincide.

### Composition of k Axes

When cutting along k axes, any permutation yields the same map (this follows as a corollary from adjacent transpositions).

This allows defining the **composite cut** σ_S for an axis set S = {i₁, ..., i_k} as a unique operation independent of order.

---

## Part 3: The Pythagorean Closed Form for the Composite Curvature Radius

### Main Theorem: Closed Form

After applying the composite cut along the axis set S, the final sphere has radius given by:

　　r_final² = r₁² − Σ_{i ∈ S} (x_i*)²

where x_i* is **the value of the axis x_i component on the sphere immediately after the first central projection**.

### Connection to the Pythagorean Theorem

A point P = (x₁, ..., xₙ) on the sphere Sⁿ⁻¹(r₁) satisfies, by definition:

　　r₁² = x₁² + x₂² + ... + xₙ²

This is precisely the n-dimensional Pythagorean theorem.

When cutting along the axis set S, the components corresponding to S are "subtracted" to form the final radius:

　　r_final² = Σ_{i ∉ S} x_i² = r₁² − Σ_{i ∈ S} (x_i*)²

In other words, **the coordinate components of a point on the sphere contribute independently as an orthogonal decomposition by the axis set**.

### Remaining Coordinates Are Invariant

The coordinate components of axes not in the cut set remain **unchanged** under the composite cut σ_S. This follows directly from the fact that "cutting along axis x_i does not affect other axis values."

### Invertibility

The composite cut σ_S is **invertible**. Given the cut-axis values {x_i*} (i ∈ S) and the post-cut image point σ_S(P), the original pre-cut sphere point P can be uniquely recovered.

---

## Part 4: The Algebraic Structure as an Abelian Semigroup

The set of axis-cut operations {σ_S | S ⊆ {1, ..., n}}, with composition ∘ as its operation, satisfies:

- **Commutativity**: σ_S ∘ σ_T = σ_T ∘ σ_S (when S and T are disjoint)
- **Associativity**: (σ_S ∘ σ_T) ∘ σ_U = σ_S ∘ (σ_T ∘ σ_U)
- **Identity**: σ_∅ = id (the identity map)

This is the structure of an **Abelian semigroup**, isomorphic to subset union (when disjoint).

There are no inverse elements (an inverse operation that increases dimension is not defined), so this is a semigroup, not a group.

---

## Part 5: Numerical Verification (7-dimensional → 5-dimensional Example)

We apply central projection with r₁ = 10 to the 7-dimensional Euclidean point P = (1, 2, 3, 4, 10, 5, 3), then cut along axes x₅ and x₆.

After central projection:
　　P' = (10/√164) × P ≈ (0.781, 1.562, 2.343, 3.123, 7.809, 3.904, 2.343)

Axis x₅ component: x₅* ≈ 3.904
Axis x₆ component: x₆* ≈ 2.343

**Order A** (x₅ → x₆):
- After x₅ cut: √(100 − 15.244) ≈ 9.206
- After x₆ cut: √(84.756 − 5.488) ≈ √79.268 ≈ 8.903

**Order B** (x₆ → x₅):
- After x₆ cut: √(100 − 5.488) ≈ 9.722
- After x₅ cut: √(94.512 − 15.244) ≈ √79.268 ≈ 8.903

**Verification by closed form**:
　　r_final² = 100 − 15.244 − 5.488 = 79.268
　　r_final = √79.268 ≈ 8.903

Order A, Order B, and the closed form all agree.

---

## Part 6: Why This Result Is Powerful

The core results of this paper, viewed by experts, are nothing more than **the Pythagorean decomposition on the sphere and the commutativity of cuts** — no new theorems are proved.

Yet the paper is extremely powerful, for the following reasons.

### The Power of Making the Implicit Explicit

The ambiguous description "successive central projections along multiple axes" has appeared in various contexts. This paper organizes that description as the algebraic operation **one central projection plus k cuts on the sphere**. The mathematical content is the same, but once the orthography is established, the logic built upon it becomes transparent.

### Commutativity Breaks the Curse of Order

The primary weakness of successive operations is order dependence. Discussions of "should we cut x_i first or x_j first?" are completely eliminated. Physical descriptions are guaranteed to depend only on the **set of axes**, and the absence of hidden path dependence is established algebraically.

### The Closed Form r² = r₁² − Σ(x_i*)² Provides a Single Computational Rule

The fact that k operations can be written in a single equation means that any n → d dimensional reduction is **computable in one line**. Reductions like 24 → 4, 16 → 4, 9 → 4, 6 → 4 are all processed as different inputs to the same equation.

### Algebraic Leverage from the Abelian Semigroup Structure

The existence of an Abelian semigroup means that algebraic tools such as representation theory, tensor product structures, homological algebra, and categorical languages become applicable. Compatibility with related fields like Wilson lattice gauge theory follows naturally.

### The Power of an Indisputable Foundation

This rests solely on the Pythagorean decomposition of the sphere and the commutativity of cuts, so it is **self-evident enough that it cannot be proved**. Precisely because of this, applied papers that depend on it preserve their foundational integrity.

---

## Part 7: Historical Parallels

Results that are "self-evident yet powerful" appear repeatedly at the turning points of mathematics and physics.

- The inner product of Hilbert space: supports the entire structure of quantum mechanics
- The alternation of differential forms dx ∧ dy = −dy ∧ dx: gateway to de Rham cohomology
- The Lagrangian form of Maxwell's equations: starting point of gauge theory
- Yoneda's lemma: pervades the entire structure of category theory

In all cases, the proofs are nearly trivial, but the formulations unleash modern structures. The axis-cut semigroup of this paper belongs in this lineage.

---

## Part 8: What This Paper Does Not Address

This paper is purely a paper on geometric algebraic structure. The following are explicitly out of scope.

- Concrete interpretation of axes x_i (specific meaning in spacetime, observers, gauges, relativity, etc.)
- Concrete meaning of axis selection S
- Physical necessity of specific dimensions (n, d)
- Generalization to non-Euclidean backgrounds (pseudo-Riemannian metrics, curved manifolds)

These are problems addressed in separate papers. The algebraic results here hold independently of those application contexts.

---

## Conclusion

The results of this paper, on the surface, amount to little more than spherical Pythagoras and the commutativity of cuts. But by organizing them as the **correct algebraic language**:

- Reduction from high to low dimensions can be uniquely described as a closed operation on the combinatorial object "set of axes"
- Confusion over order dependence is completely resolved
- Any application area (physics, statistics, machine learning, control theory, etc.) can use it as a common foundational language

There are moments when establishing the correct symbolic system is more important than proving new theorems. This paper takes that stance.

---

## References (from the paper)

- Snyder, J.P. (1987). *Map Projections—A Working Manual*. U.S. Geological Survey Professional Paper 1395.
- Howie, J.M. (1995). *Fundamentals of Semigroup Theory*. London Mathematical Society Monographs, Oxford University Press.
- Kihara, N. (2026). *A Geometric Formulation of 4-Dimensional Space via Central Projection*. Zenodo. DOI: 10.5281/zenodo.19427780.

---

## Related Links

- Full paper (PDF, LaTeX, Markdown; bilingual JA/EN): https://zenodo.org/records/20060729
- Concept DOI (auto-redirects to the latest version): https://doi.org/10.5281/zenodo.20060728
- Author ORCID: https://orcid.org/0009-0004-6753-4020
- License: CC BY 4.0

---

#CentralProjection #GnomonicProjection #ProjectiveGeometry #Geometry #Mathematics #Algebra #Semigroup #AbelianSemigroup #PythagoreanTheorem #SphericalGeometry #DimensionalReduction #HighDimensionalGeometry #Commutativity #MathematicalPhysics #TheoreticalPhysics #LatticeGaugeTheory #DifferentialGeometry #FourDimensionalGeometry #Preprint #Zenodo
