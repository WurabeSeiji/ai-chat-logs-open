# Bivector Classification of Spin from the Orientation Structure of the 5-Dimensional Hypercube

## — A Natural Derivation of Spin via 3-Cell Orientation of the Dual Polytope —

**Author:** Noriaki Kihara (WF System Co., Ltd.)

**Date:** April 2026

**DOI:** [10.5281/zenodo.19643358](https://doi.org/10.5281/zenodo.19643358)

---

## 0. Purpose of This Paper

The 5-dimensional hypercube (5-cube, penteract) is the dual polytope of the 5-dimensional orthoplex (5-orthoplex), defined by vertices at $(\pm 1, \pm 1, \pm 1, \pm 1, \pm 1)$.

In the author's prior work [6], the geometric classification of spin was derived from the orientation structure of the 5-dimensional orthoplex. In that formulation, spin was characterized through **sign reversals of vertices**. However, in physics, spin is fundamentally angular momentum — a quantity corresponding to a **bivector** (the orientation of a 2-dimensional surface).

This paper reformulates the same spin classification using the dual polytope, the 5-dimensional hypercube, expressing spin as the **normal bivector of 3-cells** (3-dimensional cubic cells). The 40 cubic 3-cells of the 5-dimensional hypercube are naturally indexed by $\binom{5}{2} = 10$ bivectors $e_i \wedge e_j$. This formulation expresses spin directly as face orientation, yielding a physically more transparent derivation.

Under the assumption that the five coordinate axes are classified into "three observable spatial axes" and "two unobservable axes" (Assumption A, Section 1.3), we derive the following:

1. The set of representable spin quantum numbers is limited to six values: $s = 0, \frac{1}{2}, 1, \frac{3}{2}, 2, \frac{5}{2}$, from the orientation degrees of freedom of the 10 bivector types of 3-cells.
2. The distinction between half-integer spin (fermions) and integer spin (bosons) is derived as a geometric difference between 3-cell centers and 4-facet centers.
3. The directionality of forces is determined by the parity of the sign-reversal structure.
4. A generational structure of particles naturally emerges from the permutation symmetry $S_3$ of the three spatial axes.
5. The spin classifications derived from the orthoplex and the hypercube are isomorphic (Duality Theorem).

All results are stated as propositions of discrete geometry and number theory.

---

## 1. Fundamental Structure of the 5-Dimensional Hypercube

### 1.1 Definition and Element Counts

The 5-dimensional hypercube is the regular polytope in $\mathbb{R}^5$ with vertices at $(\epsilon_1, \epsilon_2, \epsilon_3, \epsilon_4, \epsilon_5)$ where $\epsilon_i = \pm 1$.

The number of $k$-dimensional faces is given by $2^{5-k}\binom{5}{k}$.

| $k$ | Element | Count | Computation |
|:---:|:-------:|:-----:|:-----------:|
| 0 | Vertices | 32 | $2^5 \binom{5}{0}$ |
| 1 | Edges | 80 | $2^4 \binom{5}{1}$ |
| 2 | Square faces | 80 | $2^3 \binom{5}{2}$ |
| 3 | Cubic cells (3-cells) | 40 | $2^2 \binom{5}{3}$ |
| 4 | Tesseract facets (4-facets) | 10 | $2^1 \binom{5}{4}$ |

Compared with the element counts $10, 40, 80, 80, 32$ of the dual 5-dimensional orthoplex, the $k$-faces and $(4-k)$-faces are interchanged.

### 1.2 A Numerical Coincidence Specific to 5 Dimensions

The count of edges (80) coincides with the count of square faces (80). This follows from $2^4 \binom{5}{1} = 2^3 \binom{5}{2}$ and is specific to 5 dimensions, dual to the coincidence $\binom{5}{3} = \binom{5}{4}$ in the orthoplex.

### 1.3 Assumption on Axis Classification

The five coordinate axes of the 5-dimensional hypercube are mathematically equivalent. In this paper, we adopt the following assumption.

**Assumption A (Physical classification of axes).** The five coordinate axes are classified as follows.

- **Spatial axes** ($x_1, x_2, x_3$): Three axes corresponding to the 3-dimensional space directly observable by the observer. Since their directions are fixed from the observer's perspective, sign reversal is not permitted. The three axes are symmetric under the permutation group $S_3$.
- **Fourth axis** ($x_4$): An axis corresponding to time. It is orthogonal to 3-dimensional space and not directly observable; sign reversal is permitted.
- **Fifth axis** ($x_5$): An axis corresponding to a large-scale spatial curvature radius. It is orthogonal to 3-dimensional space and not directly observable; sign reversal is permitted.

Under this classification, only the fourth and fifth axes can undergo independent sign reversal, totaling two.

**Remark.** This assumption is naturally derived from the framework of central projection in the author's other work, but the results of this paper follow from this assumption alone. This is the same assumption as in [6].

---

## 2. Bivector Structure of 3-Cells

### 2.1 3-Cells and Normal Bivectors

Each 3-cell (cubic cell) of the 5-dimensional hypercube is defined by choosing three coordinate axes as free axes (each coordinate taking values $\pm 1$) and fixing the remaining two axes at specific values of $\pm 1$.

If the free axes are $\{x_{i_1}, x_{i_2}, x_{i_3}\}$ and the fixed axes are $\{x_j, x_k\}$, the normal direction of this 3-cell is the 2-dimensional plane spanned by the fixed axes, represented by the bivector $e_j \wedge e_k$.

There are $\binom{5}{2} = 10$ ways to choose the two fixed axes, and for each choice, $2^2 = 4$ sign combinations for the fixed values, giving a total of $10 \times 4 = 40$.

### 2.2 Physical Correspondence Between Bivectors and Spin

In physics, spin angular momentum is characterized by the orientation of the plane of rotation. The orientation of a rotation plane is mathematically expressed as a bivector (2-form).

The normal bivector naturally associated with each 3-cell of the 5-dimensional hypercube corresponds precisely to this structure:

$$\text{Spin} \longleftrightarrow \text{Normal bivector of 3-cell} = e_j \wedge e_k$$

This correspondence means that the 5-dimensional hypercube provides a more direct geometric foundation for the derivation of spin than the orthoplex.

### 2.3 Number of Independent Components of $k$-Vectors

The number of independent oriented $k$-dimensional subspaces in 5-dimensional space is $\binom{5}{k}$.

| $k$ | $\binom{5}{k}$ | Geometric meaning | Hypercube correspondence |
|:---:|:--------------:|:-----------------:|:------------------------:|
| 0 | 1 | Scalar (existence) | Whole polytope |
| 1 | 5 | Vector (edge direction) | 4-facet normal |
| 2 | 10 | Bivector (face orientation) | 3-cell normal |
| 3 | 10 | Trivector (volume chirality) | Square face normal |
| 4 | 5 | 4-vector (hypervolume orientation) | Edge normal hyperplane |
| 5 | 1 | Pseudoscalar (total orientation) | Vertex |

Total: $\sum_{k=0}^{5} \binom{5}{k} = 2^5 = 32$

This value of 32 coincides with the number of vertices of the 5-dimensional hypercube.

---

## 3. Geometric Derivation of Spin Values

### 3.1 Classification of 3-Cell Types (Application of Assumption A)

By Assumption A (Section 1.3), the 10 bivector types are classified according to how many non-spatial axes (fourth or fifth) appear in the normal bivector $e_j \wedge e_k$.

| Non-spatial axes $n$ | Normal bivector | Number of types |
|:--------------------:|:---------------:|:---------------:|
| 0 | $e_i \wedge e_{i'}$ ($i, i' \in \{1,2,3\}$) | 3 |
| 1 | $e_i \wedge e_4$ or $e_i \wedge e_5$ ($i \in \{1,2,3\}$) | 6 |
| 2 | $e_4 \wedge e_5$ | 1 |

Here, $n$ denotes the number of sign-reversible axes contained in the normal bivector.

### 3.2 Correspondence Between Bivector Sign Reversal and Spin

Sign reversal of the normal bivector $e_j \wedge e_k$ corresponds to sign reversal of axis $j$ or $k$. Under Assumption A, spatial axes cannot be sign-reversed, so only the number $n$ of non-spatial axes in the normal provides independent sign-reversal degrees of freedom.

When the coupling matrix $C$ introduces sign reversals on $n$ non-spatial axes, the bosonic spin is $s = n$.

| Non-spatial axes $n$ | Structure of $C$ | Bosonic spin | Fermionic spin | 3-cell type |
|:--------------------:|:----------------:|:------------:|:--------------:|:-----------:|
| 0 | $\mathrm{diag}(+,+,+,+,+)$ | 0 | 1/2 | Spatial–spatial normal |
| 1 | $\mathrm{diag}(+,+,+,-,+)$ | 1 | 3/2 | Spatial–non-spatial normal |
| 2 | $\mathrm{diag}(+,+,+,-,-)$ | 2 | 5/2 | Non-spatial–non-spatial normal |

**Proposition 3.1.** The spin quantum numbers representable by the 3-cell orientation structure of the 5-dimensional hypercube are limited to the six values

$$s \in \left\{ 0, \; \frac{1}{2}, \; 1, \; \frac{3}{2}, \; 2, \; \frac{5}{2} \right\}.$$

*Proof.* In the normal bivector $e_j \wedge e_k$, the only sign-reversible axes are the fourth and fifth. The number of non-spatial axes in the normal is $n = 0, 1, 2$. Bosonic spin is $s = n$, giving $s = 0, 1, 2$. Fermionic spin is $s + \frac{1}{2}$, giving $s = \frac{1}{2}, \frac{3}{2}, \frac{5}{2}$. The total is six. $\square$

### 3.3 Number of Spin States

A particle of spin $s$ has $2s + 1$ values of $s_z$.

| $s$ | Number of states $2s+1$ |
|:---:|:----------------------:|
| 0 | 1 |
| 1/2 | 2 |
| 1 | 3 |
| 3/2 | 4 |
| 2 | 5 |
| 5/2 | 6 |
| **Total** | **21** |

$$\sum_{s} (2s+1) = 1 + 2 + 3 + 4 + 5 + 6 = 21 = T(6)$$

21 is the sixth triangular number and also equals $\binom{7}{2}$.

---

## 4. Geometric Origin of Half-Integer Spin

### 4.1 Coordinate Structure of Cell Centers

The center of each cell of the 5-dimensional hypercube has zero coordinates for the free axes and $\pm 1$ for the fixed axes.

**4-facet center:** A 4-facet fixes only one axis. The center of the 4-facet with axis $j$ fixed at $\epsilon_j$ is

$$c_j = (0, \ldots, 0, \epsilon_j, 0, \ldots, 0),$$

which has **one nonzero coordinate**.

**3-cell center:** A 3-cell fixes two axes. The center of the 3-cell with axes $j, k$ fixed at $\epsilon_j, \epsilon_k$ is

$$c_{jk} = (0, \ldots, 0, \epsilon_j, 0, \ldots, 0, \epsilon_k, 0, \ldots, 0),$$

which has **two nonzero coordinates**.

**Proposition 4.1.** In the 5-dimensional hypercube, 4-facet centers have one nonzero coordinate and 3-cell centers have two nonzero coordinates. The set of 4-facet centers coincides with the vertex set of the orthoplex. The 3-cell centers do not lie on the 4-facet center lattice. This non-coincidence constitutes the geometric locus of half-integer spin (fermion) states.

*Proof.* The 4-facet center $c_j = \epsilon_j e_j$ has one nonzero coordinate and coincides with a vertex of the 5-dimensional orthoplex. The 3-cell center $c_{jk} = \epsilon_j e_j + \epsilon_k e_k$ has two nonzero coordinates and does not lie on the lattice of points with only one nonzero coordinate. $\square$

### 4.2 Interpretation via Duality

The contrast between "vertices vs. edge midpoints" in the orthoplex is mapped to "4-facet centers vs. 3-cell centers" in the hypercube.

| | Orthoplex (W6) | Hypercube (this paper) |
|:-:|:--------------:|:----------------------:|
| Integer spin | Vertices (1 nonzero coord.) | 4-facet centers (1 nonzero coord.) |
| Half-integer spin | Edge midpoints (2 nonzero coords.) | 3-cell centers (2 nonzero coords.) |
| Element counts | 10 vertices, 40 edges | 10 4-facets, 40 3-cells |

The distinction between one and two nonzero coordinates is preserved regardless of axis scales.

### 4.3 Geometric Meaning of 720° Rotation

A 4-facet center state (boson) returns to the same point after one full rotation (360°) on the lattice.

A 3-cell center state (fermion) moves to an adjacent center after one rotation and returns to the original center only after two full rotations (720°).

### 4.4 Geometric Derivation of the Pauli Exclusion Principle

**Proposition 4.2.** For half-integer spin states (3-cell centers), the value $s_z = 0$ does not exist. Therefore, the sign of the volume is always determined to be either $+$ or $-$.

*Proof.* The $z$-component of spin $s$ takes values $s_z = -s, -s+1, \ldots, s-1, s$. When $s$ is a half-integer, $s_z = 0$ does not appear in this sequence. $\square$

**Corollary 4.3.** In states where the volume sign is always determined, two states with the same sign and the same quantum numbers cannot coexist at the same spatial point. This corresponds to the Pauli exclusion principle.

---

## 5. Volume Sign and Force Directionality

### 5.1 Volume Signs of All Spin States

We classify each spin state $(s, s_z)$ by the sign of its volume.

**Bosonic states (integer spin):**

| $s$ | $s_z$ | Volume sign |
|:---:|:-----:|:----------:|
| 0 | 0 | $+$ |
| 1 | $+1$ | $+$ |
| 1 | $0$ | $\pm$ |
| 1 | $-1$ | $-$ |
| 2 | $+2$ | $+$ |
| 2 | $+1$ | $+$ |
| 2 | $0$ | $\pm$ |
| 2 | $-1$ | $-$ |
| 2 | $-2$ | $-$ |

**Fermionic states (half-integer spin):**

| $s$ | $s_z$ | Volume sign |
|:---:|:-----:|:----------:|
| 1/2 | $+1/2$ | $+$ |
| 1/2 | $-1/2$ | $-$ |
| 3/2 | $+3/2$ | $+$ |
| 3/2 | $+1/2$ | $+$ |
| 3/2 | $-1/2$ | $-$ |
| 3/2 | $-3/2$ | $-$ |
| 5/2 | $+5/2$ | $+$ |
| 5/2 | $+3/2$ | $+$ |
| 5/2 | $+1/2$ | $+$ |
| 5/2 | $-1/2$ | $-$ |
| 5/2 | $-3/2$ | $-$ |
| 5/2 | $-5/2$ | $-$ |

### 5.2 Parity of Sign Reversals and Force Directionality

**Proposition 5.1.** For bosons whose normal bivector contains an even number of non-spatial axes ($s = 0, 2$), the product of volume signs is $(-1)^n = +1$, and the force is unidirectional (attractive only). For $n$ odd ($s = 1$), the product is $(-1)^n = -1$, and the force is bidirectional (both attractive and repulsive).

| Non-spatial axes | Sign product | Bosonic spin | Force directionality |
|:----------------:|:-----------:|:------------:|:--------------------:|
| 0 (even) | $(+)(+) = +$ | 0 | None (scalar field) |
| 1 (odd) | $(+)(-) = -$ | 1 | $\pm$ (attractive and repulsive) |
| 2 (even) | $(-)(-) = +$ | 2 | $+$ only (attractive only) |

### 5.3 Qualitative Difference Between the Signs of Spin 0 and Spin 2

Spin 0 ($n = 0$, spatial–spatial normal) and spin 2 ($n = 2$, non-spatial–non-spatial normal) both yield a positive volume product, but their origins differ fundamentally.

| | Spin 0 | Spin 2 |
|:-:|:------:|:------:|
| Volume product | $+$ | $+$ |
| Normal bivector | $e_i \wedge e_{i'}$ (spatial axes) | $e_4 \wedge e_5$ (non-spatial axes) |
| Geometric property | Intra-spatial rotation (identity-like) | Time–curvature double reversal |

**Proposition 5.2.** The $(+)(+)$ type and the $(-)(-)$ type both produce a positive volume product, but their algebraic structures differ. The $(+)(+)$ type is the identity transformation, while the $(-)(-)$ type is the self-inverse of sign reversal.

---

## 6. Asymmetry of Volume Signs

### 6.1 Counts of Positive and Negative States

We classify all 21 spin states by their volume signs.

| Classification | Count |
|:--------------:|:-----:|
| Always $+$ | 10 |
| Always $-$ | 9 |
| Both $\pm$ | 2 |
| **Total** | **21** |

**Proposition 6.1.** In the 3-cell orientation structure of the 5-dimensional hypercube, the number of states with always-positive volume (10) exceeds the number with always-negative volume (9) by exactly one. This difference arises because the spin-0 state admits only the $+$ sign.

### 6.2 The Number 10 and the 4-Facet Count

The number of always-positive volume states, **10**, coincides with the number of 4-facets (tesseracts) of the 5-dimensional hypercube, as well as with the number of vertices of the dual orthoplex.

---

## 7. Classification of 3-Cell Types and Particle Hierarchy

### 7.1 Bivector Types and $S_3$ Symmetry

The 10 normal bivector types are classified into four groups under the permutation symmetry $S_3$ of the three spatial axes.

| 3-cell type | Normal bivector | Count | $S_3$ orbit |
|:-----------:|:---------------:|:-----:|:-----------:|
| Spatial–spatial normal | $e_1 \wedge e_2, \; e_2 \wedge e_3, \; e_1 \wedge e_3$ | 3 | 1 orbit |
| Spatial–fourth normal | $e_1 \wedge e_4, \; e_2 \wedge e_4, \; e_3 \wedge e_4$ | 3 | 1 orbit |
| Spatial–fifth normal | $e_1 \wedge e_5, \; e_2 \wedge e_5, \; e_3 \wedge e_5$ | 3 | 1 orbit |
| Non-spatial normal | $e_4 \wedge e_5$ | 1 | 1 orbit |

$$10 = 3 + 3 + 3 + 1$$

### 7.2 Meaning of the Threefold Repetition

Within each group, the three bivector types differ only in the spatial axis index $i = 1, 2, 3$ and are geometrically equivalent under $S_3$.

**Proposition 7.1.** Under the $S_3$ symmetry of the spatial axes, the fermionic 3-cell types form three groups of three.

**Remark.** In the orthoplex [6], the edge types decompose as $3 + 3 + 3 + 1$. In the hypercube, the same decomposition arises for bivector types. Since spin is physically a bivector quantity, this correspondence is more direct.

### 7.3 Correspondence with the Total 3-Cell Count

The total number of 3-cells is 40. This decomposes as follows.

$$40 = 10 \times 2 \times 2$$

| Factor | Geometric meaning |
|:------:|:-----------------:|
| 10 | Bivector types ($\binom{5}{2}$) |
| $\times 2$ | Sign of first fixed axis ($\pm 1$) |
| $\times 2$ | Sign of second fixed axis ($\pm 1$) |

---

## 8. Duality Theorem

### 8.1 Isomorphism of Spin Classifications

**Theorem 8.1 (Duality Theorem).** The spin classification derived from the orientation structure of the 5-dimensional orthoplex and the spin classification derived from the 3-cell orientation structure of the 5-dimensional hypercube are isomorphic under the following correspondence.

| Orthoplex | Hypercube | Common structure |
|:---------:|:---------:|:----------------:|
| Vertices (10) | 4-facets (10) | Integer spin states |
| Edges (40) | 3-cells (40) | Fermion states |
| Edge types ($\binom{5}{2} = 10$) | Bivector types ($\binom{5}{2} = 10$) | Basis of spin classification |
| $S_3$ decomposition $3+3+3+1$ | $S_3$ decomposition $3+3+3+1$ | Generation structure |
| Vertex: 1 nonzero coord. | 4-facet center: 1 nonzero coord. | Boson |
| Edge midpoint: 2 nonzero coords. | 3-cell center: 2 nonzero coords. | Fermion |

*Proof.* The 5-dimensional orthoplex and hypercube are dual polytopes, with $k$-faces in one-to-one correspondence with $(4-k)$-faces. In particular, 0-faces (vertices, 10) of the orthoplex correspond to 4-faces (4-facets, 10) of the hypercube, and 1-faces (edges, 40) correspond to 3-faces (3-cells, 40). Assumption A concerns the coordinate axes of 5-dimensional space and is independent of the choice of polytope. Therefore, all classification structures under $S_3$ symmetry are preserved. $\square$

### 8.2 Advantages of the Hypercube Formulation

The Duality Theorem establishes mathematical equivalence, but the hypercube formulation has the following advantages.

1. **Naturalness of spin**: Spin is fundamentally a bivector (face orientation), and in the hypercube it is directly expressed as the normal bivector of 3-cells. In the orthoplex, it was expressed indirectly through vertex sign reversals.
2. **32 vertices = $2^5$**: The 32 vertices of the hypercube represent all sign combinations of the 5 axes, corresponding directly to the dimension of the Clifford algebra $\mathrm{Cl}(5)$.
3. **4-facet centers = orthoplex vertices**: The coordinates of 4-facet centers are $\epsilon_j e_j$, which coincide exactly with the vertex set of the orthoplex. The spin structures of both polytopes naturally overlap within the same space.

---

## 9. Relation to Prior Work

### 9.1 Direct Precedents

To the author's knowledge, no prior work derives spin as bivectors from the 3-cell orientation structure of the 5-dimensional hypercube. The derivation from the dual orthoplex is due to the author [6].

### 9.2 Related Research

**Topological derivation of spin-statistics:**
- Finkelstein, D. & Rubinstein, J. (1968). "Connection between Spin, Statistics, and Kinks." *J. Math. Phys.* 9, 1762.

**Geometric models of matter:**
- Atiyah, M.F., Manton, N.S. & Schroers, B.J. (2012). "Geometric Models of Matter." *Proc. R. Soc. A* 468, 1252. arXiv: 1108.5151.

**Clifford algebras and the Standard Model:**
- Furey, C. (2016). "Standard Model Physics from an Algebra?" arXiv: 1611.09182.

**5-dimensional unification:**
- Kaluza, Th. (1921). "Zum Unitätsproblem der Physik." *Sitzungsber. Preuss. Akad. Wiss. Berlin* 966-972.

### 9.3 Numerical Correspondences

| Number in this paper | Known structure |
|:--------------------:|:--------------:|
| 21 (total spin states) | $\dim \mathrm{SO}(7) = 21$ |
| 32 (vertices) | $\dim \mathrm{Cl}(5) = 2^5 = 32$ |
| 10 (4-facets = bivector types) | $\dim \mathrm{SO}(5) = 10$ |
| 40 (3-cells) | $10 \times 4$ (bivector types $\times$ signs) |

---

## 10. Open Problems and Future Directions

### 10.1 Complete Construction of the Bivector Algebra

The 10 bivectors $e_i \wedge e_j$ form a basis for the Lie algebra $\mathfrak{so}(5)$. How the spin classification relates to the representation theory of $\mathfrak{so}(5)$ is unknown.

**Problem:** Does the $3 + 3 + 3 + 1$ decomposition correspond to a subalgebra decomposition of $\mathfrak{so}(5)$?

### 10.2 Relation to the Clifford Algebra $\mathrm{Cl}(5)$

The coincidence $32 = 2^5 = \dim \mathrm{Cl}(5)$ suggests that the vertex set of the hypercube is isomorphic to the basis of the Clifford algebra.

**Problem:** Can a Clifford product be naturally introduced on the vertices of the hypercube, and can the spin classification be restated in the language of Clifford algebras?

### 10.3 Unified Formulation of Orthoplex and Hypercube

The Duality Theorem (Theorem 8.1) shows that both polytopes describe the same spin structure from different perspectives.

**Problem:** Does a higher-level structure unifying both polytopes exist (e.g., the $D_5$ lattice obtained as the union of the vertex sets of the orthoplex and hypercube)?

### 10.4 Quantitative Derivation of Mass Ratios

The derivation of mass hierarchies from $S_3$ symmetry breaking remains open (a problem shared with [6]).

### 10.5 Derivation of Gauge Groups

The rigorous construction of the correspondences $S_3 \to \mathrm{SU}(3)$ and bivector structure $\to \mathrm{SU}(2) \times \mathrm{U}(1)$ is unproven.

---

## 11. Conclusion

The following results have been derived from the 3-cell orientation structure of the 5-dimensional hypercube as geometry of bivectors.

1. **Complete classification of spin values**: The six values $s = 0, \frac{1}{2}, 1, \frac{3}{2}, 2, \frac{5}{2}$ are determined by the number of non-spatial axes in the normal bivector of 3-cells.
2. **Fermion–boson distinction**: Derived as the difference in lattice structure between 4-facet centers (1 nonzero coordinate) and 3-cell centers (2 nonzero coordinates).
3. **Force directionality**: The parity of non-spatial axes in the normal bivector determines force directionality.
4. **Origin of particle hierarchy**: The $3 + 3 + 3 + 1$ decomposition of bivector types provides a geometric candidate for the generation structure.
5. **Duality Theorem**: The spin classifications from the orthoplex and hypercube are isomorphic, and the hypercube provides a more natural representation of spin as a bivector.

The physical intuition that spin is the orientation of a surface (bivector) finds its geometric foundation in the 3-cell structure of the 5-dimensional hypercube.

---

## References

[1] D. Finkelstein & J. Rubinstein, "Connection between Spin, Statistics, and Kinks," *J. Math. Phys.* 9, 1762 (1968).

[2] M.F. Atiyah, N.S. Manton & B.J. Schroers, "Geometric Models of Matter," *Proc. R. Soc. A* 468, 1252 (2012). arXiv: 1108.5151.

[3] C. Furey, "Standard Model Physics from an Algebra?" arXiv: 1611.09182 (2016).

[4] Th. Kaluza, "Zum Unitätsproblem der Physik," *Sitzungsber. Preuss. Akad. Wiss. Berlin* 966-972 (1921).

[5] H.S.M. Coxeter, *Regular Polytopes,* Dover (1973).

[6] N. Kihara, "Geometric Classification of Spin Derived from the Orientation Structure of the 5-Dimensional Orthoplex," DOI: [10.5281/zenodo.19630972](https://doi.org/10.5281/zenodo.19630972) (2026).
