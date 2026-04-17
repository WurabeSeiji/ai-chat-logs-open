# Geometric Classification of Spin Derived from the Orientation Structure of the 5-Dimensional Orthoplex

## — A Unified Derivation of Spin Values, Statistics, and Force Directionality from Integer Geometry —

**Author:** Noriaki Kihara (WF System Co., Ltd.)

**Date:** April 2026

**DOI:** [10.5281/zenodo.19630972](https://doi.org/10.5281/zenodo.19630972)

**Prerequisite Papers:**
- [Full Spherical Coverage of Central Projection in Discrete Space and the Number-Theoretic Necessity of 5-Dimensional Background Space](https://doi.org/10.5281/zenodo.19624957)
- [Formulation of Structural Requirements for a Theory of Everything](https://doi.org/10.5281/zenodo.19601592)
- [Classification and Structural Analysis of the 19 Arbitrary Parameters of the Standard Model](https://doi.org/10.5281/zenodo.19604965)

---

## 0. Purpose of This Paper

In the preceding paper [W5], we demonstrated that the 5-dimensional orthoplex (the 5-dimensional regular cross-polytope) constitutes the unit cell for dense packing of discrete space, and that under central projection it is mapped to a hyperrectangle in the 4-dimensional subjective space. We also established in Proposition 7.1 that the sign of the volume (orientation) is naturally defined through the determinant sign of the 5-dimensional edge vectors.

This paper further analyzes the orientation structure and derives the following:

1. The set of representable spin quantum numbers is limited to six values: $s = 0, \frac{1}{2}, 1, \frac{3}{2}, 2, \frac{5}{2}$, as determined by the orientation degrees of freedom of the 5-dimensional orthoplex.
2. The distinction between half-integer spin (fermions) and integer spin (bosons) is derived as a geometric difference between vertex states and edge-midpoint states.
3. The directionality of forces (attractive only versus attractive and repulsive) is determined by the parity of the sign-reversal structure of the volume.
4. A generational structure of particles naturally emerges from the permutation symmetry $S_3$ of the three spatial axes.
5. The sign asymmetry of spin-0 states (10 positive versus 9 negative) suggests a geometric origin for matter–antimatter asymmetry.

All results in this paper are stated as propositions of discrete geometry and number theory, without invoking physical interpretations. However, correspondences with known particle physics are noted as remarks.

---

## 1. Fundamental Structure of the 5-Dimensional Orthoplex

### 1.1 Definition and Element Counts

The 5-dimensional orthoplex is the regular polytope in $\mathbb{R}^5$ with vertices at $\pm e_i$ ($i = 1, 2, 3, 4, 5$) on each coordinate axis. The edge length is set to 1 (the vertex-to-vertex distance $\sqrt{2}$ satisfies $\ell^2 = 2$, which is an integer).

The number of $k$-dimensional faces is given by $2^{k+1}\binom{5}{k+1}$.

| $k$ | Element | Count | Computation |
|:---:|:-------:|:-----:|:-----------:|
| 0 | Vertices | 10 | $2^1 \binom{5}{1}$ |
| 1 | Edges | 40 | $2^2 \binom{5}{2}$ |
| 2 | Triangular faces | 80 | $2^3 \binom{5}{3}$ |
| 3 | Tetrahedral cells | 80 | $2^4 \binom{5}{4}$ |
| 4 | 5-cell facets | 32 | $2^5 \binom{5}{5}$ |

### 1.2 Integrality of Metrics

With unit edge length, the metric at each dimension is computed as follows.

| Dimension | Metric | Value per element | Rational? |
|:---------:|:------:|:-----------------:|:---------:|
| 0 | Count | 1 | Rational |
| 1 | Edge length | 1 | Rational |
| 2 | Area (equilateral triangle) | $\frac{\sqrt{3}}{4}$ | **Irrational** |
| 3 | Volume (regular tetrahedron) | $\frac{\sqrt{2}}{12}$ | **Irrational** |

Only the 0-dimensional (count) and 1-dimensional (edge length) metrics are rational. Metrics of dimension 2 and above necessarily involve irrational numbers. This indicates the privileged status of edges (1-dimensional elements) in integer geometry.

### 1.3 A Numerical Coincidence Specific to 5 Dimensions

The count of 2-dimensional faces (triangles), 80, coincides with the count of 3-dimensional cells (tetrahedra), 80. This follows from $2^3 \binom{5}{3} = 2^4 \binom{5}{4}$ and is specific to 5 dimensions. It suggests a structure in which face orientation and chirality are in one-to-one correspondence.

---

## 2. Orientation Degrees of Freedom and Independent Oriented Subspaces

### 2.1 Number of Independent Components of $k$-Vectors

The number of independent oriented $k$-dimensional subspaces in 5-dimensional space is $\binom{5}{k}$.

| $k$ | $\binom{5}{k}$ | Geometric meaning | Meaning of the sign |
|:---:|:--------------:|:-----------------:|:-------------------:|
| 0 | 1 | Scalar (existence) | Present / absent |
| 1 | 5 | Vector (edge direction) | $\pm$ direction |
| 2 | 10 | Bivector (face orientation) | $\pm$ rotational sense |
| 3 | 10 | Trivector (volume chirality) | $\pm$ chirality |
| 4 | 5 | 4-vector (hypervolume orientation) | $\pm$ orientation |
| 5 | 1 | Pseudoscalar (total orientation) | $\pm$ total chirality |

Total: $\sum_{k=0}^{5} \binom{5}{k} = 2^5 = 32$

This value of 32 coincides with the number of 4-dimensional facets (5-cells) of the 5-dimensional orthoplex.

### 2.2 Relation to Pascal's Triangle

$$1, \; 5, \; 10, \; 10, \; 5, \; 1$$

This symmetric sequence reflects the Hodge duality between $k$ and $(5-k)$. In particular, the coincidence $\binom{5}{2} = \binom{5}{3} = 10$ means that the orientation of 2-dimensional faces and the chirality of 3-dimensional volumes are in one-to-one correspondence.

---

## 3. Geometric Derivation of Spin Values

### 3.1 Classification of Axes by Physical Role

The five coordinate axes are classified based on symmetry.

- **Spatial axes**: $x_1, x_2, x_3$ (three axes, symmetric under the permutation group $S_3$)
- **Fourth axis**: $x_4$ (asymmetric with respect to spatial axes; first candidate for sign reversal)
- **Fifth axis**: $x_5$ (the additional dimension for projection; second candidate for sign reversal)

The three spatial axes are interchangeable, while the fourth and fifth axes each play an independent role.

### 3.2 Correspondence Between Sign Reversals and Spin

When the coupling matrix $C$ introduces sign reversals on $n$ non-spatial axes, the bosonic spin is $s = n$.

| Sign reversals $n$ | Structure of $C$ | Bosonic spin | Fermionic spin |
|:------------------:|:----------------:|:------------:|:--------------:|
| 0 | $\mathrm{diag}(+,+,+,+,+)$ | 0 | 1/2 |
| 1 | $\mathrm{diag}(+,+,+,-,+)$ | 1 | 3/2 |
| 2 | $\mathrm{diag}(+,+,+,-,-)$ | 2 | 5/2 |

Since there are only two non-spatial axes ($x_4$ and $x_5$), the maximum number of sign reversals is two.

**Proposition 3.1.** The spin quantum numbers representable by the orientation structure of the 5-dimensional orthoplex are limited to the six values

$$s \in \left\{ 0, \; \frac{1}{2}, \; 1, \; \frac{3}{2}, \; 2, \; \frac{5}{2} \right\}.$$

*Proof.* The three spatial axes do not contribute independent sign reversals due to their $S_3$ symmetry. The only axes that can undergo independent sign reversal are the fourth and fifth axes, totaling two. Since bosonic spin equals the number of sign reversals, the bosonic values are $s = 0, 1, 2$. Since fermionic spin is $s + \frac{1}{2}$, the fermionic values are $s = \frac{1}{2}, \frac{3}{2}, \frac{5}{2}$. The total is six. $\square$

**Remark.** The experimentally confirmed spins of fundamental particles are $0, \frac{1}{2}, 1$ only. The values $\frac{3}{2}$ and $2$ are theoretical predictions (gravitino and graviton, respectively). All six values fall within the scope of known theoretical frameworks.

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

### 4.1 Edge Midpoints and Half-Integer Coordinates

The vertices of the 5-dimensional orthoplex lie on the integer lattice (combinations of $\pm 1$ and $0$ along each axis).

The midpoint of the edge connecting two adjacent vertices, for example $e_1 = (1,0,0,0,0)$ and $e_2 = (0,1,0,0,0)$, is

$$m_{12} = \left(\frac{1}{2}, \frac{1}{2}, 0, 0, 0\right),$$

which has **half-integer coordinates**.

**Proposition 4.1.** In the 5-dimensional orthoplex defined on the integer lattice, the edge midpoints lie on a half-integer lattice. These half-integer lattice points constitute the geometric locus of half-integer spin (fermion) states.

### 4.2 Geometric Meaning of 720° Rotation

A vertex state (boson) returns to the same vertex after one full rotation (360°) on the lattice.

An edge-midpoint state (fermion) moves to an adjacent midpoint after one rotation and returns to the original midpoint only after two full rotations (720°).

This is the geometric expression of the property that a spin-$\frac{1}{2}$ particle acquires a phase factor of $(-1)$ under 360° rotation and returns to its original state only after 720° rotation.

### 4.3 Geometric Derivation of the Pauli Exclusion Principle

**Proposition 4.2.** For half-integer spin states (edge midpoints), the value $s_z = 0$ does not exist. Therefore, the sign of the volume is always determined to be either $+$ or $-$.

*Proof.* The $z$-component of spin $s$ takes values $s_z = -s, -s+1, \ldots, s-1, s$. When $s$ is a half-integer, $s_z = 0$ does not appear in this sequence. $\square$

**Corollary 4.3.** In states where the volume sign is always determined, two states with the same sign and the same quantum numbers cannot coexist at the same spatial point. This corresponds to the Pauli exclusion principle.

In contrast, for integer spin (bosons), the state $s_z = 0$ exists and admits both $\pm$ volume signs, permitting multiple occupation of the same state.

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

**Proposition 5.1.** For bosons with an even number of sign reversals ($s = 0, 2$), the product of volume signs is $(-1)^n = +1$, and the force (if mediated) is unidirectional (attractive only). For bosons with an odd number of sign reversals ($s = 1$), the product is $(-1)^n = -1$, and the force is bidirectional (both attractive and repulsive).

| Sign reversals | Sign product | Bosonic spin | Force directionality |
|:--------------:|:-----------:|:------------:|:--------------------:|
| 0 (even) | $(+)(+) = +$ | 0 | None (scalar field) |
| 1 (odd) | $(+)(-) = -$ | 1 | $\pm$ (attractive and repulsive) |
| 2 (even) | $(-)(-) = +$ | 2 | $+$ only (attractive only) |

**Remark.** The fact that the electromagnetic force (mediated by spin-1 particles) exhibits both attraction and repulsion, while gravity (mediated by spin-2 particles) is attractive only, is consistent with experiment.

### 5.3 Qualitative Difference Between the Signs of Spin 0 and Spin 2

Spin 0 (zero sign reversals) and spin 2 (two sign reversals) both yield a positive volume product, but their origins are fundamentally different.

| | Spin 0 | Spin 2 |
|:-:|:------:|:------:|
| Volume product | $+$ | $+$ |
| Origin | $(+)(+)$: no reversal | $(-)(-)$: double reversal |
| Geometric property | Identity transformation | Self-inverse of reversal |

**Proposition 5.2.** The $(+)(+)$ type and the $(-)(-)$ type both produce a positive volume product, but their algebraic structures differ. The $(+)(+)$ type is the identity transformation, while the $(-)(-)$ type is the self-inverse of sign reversal.

**Remark.** In particle physics, the Higgs particle (spin 0) is the source of mass and is itself massive, whereas the graviton (spin 2) is massless. This distinction corresponds to the algebraic difference between the $(+)(+)$ type and the $(-)(-)$ type.

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

**Proposition 6.1.** In the orientation structure of the 5-dimensional orthoplex, the number of states with always-positive volume (10) exceeds the number of states with always-negative volume (9) by exactly one. This difference arises because the spin-0 state admits only the $+$ sign.

### 6.2 The Number 10 and the Vertex Count

The number of always-positive volume states, **10**, coincides with the number of vertices of the 5-dimensional orthoplex.

**Remark.** The observational fact that the universe contains slightly more matter than antimatter (baryon asymmetry) is one of the outstanding problems in cosmology. The result of this paper suggests that the origin of this asymmetry may lie in the sign structure of spin-0 (scalar field) states.

---

## 7. Classification of Edge Types and Particle Hierarchy

### 7.1 Edge Types and $S_3$ Symmetry

Each edge of the 5-dimensional orthoplex connects two vertices on different axes. The number of ways to choose 2 axes from 5 is $\binom{5}{2} = 10$.

Under the permutation symmetry $S_3$ of the three spatial axes, the 10 edge types are classified into four groups.

| Edge type | Axis pair | Count | $S_3$ orbit |
|:---------:|:---------:|:-----:|:-----------:|
| Spatial–spatial | $(x_1, x_2), (x_2, x_3), (x_1, x_3)$ | 3 | 1 orbit |
| Spatial–fourth | $(x_1, x_4), (x_2, x_4), (x_3, x_4)$ | 3 | 1 orbit |
| Spatial–fifth | $(x_1, x_5), (x_2, x_5), (x_3, x_5)$ | 3 | 1 orbit |
| Fourth–fifth | $(x_4, x_5)$ | 1 | 1 orbit |

$$10 = 3 + 3 + 3 + 1$$

### 7.2 Meaning of the Threefold Repetition

The $S_3$ symmetry is the invariance under exchange of the three spatial axes. Within each group, the three edge types differ only in the spatial axis index $i = 1, 2, 3$ and are geometrically equivalent.

**Proposition 7.1.** Under the $S_3$ symmetry of the spatial axes, the fermionic edge types form three groups of three.

**Remark.** The three-generation structure of particle physics (e.g., $e, \mu, \tau$) is one of the outstanding unsolved problems. The result of this paper suggests that the three generations may have a geometric origin in the permutation symmetry of the three spatial axes.

### 7.3 Correspondence with the Total Edge Count

The total number of edges of the 5-dimensional orthoplex is 40. This can be decomposed as follows.

$$40 = 10 \times 2 \times 2$$

| Factor | Geometric meaning |
|:------:|:-----------------:|
| 10 | Edge types ($\binom{5}{2}$) |
| $\times 2$ | Edge direction ($\pm$) |
| $\times 2$ | Hemisphere selection (north/south) |

---

## 8. Relation to Prior Work

### 8.1 Direct Precedents

To the author's knowledge, no prior work derives spin values from the orientation structure of the 5-dimensional orthoplex.

### 8.2 Related Research

The following works are related to individual aspects of this paper.

**Topological derivation of spin-statistics:**
- Finkelstein, D. & Rubinstein, J. (1968). "Connection between Spin, Statistics, and Kinks." *J. Math. Phys.* 9, 1762.

**Geometric models of matter:**
- Atiyah, M.F., Manton, N.S. & Schroers, B.J. (2012). "Geometric Models of Matter." *Proc. R. Soc. A* 468, 1252. arXiv: 1108.5151.

**Clifford algebras and the Standard Model:**
- Furey, C. (2016). "Standard Model Physics from an Algebra?" arXiv: 1611.09182.

**5-dimensional unification:**
- Kaluza, Th. (1921). "Zum Unitätsproblem der Physik." *Sitzungsber. Preuss. Akad. Wiss. Berlin* 966–972.

### 8.3 Numerical Correspondences

We record the correspondences between the numbers arising in this paper and known mathematical structures.

| Number in this paper | Known structure |
|:--------------------:|:--------------:|
| 21 (total spin states) | $\dim \mathrm{SO}(7) = 21$ |
| 32 (4-dimensional facets) | $\dim \mathrm{Cl}(5) = 2^5 = 32$ |
| 10 (vertices) | $\dim \mathrm{SO}(5) = 10$ |

Whether these coincidences are accidental or fundamental is left as a question for future work.

---

## 9. Open Problems and Future Directions

The results of this paper are self-contained as propositions of discrete geometry, but the following problems remain open.

### 9.1 Quantitative Derivation of Mass Ratios

We suggested that the mass hierarchy of three generations may arise from the breaking of the $S_3$ symmetry of spatial axes by the choice of projection axis (Section 7). However, the specific values of mass ratios have not been derived.

**Problem:** Can the mass ratios $m_e : m_\mu : m_\tau$, etc., be derived quantitatively from the angles between the projection axis and each spatial axis?

### 9.2 Derivation of the Gauge Group $\mathrm{SU}(3) \times \mathrm{SU}(2) \times \mathrm{U}(1)$

The permutation group $S_3$ of the spatial axes is related to $\mathrm{SU}(3)$ as a Weyl group, but whether the full Standard Model gauge group can be derived from the symmetry of the 5-dimensional orthoplex remains unproven.

**Problem:** Can the correspondences $S_3 \to \mathrm{SU}(3)$ and edge-type structure $\to \mathrm{SU}(2) \times \mathrm{U}(1)$ be rigorously constructed?

### 9.3 Discrepancy Between 40 Edges and 30 Standard Model States

The number of fundamental fermionic states per generation in the Standard Model is 15 (or 30 including antiparticles), whereas the decomposition of the 5-dimensional orthoplex edges as type $\times$ direction $\times$ hemisphere yields 40. The physical meaning of the difference of 10 is unclear.

**Problem:** Does the discrepancy between 40 and 30 predict the existence of undiscovered states such as right-handed neutrinos, or does the counting correspondence require modification?

### 9.4 Meaning of $21 = \dim \mathrm{SO}(7)$

The coincidence between the total number of spin states (21) and the dimension of $\mathrm{SO}(7)$ is unexplained.

**Problem:** Is the spin state space of the 5-dimensional orthoplex isomorphic to the adjoint representation of $\mathrm{SO}(7)$ or a subspace thereof?

### 9.5 Geometric Origin of Coupling Constants

This paper derived spin values and force directionality but did not address the magnitudes of couplings (coupling constants).

**Problem:** Do values such as the electromagnetic coupling constant $\alpha \approx 1/137$ correspond to geometric quantities (angles, area ratios, etc.) of the 5-dimensional orthoplex?

---

## 10. Conclusion

The following results have been derived purely geometrically from the orientation structure of the 5-dimensional orthoplex.

1. **Complete classification of spin values**: The six values $s = 0, \frac{1}{2}, 1, \frac{3}{2}, 2, \frac{5}{2}$ are necessarily determined by the sign-reversal structure of two non-spatial axes and the half-integer structure of edge midpoints in 5 dimensions.
2. **Derivation of the fermion–boson distinction**: The existence or non-existence of the $s_z = 0$ state provides a geometric account of the differing statistical properties of integer and half-integer spins.
3. **Force directionality**: The parity of the number of sign reversals determines the directionality of the mediated force (unidirectional versus bidirectional).
4. **Origin of particle hierarchy**: The $3 + 3 + 3 + 1$ decomposition of edge types under $S_3$ symmetry provides a geometric candidate for the generation structure.
5. **Matter–antimatter asymmetry**: The difference of 1 between 10 positive-volume states and 9 negative-volume states is attributable to the sign structure of spin 0.

These results indicate the possibility that the axioms of quantum mechanics concerning spin (the spin-statistics theorem, the Pauli exclusion principle) may be re-derived as theorems of 5-dimensional integer geometry.

---

## References

[1] N. Kihara, "Geometric Formulation of 4-Dimensional Space via Central Projection," 2026. DOI: 10.5281/zenodo.19427780

[2] N. Kihara, "Geometric Symmetry of Central Projection: Mathematical Foundations of the Multi-Axis Model," 2026. DOI: 10.5281/zenodo.19434932

[3] N. Kihara, "Formulation of Structural Requirements for a Theory of Everything," 2026. DOI: 10.5281/zenodo.19601592

[4] N. Kihara, "Classification and Structural Analysis of the 19 Arbitrary Parameters of the Standard Model," 2026. DOI: 10.5281/zenodo.19604965

[5] N. Kihara, "Full Spherical Coverage of Central Projection in Discrete Space and the Number-Theoretic Necessity of 5-Dimensional Background Space," 2026. DOI: 10.5281/zenodo.19624957

[6] D. Finkelstein & J. Rubinstein, "Connection between Spin, Statistics, and Kinks," *J. Math. Phys.* 9, 1762 (1968).

[7] M.F. Atiyah, N.S. Manton & B.J. Schroers, "Geometric Models of Matter," *Proc. R. Soc. A* 468, 1252 (2012). arXiv: 1108.5151.

[8] C. Furey, "Standard Model Physics from an Algebra?" arXiv: 1611.09182 (2016).

[9] Th. Kaluza, "Zum Unitätsproblem der Physik," *Sitzungsber. Preuss. Akad. Wiss. Berlin* 966–972 (1921).

[10] H.S.M. Coxeter, *Regular Polytopes*, Dover (1973).
