# Full Spherical Coverage of Central Projection in Discrete Space and the Number-Theoretic Necessity of 5-Dimensional Background Space

## — Resolution of the Inverse Mapping Problem via the 5-Axis 10-Face Relay Model —

**Author:** Noriaki Kihara (WF System Co., Ltd.)

**Date:** April 2026

**DOI:** [10.5281/zenodo.19624957](https://doi.org/10.5281/zenodo.19624957)

**Prerequisite papers:**
- [Geometric Formulation of 4-Dimensional Space via Central Projection](https://doi.org/10.5281/zenodo.19427780)
- [Geometric Symmetries of Central Projection: Mathematical Foundations of Multi-Axis Models](https://doi.org/10.5281/zenodo.19434932)
- [Formulation of Structural Requirements for a Theory of Everything](https://doi.org/10.5281/zenodo.19601592)

---

## 0. Purpose of This Paper

When inversely mapping the entire surface of $S^n(R)$ (the $n$-dimensional hypersphere of curvature radius $R$) onto background coordinates (tangent hyperplanes) using central projection, two essential problems arise.

1. **Equatorial divergence problem**: In regions near the equator as seen from the projection center, the metric of the mapping diverges.
2. **South-hemisphere invisibility problem**: A single tangent hyperplane can only map one hemisphere of the hypersphere.

Starting from fundamental theorems of number theory, this paper derives the following:

- The minimum dimension for lattice distances in discrete space to be defined without gaps is **4** (Lagrange's four-square theorem).
- The background space for central projection onto 4-dimensional discrete space is **5-dimensional** (necessity of homogeneous coordinates).
- The orthogonal axis structure of 5-dimensional space yields a **5-axis 10-face relay model** that resolves both problems above.
- The 5-dimensional cross-polytope (5-dimensional analogue of the regular octahedron) serves as the unit cell for dense packing, and the discrete system closes exactly when the curvature radius satisfies $R^2 \in \mathbb{Z}_{>0}$.
- In the subjective space (4-dimensional), this unit cell is mapped as a **4-dimensional hyperrectangular cell**, and positive/negative volume (orientation) is naturally defined.

All results in this paper are propositions of projective geometry and discrete geometry, and no physical interpretation is introduced.

---

## 1. Discrete Space and the Sum-of-Squares Problem

### 1.1 Problem Setting

In an $n$-dimensional discrete space (integer lattice $\mathbb{Z}^n$), the distance $d$ between lattice points is given by the generalized Pythagorean theorem:

$$d^2 = a_1^2 + a_2^2 + \cdots + a_n^2, \quad a_i \in \mathbb{Z}$$

This $d^2$ takes integer values, though $d$ itself is generally irrational.

When using discrete space as the foundation for geometry, it is an essential requirement for isotropic discretization that **the square of the distance between any two lattice points can take any positive integer value**. If some integer $N$ cannot be expressed as the sum of $n$ squares, then no lattice point exists at that distance, creating **direction-dependent gaps (anisotropy)** in the discrete space.

The question is therefore formulated as follows:

> **Question**: What is the minimum $n$ required to express any positive integer $N$ as a sum of $n$ squares of integers: $N = a_1^2 + a_2^2 + \cdots + a_n^2$?

### 1.2 The Case $n = 1$

$N = a_1^2$ holds only when $N$ itself is a perfect square ($1, 4, 9, 16, \ldots$). The vast majority of positive integers cannot be expressed as a single square.

### 1.3 The Case $n = 2$

By Fermat's two-square theorem, the necessary and sufficient condition for a positive integer $N$ to be expressible as the sum of two squares $N = a_1^2 + a_2^2$ is that every prime factor of $N$ of the form $4k + 3$ appears to an even power.

For example, $N = 3$ is a prime of the form $4 \times 0 + 3$, and no integer solution exists for $3 = a_1^2 + a_2^2$. Similarly, $N = 7, 15, 21, \ldots$ cannot be expressed as the sum of two squares.

Therefore, in a 2-dimensional lattice, no lattice points exist in directions where the square of the distance takes these values, resulting in severe anisotropy in the discretization.

### 1.4 The Case $n = 3$

By Legendre's three-square theorem, a positive integer $N$ **cannot** be expressed as the sum of three squares if and only if $N$ is of the form

$$N = 4^a(8b + 7), \quad a \geq 0, \quad b \geq 0$$

Specifically, $N = 7, 15, 23, 28, 31, 39, 47, 55, 60, \ldots$ cannot be expressed as the sum of three squares. Infinitely many such values exist.

Therefore, isotropic discretization is impossible in principle even in a 3-dimensional lattice. No lattice point exists at distances with $d^2 = 7$ or $d^2 = 15$ in certain directions.

### 1.5 The Case $n = 4$: Lagrange's Four-Square Theorem

**Lagrange's four-square theorem** (1770) states the following:

> **Theorem 1.1 (Lagrange).** Every positive integer $N$ can be expressed as the sum of four squares of integers. That is, there always exist $a_1, a_2, a_3, a_4$ such that
> $$N = a_1^2 + a_2^2 + a_3^2 + a_4^2, \quad a_i \in \mathbb{Z}_{\geq 0}$$

Concrete examples:

- $1 = 1^2 + 0^2 + 0^2 + 0^2$
- $2 = 1^2 + 1^2 + 0^2 + 0^2$
- $3 = 1^2 + 1^2 + 1^2 + 0^2$
- $7 = 2^2 + 1^2 + 1^2 + 1^2$ (impossible with three squares)
- $15 = 3^2 + 2^2 + 1^2 + 1^2$ (impossible with three squares)

**4 is the necessary and sufficient minimum dimension.** As shown in §1.4, $n = 3$ is insufficient, and for $n \geq 5$, setting $a_5 = 0$ reduces to the $n = 4$ case, so 4 is sufficient.

### 1.6 Relation to Waring's Problem

Lagrange's theorem corresponds to the $k = 2$ case of Waring's problem (1770):

> **Waring's problem**: Find the minimum number of terms $g(k)$ required to express any positive integer as a sum of $k$-th powers.

| $k$ (power) | $g(k)$ (minimum terms required) |
|:---:|:---:|
| 2 | **4** (Lagrange's theorem) |
| 3 | 9 |
| 4 | 19 |

The fact that $g(2) = 4$ for squares ($k = 2$) means that **4 is the minimum dimension at which a discrete space can constitute an isotropic lattice without gaps in distance computation**.

### 1.7 The Algebraic Privilege of $n = 4$: Quaternions

The algebraic basis for Lagrange's theorem holding at $n = 4$ lies in the multiplicativity of the quaternion norm.

The norm of a quaternion $q = a + bi + cj + dk$ ($a, b, c, d \in \mathbb{Z}$) is defined as

$$N(q) = a^2 + b^2 + c^2 + d^2$$

This norm is multiplicative:

$$N(q_1 q_2) = N(q_1) \cdot N(q_2)$$

This property (**Euler's four-square identity**) ensures that the product of two four-square sums is again a four-square sum. This is the core of the proof of Lagrange's theorem and is a property intrinsic to the algebraic structure of 4 dimensions (the quaternion algebra $\mathbb{H}$).

Crucially, the quaternions form a **skew field satisfying the associative law**:

- **Real numbers** (1-dimensional): commutative field.
- **Complex numbers** (2-dimensional): commutative field.
- **Quaternions** (4-dimensional): non-commutative but satisfying the associative law; a skew field.
- **Octonions** (8-dimensional): non-commutative and **not satisfying the associative law**.

The 8-dimensional octonions are optimal in packing efficiency (E8 lattice, Viazovska 2016), but the lack of associativity means that computational consistency in projection operations (which involve division) is not guaranteed. **Since central projection is essentially a division operation (normalization of homogeneous coordinates), an algebraic system satisfying the associative law is required**, thereby excluding 8 dimensions from the background space for central projection.

### 1.8 Summary: Why 4 Dimensions

Synthesizing the above arguments, three independent grounds require the dimension of discrete space to be **4** from number-theoretic considerations:

1. **Completeness of distances**: 4 is the minimum dimension at which any positive integer can be expressed as a sum of squares (Lagrange's theorem). For $n \leq 3$, certain distances cannot be represented.
2. **Algebraic closure**: 4 dimensions possess a field structure (quaternions $\mathbb{H}$) with multiplicative norms, being the unique non-commutative dimension with this property.
3. **Projective consistency**: Central projection (a division operation) requires the associative law, and the highest-dimensional field satisfying this is the quaternions (4 dimensions).

These do not mean "4 dimensions happen to be convenient," but rather that **4 dimensions are uniquely determined by fundamental theorems of number theory and algebra**.

---

## 2. Central Projection and the Dimension of Background Space

### 2.1 Definition of Central Projection (Recap)

The operation of mapping a point $P$ on $S^n(R)$ (the $n$-dimensional hypersphere of curvature radius $R$) to the point $P'$ where the line from the center $O$ of the hypersphere through $P$ intersects the $n$-dimensional tangent hyperplane $T$ is called central projection.

In the preceding paper [1], the pullback metric and Einstein tensor were derived from central projection onto $S^4(R)$. The projection target tangent hyperplane (subjective space) is **4-dimensional**, and the background space $S^4(R)$ is a hypersphere in $\mathbb{R}^5$.

### 2.2 Homogeneous Coordinates and Dimensional Relationship

The projective geometry of $n$-dimensional space is described using $(n+1)$-dimensional homogeneous coordinates. Points of the projective space $\mathbb{P}^n$ are defined as equivalence classes of $\mathbb{R}^{n+1} \setminus \{0\}$:

$$[x_1 : x_2 : \cdots : x_{n+1}]$$

Central projection corresponds to dividing the other components by one component (e.g., $x_{n+1}$) of these homogeneous coordinates.

Therefore, **for an $n$-dimensional subjective space, the background space must be $(n+1)$-dimensional**, as required by the structure of projective geometry.

Since §1 derived that the dimension of subjective space is 4, the dimension of background space is **5**.

### 2.3 Sum of Squares in 5-Dimensional Space

The distance between lattice points in the 5-dimensional background space $\mathbb{R}^5$ is given by

$$d^2 = a_1^2 + a_2^2 + a_3^2 + a_4^2 + a_5^2, \quad a_i \in \mathbb{Z}$$

Since Lagrange's theorem showed that any positive integer is representable at $n = 4$, the $n = 5$ case subsumes Lagrange's theorem by setting $a_5 = 0$. Moreover, the additional degree of freedom from five variables increases the number of lattice point combinations realizing the same $d^2$, improving lattice isotropy beyond 4 dimensions.

### 2.4 Why No Dimension Beyond 5 Is Needed

As stated in §1.7, 8-dimensional octonions do not satisfy the associative law. Therefore, central projection from the 5-dimensional hypersphere $S^4(R) \subset \mathbb{R}^5$ to a 4-dimensional tangent hyperplane is completely definable on the algebraic structure of quaternions (4 dimensions), and consistency is guaranteed.

Even if a background space of 6 or more dimensions were introduced, the projection target would still be 4-dimensional, and the additional dimensions would be redundant degrees of freedom not contributing to the projection operation. **The dimension of background space is necessary and sufficient at subjective space dimension + 1 = 5**.

### 2.5 Summary: Why 5 Dimensions

1. **Structure of projective geometry**: Central projection onto $n$ dimensions requires $(n+1)$-dimensional homogeneous coordinates.
2. **Algebraic constraint**: The projection operations in the background space are completed within the field structure of quaternions (associative law).
3. **Elimination of redundancy**: 6 or more dimensions are redundant degrees of freedom that do not contribute to projection, violating the minimality of structural requirements.

---

## 3. Inverse Mapping Problems of Central Projection

### 3.1 Limitations of a Single Tangent Hyperplane

Consider central projection onto the tangent hyperplane $T_N$ at the north pole $N$ of $S^4(R)$. A point $P$ on $S^4(R)$ is mapped to the point $P'$ where the line from the center $O$ (center of the hypersphere) through $P$ intersects $T_N$.

This mapping has two essential limitations:

**Limitation 1 (Equatorial divergence)**: As $P$ approaches the equator (the set of points at geodesic distance $\pi R / 2$ from the north pole), the angle between line $OP$ and $T_N$ approaches 0, so $P'$ diverges to infinity on $T_N$. That is, the metric (Jacobian) of the mapping diverges near the equator.

$$\lim_{\theta \to \pi/2} |P'| = \infty$$

where $\theta$ is the colatitude (angle from the north pole) of $P$.

**Limitation 2 (South-hemisphere invisibility)**: For points $P$ in the south hemisphere beyond the equator, the line $OP$ does not intersect $T_N$ (or intersects on the opposite side, making the mapping non-injective). Therefore, **a single tangent hyperplane can only map one hemisphere of the hypersphere**.

### 3.2 The Essence of the Problem

These limitations are inherent in the geometric structure of central projection and cannot be resolved by any modification of the coordinate system on a single tangent hyperplane.

Symmetry II (equivalence of axes) demonstrated in the preceding paper [2] states that "projection axes can be freely chosen in central projection." That is, constructing a tangent hyperplane at any point other than the north pole yields an equivalent central projection.

By maximally exploiting this property, full spherical coverage—impossible with a single tangent hyperplane—is realized through a **combination of multiple tangent hyperplanes**. This is the 5-axis 10-face relay model of the next section.

---

## 4. The 5-Axis 10-Face Relay Model

### 4.1 Orthogonal Axis Structure of 5-Dimensional Space

In 5-dimensional space $\mathbb{R}^5$, there exist 5 orthogonal axes. Each axis has positive and negative directions, yielding a total of **10 directions** ($\pm e_1, \pm e_2, \pm e_3, \pm e_4, \pm e_5$).

On $S^4(R) \subset \mathbb{R}^5$, each of these 10 directions corresponds to a point on the hypersphere. These 10 points are called **reference points**, and the tangent hyperplanes at each reference point are denoted $T_1^+, T_1^-, T_2^+, T_2^-, \ldots, T_5^+, T_5^-$.

### 4.2 Symmetry of the 10-Face Configuration

The 10 reference points have the following properties on $S^4(R)$:

1. **90-degree spacing**: The geodesic distance between any two reference points is $\pi R / 2$ (90 degrees). This follows immediately from the fact that the orthogonal axes are at right angles.

2. **Equivalence**: The symmetry group of $S^4(R)$ is $O(5)$, and there exists a rotation mapping any reference point to any other. Therefore, all 10 tangent hyperplanes are geometrically equivalent.

3. **Antipodality**: $T_k^+$ and $T_k^-$ correspond to the positive and negative directions of the same axis and are located at antipodal points on the hypersphere.

### 4.3 Construction of Coverage

The central projection from each reference point $p_k$ can regularly map the hemisphere centered at $p_k$ (the region within geodesic distance $\pi R / 2$ from $p_k$). Denote this hemisphere as $H_k$.

**Proposition 4.1.** The 10 hemispheres $H_1^+, H_1^-, H_2^+, H_2^-, \ldots, H_5^+, H_5^-$ completely cover $S^4(R)$. That is,

$$S^4(R) = \bigcup_{k=1}^{5} \left( H_k^+ \cup H_k^- \right)$$

**Proof.** Take any point $P$ on $S^4(R)$. Let the 5-dimensional coordinates of $P$ be $(x_1, x_2, x_3, x_4, x_5)$ ($\sum x_i^2 = R^2$). At least one component satisfies $|x_k| \geq R/\sqrt{5}$ (otherwise $\sum x_i^2 < 5 \cdot R^2/5 = R^2$, a contradiction).

When $|x_k| \geq R/\sqrt{5}$, the geodesic distance from $P$ to the reference point $(\text{sgn}(x_k)) \cdot R \cdot e_k$ is $\arccos(|x_k|/R) \leq \arccos(1/\sqrt{5}) < \pi/2$, so $P \in H_k^{\text{sgn}(x_k)}$. $\square$

### 4.4 The Relay Principle

A point $P$ on $S^4(R)$ generally belongs to multiple hemispheres (the coverage has overlaps). When mapping $P$, one selects the tangent hyperplane where $P$ is mapped most regularly (i.e., with minimum metric distortion).

Specifically, one computes the geodesic distance $\theta_k$ between $P$ and each reference point $p_k$, and selects the tangent hyperplane of the reference point for which $\theta_k$ is smallest. This is equivalent to "choosing the axis corresponding to the coordinate component with the largest absolute value among $P$'s coordinates."

### 4.5 Transformation Between Tangent Hyperplanes

All tangent hyperplanes **share the same center $O$ (origin)** and **have the same curvature radius $R$**. Only the **axis direction** differs. Therefore, coordinate transformations between tangent hyperplanes are **orthogonal transformations** (combinations of rotations and reflections) in 5-dimensional space, and are extremely simple.

Specifically, the transformation between the $k$-th axis and the $l$-th axis is a 90-degree rotation in the $(k, l)$ plane:

$$\begin{pmatrix} x_k' \\ x_l' \end{pmatrix} = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix} \begin{pmatrix} x_k \\ x_l \end{pmatrix}$$

with all other components unchanged.

### 4.6 Resolution of the Equatorial Divergence Problem

A point $P$ near the equator as seen from $T_k$ ($\theta_k \to \pi/2$) satisfies $\theta_l < \pi/2$ for a different axis $l$ ($l \neq k$), and is regularly mapped from the tangent hyperplane $T_l$.

**Theorem 4.2.** In the 5-axis 10-face relay model, for any point $P$ on $S^4(R)$, there exists at least one reference point within geodesic distance $\arccos(1/\sqrt{5}) \approx 63.4°$ ($< 90°$) from $P$. Therefore, metric distortion is bounded, and equatorial divergence does not occur.

**Proof.** This follows immediately from the estimate $|x_k| \geq R/\sqrt{5}$ in the proof of Proposition 4.1. $\square$

### 4.7 Resolution of the South-Hemisphere Problem

Since the tangent hyperplane $T_k^-$ corresponding to the antipodal point $p_k^-$ exists, the region constituting the "south hemisphere" as seen from $T_k^+$ is regularly mapped from $T_k^-$. Furthermore, tangent hyperplanes of other axes also contribute to coverage, so the combination of all 10 faces covers the entire sphere without exception.

### 4.8 Summary

The 5-axis 10-face relay model maximally exploits the symmetry of central projection (equivalence of axes) to achieve full spherical coverage with the following properties:

1. **Completeness**: The entire surface of $S^4(R)$ is mapped.
2. **Regularity**: At any point, metric distortion is bounded (within $\arccos(1/\sqrt{5}) \approx 63.4°$).
3. **Simplicity**: Transformations between tangent hyperplanes are 90-degree orthogonal transformations, sharing the same origin and curvature radius.
4. **Symmetry**: All 10 tangent hyperplanes are geometrically equivalent.

---

## 5. Dense Packing by the 5-Dimensional Cross-Polytope

### 5.1 Definition of the 5-Dimensional Cross-Polytope

The $n$-dimensional cross-polytope (regular hyperoctahedron) $\beta_n$ is the regular polytope in $\mathbb{R}^n$ whose vertices are the $\pm 1$ points on each coordinate axis. The number of vertices is $2n$.

In the 5-dimensional case, $\beta_5$ has the following 10 vertices:

$$(\pm 1, 0, 0, 0, 0), \quad (0, \pm 1, 0, 0, 0), \quad (0, 0, \pm 1, 0, 0), \quad (0, 0, 0, \pm 1, 0), \quad (0, 0, 0, 0, \pm 1)$$

### 5.2 Correspondence Between the Cross-Polytope and the 10-Face Relay

The 10 reference points constructed in §4 are the points on $S^4(R)$:

$$\pm R \cdot e_1, \quad \pm R \cdot e_2, \quad \pm R \cdot e_3, \quad \pm R \cdot e_4, \quad \pm R \cdot e_5$$

These are none other than the vertices of the 5-dimensional cross-polytope $\beta_5$ scaled by $R$.

That is, **the configuration of reference points in the 5-axis 10-face relay model coincides exactly with the vertex structure of the 5-dimensional cross-polytope.** This correspondence is not coincidental but follows necessarily from the orthogonal axis structure of 5-dimensional space.

### 5.3 Partition of Space by the Cross-Polytope

The line segments from the center (origin) of the 5-dimensional cross-polytope to each vertex partition $\mathbb{R}^5$ into $2^5 = 32$ sectors (orthants). On $S^4(R)$, a corresponding partition into $32$ spherical simplices (sections of the hypersphere) is obtained.

Each sector corresponds to one simplex of the cross-polytope, and adjacent sectors share faces. This partition is **isotropic**, with all sectors being congruent.

### 5.4 Condition for Dense Packing

Consider densely packing $S^4(R)$ with a discrete lattice. Let the minimal discrete unit (cell) be the 5-dimensional cross-polytope with edge length (lattice constant) $\ell$.

**Proposition 5.1.** The number $M$ of cells along a great circle (1-dimensional geodesic) of $S^4(R)$ is

$$M = \frac{2\pi R}{\ell}$$

When $M$ is an integer, the geodesic closes exactly at an integer multiple of discrete steps.

In the case of 4-axis symmetry, all four axes share the same $M$, so

$$R = \frac{M \ell}{2\pi}$$

Taking the lattice constant $\ell$ as the unit of length ($\ell = 1$),

$$R = \frac{M}{2\pi}$$

### 5.5 The Integer Condition on $R^2$

For the discrete lattice to **close exactly without phase drift** on $S^4(R)$, the number of steps along geodesics must be consistent in all directions.

The distance $d$ in the 4-dimensional subjective space is given by the Pythagorean theorem:

$$d^2 = n_1^2 + n_2^2 + n_3^2 + n_4^2, \quad n_i \in \mathbb{Z}$$

Requiring that the step count $M$ be an integer when traversing a great circle, and that distances in all directions be expressible as sums of squares of integers by Lagrange's theorem, yields the following condition:

**Theorem 5.2.** A necessary condition for a 4-dimensional discrete lattice to close isotropically and exactly on $S^4(R)$ is that $R^2$ (the square of the curvature radius in units of the lattice constant) is a **positive integer**.

Under this condition, all steps of the discrete lattice are managed as integer values, and the phase drift (cumulative error) upon completing one circuit is zero. **The curvature radius takes not continuous real values but discrete values ($R^2 \in \mathbb{Z}_{>0}$).**

---

## 6. Mapping Shape in the Subjective Space (4 Dimensions)

### 6.1 Mapping from 5-Dimensional Cross-Polytope to 4-Dimensional Hyperrectangular Cell

When a 5-dimensional cross-polytope cell on $S^4(R)$ is centrally projected onto a tangent hyperplane (4-dimensional subjective space), what shape does it take in the projection target?

The projection of a 5-dimensional cross-polytope into 4 dimensions depends on the projection direction; for projection along an axis direction of the cross-polytope (choosing one axis as the projection direction), the shadow is a **4-dimensional regular octahedron** ($\beta_4$, 8 vertices).

On the other hand, in the relay model, the observer in the subjective space seamlessly traverses multiple tangent hyperplanes. In this case, the local discrete lattice is perceived as a lattice with 4-axis symmetry.

**Proposition 6.1.** To an observer in the 4-dimensional subjective space, the minimal cell of the discrete lattice is perceived as a **4-dimensional hyperrectangular cell**.

This is for the following reasons:

1. Since the four axes are symmetric, lattice spacings in each axis direction are equal.
2. Central projection maps straight lines to straight lines (preservation of geodesics), so the orthogonal structure of the background lattice is preserved as orthogonal structure in the subjective space.
3. When a gravitational field (curvature) is present, lattice spacings vary with location, but locally the hyperrectangular shape is maintained.

### 6.2 Maintenance of Density

When 5-dimensional cross-polytope cells cover $S^4(R)$ without gaps, their central projection images—4-dimensional hyperrectangular cells—also cover the subjective space without gaps.

- **Near the projection center**: The lattice is observed as an orderly orthogonal lattice.
- **Regions near the equator**: The lattice is radially stretched, but relay to another tangent hyperplane prevents excessive distortion.
- **Boundaries between adjacent tangent hyperplanes**: They interweave at 45-degree angles, but the observer in the subjective space perceives a continuous lattice.

### 6.3 Lattice Deformation by Gravitational Fields

When a gravitational field exists in the subjective space, the lattice deforms from a perfect hypercube into a hyperrectangular cell (a rectangular solid with unequal edge lengths in each direction). However, two conditions are maintained:

1. **Topological connectivity**: Adjacent cells share faces, and the topology of the lattice is invariant.
2. **Hypervolume preservation**: The 5-dimensional hypervolume of each cell (computed in background coordinates as the hypervolume of the cross-polytope) is invariant. The "distortion" of the lattice observed in the subjective space is merely a reflection of the change in central projection angle.

---

## 7. Volume Orientation and Sign

### 7.1 Orientation in 5 Dimensions

Hypervolume in $\mathbb{R}^5$ has a defined **orientation**. When the order of the vertices of a 5-dimensional cross-polytope is fixed, the hypervolume takes a positive or negative value depending on whether the ordering is right-handed or left-handed.

Mathematically, the sign of the determinant of five edge vectors $v_1, v_2, v_3, v_4, v_5$

$$V = \det(v_1, v_2, v_3, v_4, v_5)$$

determines the orientation, and $|V|$ gives the hypervolume.

### 7.2 Orientation Reversal and Sign of Volume

When centrally projecting a cross-polytope cell on $S^4(R)$, the sign of the Jacobian of the mapping changes depending on the direction of projection (positive or negative of the projection axis).

**Proposition 7.1.** The mapping with reversed orientation of the 5-dimensional cross-polytope (the "inside-out" mapping) is naturally defined as a **hyperrectangular cell with negative hypervolume** in the subjective space (4 dimensions).

This is because the extra dimension of 5-dimensional space provides the distinction between "front and back," enabling a geometric sign assignment to volume that is impossible in 4-dimensional space alone, realized through projection from 5 dimensions.

### 7.3 Coexistence of Positive and Negative Volumes

Cells with different orientations (positive hypervolume) and (negative hypervolume) can coexist at the same coordinates in the subjective space. From the 5-dimensional perspective, these are separate cells belonging to different sectors, but they are mapped to the same coordinate region in the 4-dimensional projection.

In this case, the algebraic sum of positive and negative hypervolumes can be zero. That is, **through the orientation structure of 5 dimensions, the creation and annihilation of volume are geometrically expressible**.

---

## 8. Equivalence of Curvature

### 8.1 5-Dimensional Curvature and 4-Dimensional Closed Space

$S^4(R) \subset \mathbb{R}^5$ is a 4-dimensional hypersurface with curvature radius $R$ in 5-dimensional space. For a 4-dimensional observer living on this hypersurface, the "curvature in the 5th dimension" is not directly observable.

However, such an observer can indirectly perceive the existence of curvature through the following means:

1. **Sum of angles in geodesic triangles**: In flat 4-dimensional space this equals $\pi$ (180 degrees), but on $S^4(R)$ it exceeds $\pi$.
2. **Closure of geodesics**: Traveling straight in any direction, one returns to the starting point after a finite distance ($2\pi R$).
3. **Finiteness of space**: The volume is finite (hypersurface area of $S^4(R) = 8\pi^2 R^4/3$), but no "boundary" exists.

### 8.2 Subjective Equivalence

**Theorem 8.1.** For an observer in the 4-dimensional subjective space, the following two descriptions are observationally indistinguishable (equivalent):

(a) The subjective space is constructed as central projection onto $S^4(R)$ in $\mathbb{R}^5$, with curvature radius $R$ in the 5th-dimensional direction.

(b) The subjective space is a closed 4-dimensional space with finite volume and no boundary, whose curvature radius computed from the deviation of angle sums and closure of geodesics is $R$.

**Sketch of proof.** Central projection maps geodesics to straight lines (preceding paper [1]), and the sum of angles is determined by the curvature of $S^4(R)$. The quantities measurable by the subjective-space observer (angle sums, geodesic periods, volume) are all functions of $R$ alone and contain no additional information about the structure in the 5th dimension. $\square$

This equivalence is the fundamental justification for using central projection as a geometric framework for gravitational fields. The 5-dimensional structure of the background space is completely encoded as curvature in the subjective space, and the observer perceives the background geometry solely through curvature.

---

## 9. Conclusions

### 9.1 Results of This Paper

Starting from fundamental theorems of number theory, this paper has derived the following:

1. **Necessity of 4 dimensions**: By Lagrange's four-square theorem, the minimum dimension for lattice distances in discrete space to be defined isotropically without gaps is 4. This corresponds to the algebraic structure of quaternions (multiplicativity of norms, associative law).

2. **Necessity of 5 dimensions**: By the structure of projective geometry, the background space for central projection onto 4 dimensions is 5-dimensional. Since octonions (8 dimensions) do not satisfy the associative law, they cannot serve as the background space for central projection (a division operation).

3. **5-axis 10-face relay model**: The 10 tangent hyperplanes corresponding to the positive and negative directions of the 5 orthogonal axes in 5-dimensional space (= the 10 vertices of the 5-dimensional cross-polytope) constitute a full spherical coverage of $S^4(R)$. This simultaneously resolves the equatorial divergence problem and the south-hemisphere invisibility problem of central projection.

4. **Dense packing and the integer condition on $R^2$**: A discrete lattice with the 5-dimensional cross-polytope as its unit cell constitutes an exact closed dense packing on $S^4(R)$ when the curvature radius satisfies $R^2 \in \mathbb{Z}_{>0}$.

5. **Mapping to 4-dimensional hyperrectangular cells**: The 5-dimensional cross-polytope cells are mapped as 4-dimensional hyperrectangular cells in the subjective space (4 dimensions). The orientation structure of 5 dimensions naturally defines positive and negative volumes.

6. **Equivalence of curvature**: Curvature in the 5th-dimensional direction and the intrinsic curvature of a closed 4-dimensional space are indistinguishable to an observer in the subjective space.

### 9.2 Scope of This Paper

All results in this paper are stated as propositions of projective geometry, discrete geometry, and number theory, without invoking physical interpretations (gravitational fields, quantum mechanics, etc.). The physical meaning of discrete lattices and the correspondence between hypervolumes of cross-polytopes and physical quantities (energy, etc.) are left as future work in separate papers.

---

## References

[1] N. Kihara, "Geometric Formulation of 4-Dimensional Space via Central Projection," 2026. DOI: [10.5281/zenodo.19427780](https://doi.org/10.5281/zenodo.19427780)

[2] N. Kihara, "Geometric Symmetries of Central Projection: Mathematical Foundations of Multi-Axis Models," 2026. DOI: [10.5281/zenodo.19434932](https://doi.org/10.5281/zenodo.19434932)

[3] N. Kihara, "Formulation of Structural Requirements for a Theory of Everything," 2026. DOI: [10.5281/zenodo.19601592](https://doi.org/10.5281/zenodo.19601592)

[4] J.-L. Lagrange, "Démonstration d'un théorème d'arithmétique," *Nouveaux Mémoires de l'Académie royale des Sciences et Belles-Lettres de Berlin*, 1770.

[5] A. Hurwitz, "Über die Composition der quadratischen Formen von beliebig vielen Variablen," *Nachrichten von der Gesellschaft der Wissenschaften zu Göttingen*, 309–316 (1898).

[6] M. Viazovska, "The sphere packing problem in dimension 8," *Annals of Mathematics* **185**, 991–1015 (2017).

[7] J. H. Conway and N. J. A. Sloane, *Sphere Packings, Lattices and Groups*, 3rd ed., Springer (1999).
