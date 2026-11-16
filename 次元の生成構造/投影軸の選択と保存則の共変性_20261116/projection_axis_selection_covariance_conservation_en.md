# Projection Axis Selection and Covariance of Conservation Laws
## ― Two Radii Carried by All Waves, Projection Asymmetry of 1-, 2-, and 4-Dimensional Subalgebras, and Reduction to Minimal Systems ―

**Author:** Noriaki Kihara (WF System Co., Ltd.)  ORCID 0009-0004-6753-4020
**Date:** 2026-11-16    **Version:** v1.0
**Version DOI:** 10.5281/zenodo.24273724    **Concept DOI:** 10.5281/zenodo.24273723
**License:** CC BY 4.0

---

## Abstract

This manuscript proposes and predicts a unified principle—**projection axis selection** and **covariance of conservation laws**—in response to the open question left in the preceding thought-experiment paper [S10]: how spacetime and internal spaces separate and reconstruct from a bounded Hilbert space without assuming a background spacetime.

By connecting the axis-selection symmetry of the equation $\sum x_n^2 = 0$ found in our earliest April central-projection models [S1][S2] and the geometry of the three-layer model $R$–$R_1$–$R_0$ [S3] with the algebraic structures of quaternions and octonions (as well as our analysis of "two $R$'s" in Thought Experiment 7), this paper advances the following four claims and predictions:
1. Rather than a single specific wave (such as a scalar wave) exclusively residing on the phase circle $R_0$, **all waves in the Hilbert space simultaneously reside on both the projection-axis side ($R_0$) and the spatial side ($R_1$)**.
2. If a multiplicatively closed subalgebra is chosen as the projection axis $R_0$, Hurwitz's theorem [X4] restricts the selectable dimensions to **$1, 2, 4$ ($\mathbb{R}, \mathbb{C}, \mathbb{H}$)**, and depending on which dimension is chosen as the projection axis, the outwardly visible space and symmetry ($S^2, \mathbb{C}^3/\mathrm{SU}(3), S^4$) switch accordingly. That is, the distinction between "spacetime" and "internal state space" is itself a matter of projection-axis selection.
3. **Covariance of conservation laws under axis selection:** Whenever a conserved quantity holds under one choice of projection axis, a corresponding quantity is also conserved when another axis (or another 1-, 2-, or 4-dimensional subalgebra) is re-selected as the projection axis.
4. **Justification of reduction to a minimal system:** By this principle, without solving all 8 dimensions simultaneously from the outset, investigating the conservation law of $R'^2$ that holds in the simple minimal system $x^2 + y^2 + z^2 = R'^2$ guarantees that, upon opening the higher dimensions, it is inherited as a **composite conserved quantity over all compressed dimensions**, $R'^2 = t^2 + R^2 + Q_1^2 + Q_2^2 + Q_3^2$.

---

## Introduction: Motivation, Positioning, and Connection to Prior Work

In the preceding paper [S10], we surveyed a line of thought running from the rest property of photons in four-dimensional spacetime to central projection, quaternions, and octonions. However, that paper served primarily as a map of the overall research program (an organization of open problems) rather than directly offering physical "predictions."

The motivation for the present paper arises from a realization during our subsequent examination of the "two $R$'s ($R_0$ and $R_1$)" in quaternions. Initially, we envisioned a division of labor in which a scalar $\cos$ wave (a Higgs-like uniform field) is placed on $R_0$ and three-dimensional space on $R_1$. Closer examination revealed, however, that **there is no reason for only a specific wave to reside on $R_0$; rather, all waves in the Hilbert space simultaneously reside on both $R_0$ and $R_1$**, and furthermore that **the dimension of the projection axis chosen as $R_0$ is itself selectable from $1, 2, 4$, connecting directly to the axis-selection symmetry of $\sum x_n^2 = 0$ from our earliest April models**.

This perspective connects directly to three research currents and fundamental theorems in contemporary theoretical physics:

1. **Quantum Reference Frames (QRF) and Extended Covariance:** Sharing the same problem consciousness as studies showing that physical laws and relational conservation laws are covariant under transformations that switch which quantum subsystem is chosen as the reference frame (clock or coordinate axes) and which is treated as the observed system [X8], this paper generalizes that idea from "particle reference frames" to "the selection of projection axes ($1, 2, 4$ dimensions) of division algebras and light cones."
2. **Quantum Mereology and Factorization Choice:** For arguments that the emergent spacetime and locality depend on how a background-free Hilbert space is factorized into subsystems [X9], this paper provides a concrete algebraic mechanism via partitioning by subalgebras of division algebras.
3. **The Coleman–Mandula Theorem and Division-Algebra Unification:** In standard quantum field theory assuming an $S$-matrix on a background spacetime, a no-go theorem [X7] states that spacetime symmetries and internal symmetries cannot be non-trivially mixed (only a direct product is allowed). In contrast, this paper adopts the stance that **at the pre-spacetime stage without a background spacetime ($\sum x_n^2 = 0$ and all waves over division algebras), the two are on equal footing, and only after a projection axis is chosen do they separate into spacetime ($R_1$) and an internal/phase reference ($R_0$)**, thereby explaining a unified origin of both symmetries from outside the theorem's premises (consistent also with division-algebraic descriptions of elementary particles [X5][X6]).

This paper is not a report of previously published experimental or numerical verification results, but a **prediction paper that clearly formulates the structural principles reached and the inheritance of conservation laws**. The boundaries between established mathematical identities, the working hypotheses/predictions of this paper, and open problems are strictly specified in the "Classification" table at the end.

| Chapter | Theme | Mathematical Structure | Physical Meaning / Prediction |
| :--- | :--- | :--- | :--- |
| Chapter 1 | Simultaneous Membership of All Waves | $o = A(\cos\theta + \mathbf{n}\sin\theta)$ and product $o_1 o_2$ | Not just specific waves, but all waves reside on both $R_0$ and $R_1$ and mutually interfere |
| Chapter 2 | Projection Asymmetry | Hopf fibration $S^1 \hookrightarrow S^3 \to S^2$ | Surjection onto $R_1$, zero projection from $R_0$ (connection to April three-layer model [S3]) |
| Chapter 3 | $1, 2, 4$-Dimensional Axis Selection | $\mathbb{R}, \mathbb{C}, \mathbb{H} \subset \mathbb{O}$ (Hurwitz's theorem) | Choice of axis switches $S^2, \mathbb{C}^3/\mathrm{SU}(3), S^4$ |
| Chapter 4 | Covariance of Conservation Laws | Axis-selection symmetry of $\sum x_n^2 = 0$ | Corresponding conserved quantities are conserved regardless of which axis is chosen (higher symmetry) |
| Chapter 5 | Reduction to Minimal Systems | Composite conservation of $x^2+y^2+z^2=R'^2$ | The $R'^2$ conservation law of the minimal system is inherited as a composite conserved quantity in higher dimensions |

---

## Chapter 1: Two $R$'s of Quaternions and Simultaneous Membership of All Waves

### 1.1 Uniformity of Scalar $\cos$ Waves and the "Center Problem"
In a closed system without externally assumed time evolution (dynamics), consider first a single-component real standing wave (scalar $\cos$ wave) as the simplest element. By the half-angle formula, its squared amplitude (power) is
$$\cos^2\theta = \frac{1}{2} + \frac{1}{2}\cos(2\theta),$$
which carries a uniform offset of $\frac{1}{2}$.

However, if one attempts to place this scalar $\cos$ wave directly onto pure three-dimensional space $(x, y, z)$, a plane wave introduces a preferred wavevector direction (anisotropy), whereas a spherical wave introduces a preferred origin (a center of the universe), breaking spatial equivalence (anonymity) [S6].

### 1.2 Distinguishing the "Two $R$'s" and Simultaneous Membership of All Waves
This problem is avoided by re-reading the geometry of quaternions (the polar form of a quaternion):
$$o = A(\cos\theta + \mathbf{n}\sin\theta) \quad (\mathbf{n} \in S^2,\ \mathbf{n}^2 = -1).$$
One can read this expression as containing **two radii $R$** of different character:

- **$R_0$ (Projection axis; radius of the 1D phase circle $S^1$):** The radius of the circle in the 2D plane on which the real part $\cos\theta$ and the magnitude of the imaginary part $\sin\theta$ reside.
- **$R_1$ (Spatial side; curvature radius of the 2D directional sphere $S^2$ inside the pure-imaginary 3D space $\mathrm{Im}(\mathbb{H}) \cong \mathbb{R}^3$):** The radius on the three-dimensional spatial side carrying direction and extension spanned by the pure-imaginary unit vector $\mathbf{n}$.

Crucially, **there is no reason why only a scalar $\cos$ wave (a Higgs-like uniform field) should reside on $R_0$; rather, every wave in the Hilbert space (all waves) simultaneously resides on both $R_0$ and $R_1$ from the outset**. Any wave represented by a quaternion $o$ necessarily possesses both a real component (on the $R_0$ side) and an imaginary component (on the $R_1$ side) simultaneously. A uniform scalar background is not a separate particle or field, but emerges as a result of all waves sharing the $R_0$ side (or in the limit where the 3D components cancel or average out).

### 1.3 Mutual Interference Between 1D and 3D via Quaternion Multiplication
Because all waves span both $R_0$ and $R_1$, the interference (quaternion multiplication) of two states $o_1, o_2$ directly couples the two sides:
$$o_1 o_2 = A_1 A_2 \Big[ \big(\cos\theta_1\cos\theta_2 - (\mathbf{n}_1\cdot\mathbf{n}_2)\sin\theta_1\sin\theta_2\big) + \big(\mathbf{n}_1\sin\theta_1\cos\theta_2 + \mathbf{n}_2\cos\theta_1\sin\theta_2 + (\mathbf{n}_1\times\mathbf{n}_2)\sin\theta_1\sin\theta_2\big) \Big]$$

- **3D $\to$ 1D (real part):** The inner product $\mathbf{n}_1 \cdot \mathbf{n}_2$ of directions on the 3D spatial side alters the scalar amplitude on the $R_0$ side.
- **1D $\to$ 3D (imaginary part):** The phases on the $R_0$ side ($\cos\theta, \sin\theta$) control the magnitude of the linear combination and twist (cross product $\mathbf{n}_1 \times \mathbf{n}_2$) of the 3D directions.

---

## Chapter 2: Projection Asymmetry in the Hopf Fibration and Correspondence with the April "Three-Layer Model"

### 2.1 Full Projection vs. Zero / Partial Projection
When all waves reside simultaneously on both $R_0$ and $R_1$, the mapping between them is not symmetric; due to the structure of the Hopf fibration $S^1 \hookrightarrow S^3 \to S^2$ [X1], it exhibits a pronounced **projection asymmetry**:

1. **"Full projection (surjection)" onto the $R_1$ side (the 2D directional sphere $S^2$ inside pure-imaginary 3D space $\mathrm{Im}(\mathbb{H}) \cong \mathbb{R}^3$):**
   Every point on the state space $S^3$ of unit quaternions is mapped without omission onto a point of $S^2$ by the Hopf map $\pi : S^3 \to S^2$.
2. **"Zero projection" from the $R_0$ side (1D phase $S^1$ fiber):**
   Rotating the phase along the $S^1$ fiber (changing the global phase) leaves the image on $S^2$ collapsed to a single point without moving. That is, a pure 1D phase change on the $R_0$ side alone cannot determine the direction $\mathbf{n}$ on the $R_1$ side, yielding **zero projection** as spatial displacement.
3. **"Partial projection" as amplitude allocation:**
   Viewed in the quaternion form $\cos\theta + \mathbf{n}\sin\theta$ ($\cos^2\theta + \sin^2\theta = 1$), the phase $\theta$ on the $R_0$ side controls only the weight magnitude $\sin\theta$ appearing on the $R_1$ side without fixing the direction $\mathbf{n}$ (a **partial projection**, becoming zero projection at $\theta = 0, \pi$).

Through this asymmetry, to an observer reading the world from the $R_1$ (spatial) side, the degrees of freedom on the $R_0$ side are invisible as spatial directions and appear instead as a **uniform background (scalar reference, phase, and temporal progression)**.

Note also that in the standard Riemannian submersion of the Hopf fibration $S^1 \hookrightarrow S^3 \to S^2$ (as well as $S^3 \hookrightarrow S^7 \to S^4$ in the next chapter), if the curvature radius of the total space is $R$, the fiber radius is $R_0 = R$, whereas the curvature radius of the base space ($S^2$ or $S^4$) is geometrically fixed to $R_1 = R/2$, i.e., **$R_0 : R_1 = 2 : 1$**. Thus, under the standard metric, the purely geometric ratio is already fixed as a theorem.

### 2.2 Correspondence with and Advance Beyond the April "Three-Layer Model ($R$–$R_1$–$R_0$)"
In the three-layer central-projection model [S3] published by the author in April, regarding concentric subjective spaces $S^n(R_1), S^n(R_0)$ sharing the origin $O$ of a background space $\mathbb{R}^{n+1}$, the following fact was pointed out:
> "A displacement in the $R$ direction connecting corresponding points on different subjective spaces is orthogonal to the tangent space of every subjective space, and therefore cannot be measured by the intrinsic metric of any subjective space."

This connects geometrically directly to the structure in the present paper in which "changes on the $R_0$ side have zero projection from the $R_1$ side and cannot be measured as spatial displacements." However, the present paper advances beyond [S3] in two respects:
1. **Elimination of background space and cross-dimensional decomposition:** Moving beyond similarity transformations between concentric spheres of the same dimension ($R_0 / R_1$), it advances to a framework that decomposes all waves over division algebras—without assuming a background space—into a "projection axis $R_0$ ($1, 2, 4$ dimensions)" and a "spatial side $R_1$."
2. **Acquisition of mutual interference via multiplication:** Rather than merely separating two mutually orthogonal, unmeasurable layers, the layers possess an interference mechanism that rewrites each other's values through the product $o_1 o_2$.

(*Notation remark:* In the April paper [S3], $R_0$ was assigned to large-scale 4D spacetime and $R_1$ to the particle scale; in the present paper, following the notation of Thought Experiment 7, $R_0$ is assigned to the projection axis (phase / internal reference side) and $R_1$ to the outwardly visible spatial side.)

---

## Chapter 3: Generalizing the Projection Axis ($1, 2, 4$ Dimensions) and Resolving the Misreading of "Spacetime vs. Internal Space"

### 3.1 Common Structure Across $\mathbb{R} \to \mathbb{C} \to \mathbb{H} \to \mathbb{O}$
The representation $o = A(\cos\theta + \mathbf{u}\sin\theta)$ in Chapter 1 preserves the exact same form across every level of the real normed division algebras ($\mathbb{R}, \mathbb{C}, \mathbb{H}, \mathbb{O}$) constructed by the Cayley–Dickson process [X3]. This is effectively nothing other than **selecting the single axis of the identity element $1$ (the real axis) from the division algebra as the projection axis $R_0$**, while treating the remaining $(d-1)$-dimensional pure-imaginary space ($0, 1, 3, 7$ dimensions) as the spatial side $R_1$ (see [S5][S7]).

### 3.2 Selecting the Projection Axis from $1, 2, 4$ Dimensions
However, there is no reason to fix the dimension of the projection axis $R_0$ to $1$. If we choose a subalgebra closed under multiplication (interference) as the projection axis $R_0$, Hurwitz's theorem (real normed division algebras exist only in dimensions $1, 2, 4, 8$) [X4] restricts the selectable dimensions within the largest algebra, the octonions $\mathbb{O}$ (8D), to three possibilities: **$1, 2, 4$ ($\mathbb{R}, \mathbb{C}, \mathbb{H}$)**. In each case, the outwardly visible space and symmetry $R_1$ switch as follows:

1. **Projection axis $= 1$D (using $\mathbb{R}$ as the reference):**
   - Quaternions $\mathbb{H}$ ($4 = 1 + 3$): 1D phase/scalar reference $R_0$ and 3D directional space $S^2$.
   - Octonions $\mathbb{O}$ ($8 = 1 + 7$): 1D reference $R_0$ and 7D space $S^6$.
2. **Projection axis $= 2$D (using $\mathbb{C}$ as the reference):**
   - Quaternions $\mathbb{H}$ ($4 = 2 + 2$): For the complex axis $\mathbb{C}$ ($S^1$ fiber), the remainder is the complex projective line $\mathbb{C}P^1 \cong S^2$ (Hopf fibration $S^1 \hookrightarrow S^3 \to S^2$) [X1].
   - Octonions $\mathbb{O}$ ($8 = 2 + 6$): Choosing one complex structure $\mathbb{C}$ as the projection axis decomposes the remaining 6 dimensions into $\mathbb{C}^3$, and as the subgroup of the automorphism group $G_2$ of $\mathbb{O}$ preserving that complex axis, **$\mathrm{SU}(3)$ (the 3-color symmetry)** [X5][X6] emerges. Moreover, quotienting $S^7$ by this $S^1$ phase yields $\mathbb{C}P^3$, corresponding to twistor space.
3. **Projection axis $= 4$D (using $\mathbb{H}$ as the reference):**
   - Octonions $\mathbb{O}$ ($8 = 4 + 4$): Choosing a quaternion axis $\mathbb{H}$ ($S^3$ fiber) as the projection axis makes the remainder the quaternionic projective line $\mathbb{H}P^1 \cong S^4$ (Hopf fibration $S^3 \hookrightarrow S^7 \to S^4$) [X2], i.e., **closed 4D spacetime** [S10].

### 3.3 Resolving the Misreading of "Spacetime vs. Internal State Space"
Conventional intuition tends to read "$R_0$ side = spacetime (the stage), $R_1$ side = internal state space." What the above structure reveals, however, is rather the opposite: **depending on which dimension ($1, 2, 4$) is chosen as the projection axis $R_0$ from the division algebra carrying all waves, which side collapses into the background/internal reference ($R_0$) and which side opens outward to be read as geometric space ($R_1$: $S^2, \mathbb{C}^3, S^4$) swap roles**. This shares the same root as the exchange between observed space and internal space ($1 \Leftrightarrow 3$ duality) found in "Observer Resolution Is a Gauge" [S9].

---

## Chapter 4: Axis-Selection Symmetry of $\sum x_n^2 = 0$ and "Covariance of Conservation Laws Under Axis Selection"

### 4.1 Emergence of Metric Signature by Transposition and Reconvergence with April Axis-Selection Symmetry
While the norms of the division algebras ($\mathbb{R}, \mathbb{C}, \mathbb{H}, \mathbb{O}$) treated in Chapters 1–3 are positive-definite (Euclidean signature: $a^2 + b^2 + \dots = R^2$), the null condition of central projection carries a Lorentzian-type signature. These two are not separate: **starting from an equation in which the sum of squares of all components is zero (or constant), transposing the components chosen as the projection axis $R_0$ to the right-hand side is itself what generates the relative minus sign (Lorentzian metric signature)**.

The conclusion that "how spacetime is read is a problem of projection-axis selection" is a reaffirmation of the principle found in our earliest April central-projection and multi-axis models [S1][S2]:

For example, in an 8-component null condition (complexified light cone),
$$x^2 + y^2 + z^2 + (it)^2 + (iR)^2 + (iQ_1)^2 + (iQ_2)^2 + (iQ_3)^2 = 0,$$
which axis (or set of axes) is kept on the left-hand side (the outwardly visible space) and which is grouped onto the right-hand side (radius $R$ / projection axis) are completely symmetric in the original equation:
- Reading that groups only $R^2$ onto the right-hand side:
  $$x^2 + y^2 + z^2 + (it)^2 + (iQ_1)^2 + (iQ_2)^2 + (iQ_3)^2 = R^2$$
- Reading that groups 5 axes onto the right-hand side as $R'^2 = t^2 + R^2 + Q_1^2 + Q_2^2 + Q_3^2$:
  $$x^2 + y^2 + z^2 = R'^2$$

### 4.2 Central Principle: Covariance of Conservation Laws Under Axis Selection
What this paper predicts as its strongest physical principle is not merely a "freedom to redraw geometry," but the following **covariance of conservation laws**:

> **[Covariance of Conservation Laws Under Axis Selection (Working Hypothesis / Prediction)]**
> Whenever a conserved quantity holds under one choice of projection axis (or subalgebra), because the underlying system is symmetric with respect to axis selection, a corresponding quantity is likewise conserved when another axis (or another 1-, 2-, or 4-dimensional subalgebra) is re-selected as the projection axis.

In the Standard Model, Lorentz symmetry consists of "transformations with the roles of time and space fixed," and gauge symmetry consists of "transformations in an internal space separate from spacetime." In contrast, the present principle—that **conservation laws continue to hold even when the very partition between which axes are chosen as the spacetime/light-cone projection axis and which are placed on the internal/$R$ side is exchanged**—is a higher-level symmetry encompassing both spacetime symmetry and internal symmetry on equal footing. Note that the Coleman–Mandula theorem [X7] is a constraint at the stage after an $S$-matrix on a background spacetime has crystallized; at the pre-spacetime stage ($\sum x_n^2 = 0$ and all waves over division algebras) treated here, this symmetry is not obstructed.

---

## Chapter 5: Justification of Reduction to $x^2 + y^2 + z^2 = R'^2$ and Prediction of Composite Conserved Quantities

### 5.1 Why Studying the Simple Minimal System Is Justified
Covariance of conservation laws under axis selection has a highly practical consequence for theoretical analysis.

Rather than solving all 8 dimensions ($x, y, z, t, R, Q_1, Q_2, Q_3$) simultaneously from the start, it is thoroughly meaningful to first examine the conservation laws holding in the minimal system with 5 axes compressed into $R'^2$ on the right-hand side,
$$x^2 + y^2 + z^2 = R'^2$$
(and the $1 + 3$ structure of quaternions).

This is because "compression" here is not an operation of discarding and erasing degrees of freedom, but an operation of **compositing and bundling the remaining axes into $R'^2$ on the right-hand side** (*of course, a single isolated component such as $x$ in the minimal system is not conserved under higher-dimensional rotations, whereas the **composite quantity over all compressed dimensions** does not break conservation*).

### 5.2 Inheritance of Composite Conserved Quantities and the Role of Octonionic Non-Associativity
If a conservation law for $R'^2$ (or the quaternion norm product $|o_1 o_2| = |o_1||o_2|$) holds in the minimal system $x^2 + y^2 + z^2 = R'^2$, then when the higher dimensions are opened, it is preserved intact as
$$R'^2 = t^2 + R^2 + Q_1^2 + Q_2^2 + Q_3^2,$$
a **composite conserved quantity over all 5 remaining compressed dimensions**.

Likewise, in the octonions $\mathbb{O}$, multiplicativity of the norm $|o_1 o_2| = |o_1||o_2|$ (the Degen–Graves eight-square identity) holds strictly, so conservation of the composite conserved quantity (the total norm product) is never broken even when opened to 8 dimensions. Meanwhile, octonionic **non-associativity** $(o_1 o_2) o_3 \neq o_1 (o_2 o_3)$ plays the role of driving "inter-axis energy redistribution and oscillation" within that conserved total [S10]. Consequently, it is predicted that the allocation rules and invariants obtained at the quaternionic level in the product between 1D phase $R_0$ and 3D direction $R_1$ are inherited, when expanded to the octonionic level, as composite conserved quantities between the projection-axis side $R_0$ (the 4D subalgebra $\mathbb{H}$) and the spatial side $R_1$ ($S^4$).

---

## References

### Self-Citations (Concept DOI, `10.5281/zenodo.<number>`)

[S1] Noriaki Kihara, "Geometric Formulation of 4-Dimensional Spacetime by Central Projection," 2026-04-05. Concept DOI 10.5281/zenodo.19427779.  
[S2] Noriaki Kihara, "Geometric Symmetry of Central Projection: Mathematical Foundation of the Multi-Axis Model," 2026-04-06. Concept DOI 10.5281/zenodo.19434931.  
[S3] Noriaki Kihara, "Three-Layer Structure Model of Central Projection: Nested Subjective Spaces via $R$–$R_1$–$R_0$ and Middleware Geometry," 2026-04-22. Concept DOI 10.5281/zenodo.19691712.  
[S4] Noriaki Kihara, "An Axiomatic Framework for Physical Symmetry Selection by Anonymous Partial Projection," 2026-07-24. Concept DOI 10.5281/zenodo.21521899.  
[S5] Noriaki Kihara, "Mapping from a Set of Complex Numbers to Spacetime (1+3) — Readout and Partial Projection via Quaternions —," 2026-10-18. Concept DOI 10.5281/zenodo.23451494.  
[S6] Noriaki Kihara, "Relational Readout of Quaternions and Inflation as Crystallization — Autonomous Structuring Preserving Anonymity and Partial Projection —," 2026-10-21. Concept DOI 10.5281/zenodo.23523237.  
[S7] Noriaki Kihara, "Did We Truly Understand Complex Numbers? — Where Amplitude Lives, Imaginary-Axis Reading, and Attribution to $\mathbb{C} \to \mathbb{H} \to \mathbb{O}$ —," 2026-10-23. Concept DOI 10.5281/zenodo.23585436.  
[S8] Noriaki Kihara, "Observation Is the Simultaneous Realization of Emission and Absorption — Deriving the Born Rule, Bell Correlations, and the Tsirelson Bound from Viewing-Angle and Magnitude Readout Rules —," 2026-11-01. Concept DOI 10.5281/zenodo.23857884.  
[S9] Noriaki Kihara, "Observer Resolution Is a Gauge — A Thought Experiment on Reverse-Reading Mode Orders and Predictions of a Theory Assembled from Known Symmetries —," 2026-11-10. Concept DOI 10.5281/zenodo.24097611.  
[S10] Noriaki Kihara, "Emergence of Spacetime and Elementary Particles from a Bounded Hilbert Space and Octonions ― From Photon Geometry and Gnomonic Central Projection to Open Problems in Background-Free Quaternionic and Octonionic Dynamics ―," 2026-11-16. Concept DOI 10.5281/zenodo.24267347.

### External References

[X1] H. Hopf, "Über die Abbildungen der dreidimensionalen Sphäre auf die Kugelfläche," *Math. Ann.* 104, 637–665 (1931). DOI 10.1007/BF01457962  
[X2] H. Hopf, "Über die Abbildungen von Sphären auf Sphären niedrigerer Dimension," *Fund. Math.* 25, 427–440 (1935).  
[X3] L. E. Dickson, "On Quaternions and Their Generalization and the History of the Eight Square Theorem," *Ann. Math.* 20, 155–171 (1919). DOI 10.2307/1967865  
[X4] A. Hurwitz, "Über die Composition der quadratischen Formen von beliebig vielen Variabeln," *Nachr. Ges. Wiss. Göttingen*, 309–316 (1898).  
[X5] J. C. Baez, "The Octonions," *Bull. Amer. Math. Soc.* 39, 145–205 (2002). arXiv:math/0105155  
[X6] C. Furey, "Standard model physics from an algebra?" Ph.D. thesis, University of Waterloo (2016). arXiv:1611.09182  
[X7] S. R. Coleman and J. Mandula, "All Possible Symmetries of the S Matrix," *Phys. Rev.* 159, 1251–1256 (1967).  
[X8] F. Giacomini, E. Castro-Ruiz, and Č. Brukner, "Quantum mechanics and the covariance of physical laws in quantum reference frames," *Nature Communications* 10, 494 (2019).  
[X9] S. M. Carroll and A. Singh, "Quantum Mereology: Factorizing Hilbert Space into Subsystems with Quasi-Classical Dynamics," *Phys. Rev. A* 103, 022213 (2021).

---

## Classification

| Category | Content |
| :--- | :--- |
| Identities / Theorems | $\cos^2\theta = \frac{1}{2} + \frac{1}{2}\cos(2\theta)$ / norm multiplicativity of quaternion and octonion products $|o_1 o_2| = |o_1||o_2|$ / Hopf fibrations $S^1 \hookrightarrow S^3 \to S^2$ [X1] and $S^3 \hookrightarrow S^7 \to S^4$ [X2] and geometric ratio $R_0 : R_1 = 2 : 1$ under standard Riemannian submersion / Cayley–Dickson construction [X3] / Hurwitz's theorem (real normed division algebras only in dimensions $1, 2, 4, 8$) [X4] / transposition of $\sum x_n^2 = 0$ to $x^2+y^2+z^2 = R'^2$ |
| Readings / Interpretations | Distinguishing $R_0$ (phase circle) and $R_1$ (spatial side spanning the 2D directional sphere $S^2$ in $\mathrm{Im}(\mathbb{H})$) in quaternions $o = A(\cos\theta + \mathbf{n}\sin\theta)$ / full projection ($R_1$), zero projection ($R_0$ fiber), and partial projection ($\sin\theta$ weight) in the Hopf fibration / correspondence with $R$-direction orthogonality in the April three-layer model [S3] / swapping of spacetime and internal space via choice of $1, 2, 4$-dimensional subalgebra |
| Hypotheses / Predictions | **All waves** in the Hilbert space simultaneously reside on both $R_0$ and $R_1$ / **covariance of conservation laws under axis selection** (corresponding conserved quantities are conserved regardless of which projection axis is chosen) / the $R'^2$ conservation law of the minimal system $x^2+y^2+z^2=R'^2$ is inherited in higher dimensions as a **composite conserved quantity** over all remaining compressed dimensions, with octonionic non-associativity driving internal inter-axis redistribution |
| Open Problems | Dynamical mechanism by which a specific projection axis (and Lorentzian signature) is selected in a concrete observer / dynamical update rule rewriting internal phases rather than $o_1 o_2$ (readout) / concrete numerical derivation of anisotropic scale ratios beyond the standard $2:1$ ratio (Berger-sphere-like ratios), harmonic circulation, and mass hierarchies [S10] |

---

## Version History

- **v1.0 (2026-11-16):** Initial release (Zenodo publication version).
