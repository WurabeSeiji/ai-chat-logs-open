# Geometric Formulation of Spacetime via Central Projection: Overview of 10 Papers

**Author:** Noriaki Kihara  
**Affiliation:** WF System Co., Ltd. (Osaka University, Faculty of Engineering Science, graduated)  
**ORCID:** 0009-0004-6753-4020  
**Date:** April 2026

---

## Paper List

| # | Title | DOI |
|---|-------|-----|
| [1] | Geometric Formulation of 4-Dimensional Space via Central Projection | [10.5281/zenodo.19427780](https://doi.org/10.5281/zenodo.19427780) |
| [2] | Geometric Symmetries of Central Projection: Mathematical Foundations of the Multi-Axis Model | [10.5281/zenodo.19434932](https://doi.org/10.5281/zenodo.19434932) |
| [3] | Relativity of Observation in Multiple Subjective Spaces: Geometric Consequences of the Symmetries of Central Projection | [10.5281/zenodo.19435162](https://doi.org/10.5281/zenodo.19435162) |
| [4] | On the Limits of the Curvature Radius in Central Projection | [10.5281/zenodo.19526549](https://doi.org/10.5281/zenodo.19526549) |
| [5] | On Dimension Addition to Background Spaces and Subjective Spaces | [10.5281/zenodo.19526913](https://doi.org/10.5281/zenodo.19526913) |
| [6] | Subjective Spaces from Zero to Four Dimensions | [10.5281/zenodo.19533292](https://doi.org/10.5281/zenodo.19533292) |
| [7] | On the Continuity of Geodesics in Subjective Spaces | [10.5281/zenodo.19533299](https://doi.org/10.5281/zenodo.19533299) |
| [8] | Diameter of the Circumscribed Hypersphere of a Unit Four-Dimensional Hyperrectangle | [10.5281/zenodo.19533313](https://doi.org/10.5281/zenodo.19533313) |
| [9] | Schwarzschild–de Sitter Exact Solution in the Central Projection Framework | [10.5281/zenodo.19538098](https://doi.org/10.5281/zenodo.19538098) |
| [10] | Dimensional Interpretation of Geodesic Structure at the $R \to 0$ Limit | [10.5281/zenodo.19538106](https://doi.org/10.5281/zenodo.19538106) |

---

## Purpose of This Research

From purely geometric properties of central projection from an $n$-dimensional tangent hyperplane to a hypersphere $S^n(R)$, we construct a unified framework in which multiple independent axes that generate acceleration can coexist without contradiction. We refrain from physical identifications and confine ourselves to geometric propositions and their proofs.

---

## Proven Geometric Propositions

### Paper [1]: Basic Structure of Central Projection

**1. Isomorphism with de Sitter Spacetime**  
The induced metric of the Lorentzian central projection coincides numerically with the metric in the Beltrami coordinate system of de Sitter spacetime, establishing complete equivalence as intrinsic geometry.  
(§1.1, Appendix C, Proposition C.1)

**2. Regularity of Forward and Inverse Projection**  
Since $l = \sqrt{R^2 + r^2} \geq R > 0$, the central projection $\Phi$ is $C^\infty$-regular on its entire domain. The inverse $x^\mu = RY^\mu/Y^{n+1}$ is regular on the upper hemisphere $H^+$. The mapping is bijective.  
(Appendix B, Proposition B.1)

**3. Induced Metric and Einstein Tensor**  
The pullback metric $g_{\mu\nu} = (R^2/l^2)(\delta_{\mu\nu} - x_\mu x_\nu/l^2)$ is rigorously derived, and the $n$-dimensional Einstein tensor satisfies $G_{\mu\nu} + \Lambda_n g_{\mu\nu} = 0$ with $\Lambda_n = (n{-}1)(n{-}2)/(2R^2)$.  
(§3, §6, Proposition 6.1)

**4. Flat Limit**  
As $R \to \infty$: $g_{\mu\nu} \to \delta_{\mu\nu}$, $\Lambda \to 0$. The central projection degenerates to the identity map, recovering Minkowski spacetime.  
(Appendix A)

---

### Paper [2]: Four Geometric Symmetries

**5. Symmetry I (Discrete Stability)**  
For any coordinate value — zero, positive integer, or negative integer — $l_A > 0$ holds and the central projection is globally regular without breakdown.  
(§3, Theorem 3.1)

**6. Symmetry II (Axis Equivalence)**  
Regardless of which axis in $\mathbb{R}^{n+1}$ is chosen as the projection center, the metric, curvature tensor, and Einstein tensor have identical mathematical structure.  
(§4, Theorem 4.1)

**7. Symmetry III (Geodesic Deviation and Indistinguishability)**  
The geodesic deviation $|\xi|^2 = g_{\mu\nu}\xi^\mu\xi^\nu$ is a scalar invariant of dimension $[L^2]$ and the sole means by which an internal observer recognizes curvature. However, neither the sign of coordinates nor the direction of curvature can be distinguished. The direction of acceleration is defined only in background coordinates.  
(§5, Theorem 5.2)

**8. Symmetry IV (Transformability of Subjective Coordinate Systems)**  
Subjective coordinate systems centered on different axes are mutually transformable via $T_{A \to B} = \Phi_B^{-1} \circ \Phi_A$, a $C^\infty$ diffeomorphism. The only parameters are the curvature radius $R$ and coordinate values. Background coordinates are preserved.  
(§6, Theorem 6.1, Proposition 6.5)

**9. Curvature Radius as Scale Factor**  
The Jacobian of coordinate transformations between subjective systems is $\partial x'^A/\partial x^B = -R^2/(x^B)^2$, with a dimensionless ($[L]/[L]$) first-order scale factor. The curvature radius $R$ serves as the sole transformation parameter.  
(§6.4, Proposition 6.2)

---

### Paper [3]: Multiple Subjective Spaces and Relativity of Observation

**10. Simultaneous Existence of Multiple Subjective Spaces**  
From the ambient space $\mathbb{R}^{n+1}$, $n+1$ subjective spaces $\{H_A^+\}$ can be simultaneously constructed by the choice of the projection center axis. Each subjective space possesses structurally identical geometric properties by Symmetry II.  
(§2, Proposition 2.1)

**11. Limits of Information Obtained by Internal Observers**  
An internal observer obtains only quantities constructed from the metric of their own subjective space and cannot access the value of the projection center axis, absolute position in the ambient space, or the existence of other subjective spaces.  
(§3, Proposition 3.1)

**12. Exchange of Axis Roles and Relativity of Observation**  
What is the "projection center axis" (an external axis not directly accessible) for observer $A$ appears as one of the subjective coordinates for observer $B$, and vice versa. The roles of axes that are equivalent in the ambient space are exchanged by the choice of the observer's standpoint. There exists no geometric basis for determining which observer's standpoint is "correct."  
(§4, Proposition 4.1; §5, Theorem 5.1)

---

### Paper [4]: Limits of the Curvature Radius

**13. Regularity of the $R \to \infty$ Limit**  
As $R \to \infty$, the induced metric asymptotes to flat Euclidean space and $\Lambda = 3/R^2 \to 0$. All terms of the Einstein equation exhibit mathematically well-defined asymptotics.  
(§2.1)

**14. Regularity and Structural Tension at $R \to 0$**  
As $R \to 0$, $\Lambda = 3/R^2$ diverges, but the induced metric, curvature tensor, and all terms of the Einstein equation remain well-defined for any $R > 0$. The upper bound of geodesic deviation $|\xi|^2_{\max}(R)$ approaches zero while the local curvature of the subjective space diverges at order $1/R^2$, giving rise to structural tension.  
(§3.1, §3.2)

**15. Existence of Nested Structure**  
For any point $p$ in any subjective space $\Pi_R$, a new subjective space $\Pi_{R'}$ with $p$ as the elevation point exists for any $R' > 0$. This construction is directly realizable from the existing degrees of freedom in Papers [1] and [2].  
(§4.2, Proposition 4.1)

---

### Paper [5]: Dimension Addition

**16. Two Categories of Dimension Addition**  
Dimension addition is classified into two types: Type I (addition under the governance of an existing curvature radius $R$) and Type II (addition involving the introduction of a new independent curvature radius $R'$).  
(§2.2)

**17. Dimension Addition to Finite-Curvature Subjective Spaces**  
Type I addition is freely executable due to axis equivalence. For Type II, unless the axis corresponding to the new curvature radius and the governed additional dimensions are added simultaneously, the resulting subjective space degenerates to zero dimensions.  
(§4.2, Observation 4.1; §4.3, Observations 4.2–4.3)

---

### Paper [6]: Subjective Spaces from Zero to Four Dimensions

**18. Characteristics of Subjective Spaces at Each Dimension**  
0D: No internal structure, geodesic deviation, or curvature recognition. 1D: No concept of direction; geodesic deviation not recognizable. 2D: Area, angle, and geodesic deviation become definable for the first time, but the choice of time axis is undetermined. 3D: Volume becomes definable; interpreting one axis as time yields two spatial axes. 4D: Hypervolume becomes definable, but the number of time axes is not uniquely determined.  
(§2–§6)

---

### Paper [7]: Continuity of Geodesics

**19. Independence of Parent-Child Geodesics**  
A geodesic $\gamma$ on the parent subjective space $\Pi_R$ and a geodesic $\gamma'$ on the child subjective space $\Pi_{R'}$ (elevated from a point on $\gamma$) are independently defined geometric objects within their respective subjective spaces, while sharing the same ambient space.  
(§4.3, direct consequence of Theorem 6.1)

---

### Paper [9]: Schwarzschild–de Sitter Exact Solution

**20. Unique Construction of the Exact Solution**  
The spherically symmetric static exact solution of the vacuum Einstein equation $G_{\mu\nu} + \Lambda g_{\mu\nu} = 0$ ($\Lambda = 3/R^2$) from Paper [1] is uniquely determined as the Schwarzschild–de Sitter metric $f(r) = 1 - 2M/r - r^2/R^2$ by Birkhoff's theorem (with $\Lambda$ extension). No additional assumptions are required.  
(§3, Eq. 3.4)

**21. Disappearance of Horizon Structure at $R \to 0$**  
As $R \to 0$, the cosmological horizon $r_{\mathrm{C}} \to 0$, merging with the black hole horizon and eliminating the region where an internal observer can exist. This is the exact-solution-level counterpart of the structural tension in Paper [4] §3.2.  
(§4.2–§4.3)

**22. Recovery of Horizon Structure in the Child Subjective Space**  
By Proposition 4.1 of Paper [4], a child subjective space $\Pi_{R'}$ can be constructed from the $R \to 0$ limit region, with $R'$ independent of $R$. In the child's exact solution $f'(r') = 1 - 2M'/r' - r'^2/R'^2$, the inter-horizon observer region recovers when $R'$ is finite.  
(§5.2–§5.3)

---

### Paper [10]: Dimensional Interpretation of the $R \to 0$ Limit

**23. Correspondence between $R \to 0$ and Zero-Dimensional Subjective Space**  
The structure of a subjective space at the $R \to 0$ limit geometrically corresponds to the description of a zero-dimensional subjective space in Paper [6]: disappearance of geodesic deviation, disappearance of internal structure, and disappearance of curvature recognition.  
(§3, Observation 3.1)

**24. Description as Dimensional Transition**  
Tracking a geodesic of the parent subjective space toward $R \to 0$ is described in two stages. Stage 1 (Dimensional collapse): the parent's dimensional structure is lost, approaching zero dimensions. Stage 2 (Dimensional reconstruction): a child subjective space is elevated via Proposition 4.1, and dimensional structure unfolds. The "terminus" of the geodesic is not the annihilation of structure but a transition point of dimensional structure.  
(§5, Observation 5.1)

---

### Paper [8]: Singularity of Four Dimensions

**25. Integer Property of the Circumscribed Hypersphere Diameter**  
The diameter $2\sqrt{n}$ of the circumscribed hypersphere of a unit $n$-dimensional hypercube is an integer if and only if $n$ is a perfect square. Excluding the trivial $n = 1$, $n = 4$ is the smallest non-trivial dimension satisfying the integer property.  
(§3)

**26. Settling at Four Dimensions**  
The dimension that simultaneously satisfies symmetry preservation, integer property, circumscribed sphere diameter minimization, and hypervolume minimization — the "settling dimension" — is $n = 4$. Starting from zero dimensions and sequentially releasing degrees of freedom, dimensions 1–3 cannot satisfy the stability conditions, and it is only at four dimensions that all rules are simultaneously satisfied.  
(§12)

**27. Origin of the Time Axis**  
An odd number of axis inversions is required to realize negative volume. The most economical single-axis inversion generates $(-1)^1 = -1$. This "1-to-3 asymmetry" gives rise to the structure of one time axis and three spatial axes in the Minkowski metric $(-,+,+,+)$.  
(§13)

---

## Structure of the Framework

```
Central Projection Φ: Π_R → Sⁿ(R)
│
├── Basic Structure (Paper [1])
│    ├── Isomorphism with de Sitter spacetime
│    ├── Regularity of forward and inverse projection
│    ├── Induced metric and Einstein tensor
│    └── Flat limit (R → ∞)
│
├── Four Geometric Symmetries (Paper [2])
│    ├── Symmetry I: Discrete stability
│    ├── Symmetry II: Axis equivalence
│    ├── Symmetry III: Indistinguishability
│    └── Symmetry IV: Transformability of subjective coordinates
│
├── Multiple Subjective Spaces and Relativity of Observation (Paper [3])
│    ├── n+1 subjective spaces simultaneously constructible
│    ├── Internal observers access only their own subjective space
│    └── Axis roles are exchanged by the observer's standpoint
│
├── Limits of Curvature Radius and Nested Structure (Paper [4])
│    ├── Regularity of the flat limit (R → ∞)
│    ├── Structural tension at R → 0
│    └── Nested structure of subjective spaces
│
├── Dimension Addition (Paper [5])
│    ├── Type I (under existing R) and Type II (new R')
│    └── Type II requires simultaneous addition to avoid 0D degeneration
│
├── Subjective Spaces at Each Dimension (Paper [6])
│    └── Geometric characteristics from 0D to 4D
│
├── Continuity of Geodesics (Paper [7])
│    └── Parent-child geodesics are independent geometric objects
│
├── Schwarzschild–de Sitter Exact Solution (Paper [9])
│    ├── Uniquely determined from Paper [1]'s equation without additional assumptions
│    ├── Horizons merge and observer region vanishes at R → 0
│    └── Horizon structure recovers in child subjective space
│
├── Dimensional Interpretation of R → 0 Limit (Paper [10])
│    ├── R → 0 corresponds to zero-dimensional subjective space
│    └── Dimensional collapse → Dimensional reconstruction (two-stage description)
│
└── Singularity of Four Dimensions (Paper [8])
     ├── Integer property of circumscribed hypersphere diameter
     ├── Smallest non-trivial dimension satisfying stability conditions
     └── 1-to-3 asymmetry and origin of the time axis
```

---

## What This Research Does Not Claim

This research is confined to purely geometric propositions and their proofs. The following claims are not made:

- That the internal observers of subjective spaces physically exist
- That the "relativity of observation" corresponds to any specific theory of relativity in physics
- That a particular subjective space corresponds to our universe
- That physical laws are derivable from the structure of subjective spaces

These physical interpretations require additional assumptions beyond the geometric propositions of this research (Paper [3], §7).

---

## References

![QR Paper 1](qr_paper1.png)
[1] Kihara, N. (2026). *Geometric Formulation of 4-Dimensional Spacetime via Gnomonic Projection*. DOI: 10.5281/zenodo.19427780.

&nbsp;

![QR Paper 2](qr_paper2.png)
[2] Kihara, N. (2026). *Geometric Symmetries of Central Projection: Mathematical Foundations of the Multi-Axis Model*. DOI: 10.5281/zenodo.19434932.

&nbsp;

![QR Paper 3](qr_paper3.png)
[3] Kihara, N. (2026). *Relativity of Observation in Multiple Subjective Spaces: Geometric Consequences of the Symmetries of Central Projection*. DOI: 10.5281/zenodo.19435162.

&nbsp;

![QR Paper 4](qr_paper4.png)
[4] Kihara, N. (2026). *On the Limits of the Curvature Radius in Central Projection*. DOI: 10.5281/zenodo.19526549.

&nbsp;

![QR Paper 5](qr_paper5.png)
[5] Kihara, N. (2026). *On Dimension Addition to Background Spaces and Subjective Spaces*. DOI: 10.5281/zenodo.19526913.

&nbsp;

![QR Paper 6](qr_paper6.png)
[6] Kihara, N. (2026). *Subjective Spaces from Zero to Four Dimensions*. DOI: 10.5281/zenodo.19533292.

&nbsp;

![QR Paper 7](qr_paper7.png)
[7] Kihara, N. (2026). *On the Continuity of Geodesics in Subjective Spaces*. DOI: 10.5281/zenodo.19533299.

&nbsp;

![QR Paper 8](qr_paper8.png)
[8] Kihara, N. (2026). *Diameter of the Circumscribed Hypersphere of a Unit Four-Dimensional Hyperrectangle*. DOI: 10.5281/zenodo.19533313.

&nbsp;

[9] Kihara, N. (2026). *Schwarzschild–de Sitter Exact Solution in the Central Projection Framework*. DOI: 10.5281/zenodo.19538098.

&nbsp;

[10] Kihara, N. (2026). *Dimensional Interpretation of Geodesic Structure at the R→0 Limit*. DOI: 10.5281/zenodo.19538106.
