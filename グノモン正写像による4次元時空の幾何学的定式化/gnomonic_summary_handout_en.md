# Geometric Formulation of Spacetime via Central Projection: Overview of Three Papers

**Author:** Noriaki Kihara  
**Affiliation:** WF System Co., Ltd. (Osaka University, Faculty of Engineering Science, B.Eng.)  
**ORCID:** 0009-0004-6753-4020  
**Date:** April 2026

---

## Papers

| # | Title | DOI |
|---|-------|-----|
| [1] | Geometric Formulation of 4-Dimensional Space via Central Projection | [10.5281/zenodo.19427780](https://doi.org/10.5281/zenodo.19427780) |
| [2] | Geometric Symmetries of Central Projection: Mathematical Foundations of the Multi-Axis Model | [10.5281/zenodo.19434932](https://doi.org/10.5281/zenodo.19434932) |
| [3] | Geometric Role of the Fifth Dimension in Discrete Central Projection: Charge Quantization and the $1/r^2$ Law | [10.5281/zenodo.19435162](https://doi.org/10.5281/zenodo.19435162) |

| Paper [1] | Paper [2] | Paper [3] |
|:-:|:-:|:-:|
| ![QR Paper 1](qr_paper1.png) | ![QR Paper 2](qr_paper2.png) | ![QR Paper 3](qr_paper3.png) |

---

## Purpose of This Research

From purely geometric properties of central projection from an $n$-dimensional tangent hyperplane to a hypersphere $S^n(R)$, we construct a unified framework in which multiple independent axes that generate acceleration can coexist without contradiction. We refrain from physical identifications and confine ourselves to geometric propositions and their proofs.

---

## Proven Geometric Propositions

### Paper [1]: Basic Structure of Central Projection

**1. Isomorphism with de Sitter Spacetime**  
The induced metric of the Lorentzian central projection coincides numerically with the metric in the Beltrami coordinate system of de Sitter spacetime, establishing complete equivalence as intrinsic geometry.  
(Section 1.1, Appendix C, Proposition C.1)

**2. Regularity of Forward and Inverse Projection**  
Since $l = \sqrt{R^2 + r^2} \geq R > 0$, the central projection $\Phi$ is $C^\infty$-regular on its entire domain. The inverse $x^\mu = RY^\mu/Y^{n+1}$ is regular on the upper hemisphere $H^+$. The mapping is bijective.  
(Appendix B, Proposition B.1)

**3. Induced Metric and Einstein Tensor**  
The pullback metric $g_{\mu\nu} = (R^2/l^2)(\delta_{\mu\nu} - x_\mu x_\nu/l^2)$ is rigorously derived, and the $n$-dimensional Einstein tensor satisfies $G_{\mu\nu} + \Lambda_n g_{\mu\nu} = 0$ with $\Lambda_n = (n{-}1)(n{-}2)/(2R^2)$.  
(Sections 3, 6, Proposition 6.1)

**4. Flat Limit**  
As $R \to \infty$: $g_{\mu\nu} \to \delta_{\mu\nu}$, $\Lambda \to 0$. The central projection degenerates to the identity map, recovering Minkowski spacetime.  
(Appendix A)

---

### Paper [2]: Five Geometric Symmetries

**5. Symmetry I (Discrete Stability)**  
For any coordinate value — zero, positive integer, or negative integer — $l_A > 0$ holds and the central projection is globally regular without breakdown.  
(Section 3, Theorem 3.1)

**6. Symmetry II (Axis Equivalence)**  
Regardless of which axis in $\mathbb{R}^{n+1}$ is chosen as the projection center, the metric, curvature tensor, and Einstein tensor have identical mathematical structure.  
(Section 4, Theorem 4.1)

**7. Symmetry III (Geodesic Deviation and Indistinguishability)**  
The geodesic deviation $|\xi|^2 = g_{\mu\nu}\xi^\mu\xi^\nu$ is a scalar invariant of dimension $[L^2]$ and the sole means by which an internal observer recognizes curvature. However, neither the sign of coordinates nor the direction of curvature can be distinguished. The direction of acceleration is defined only in background coordinates.  
(Section 5, Theorem 5.2)

**8. Symmetry IV (Transformability of Subjective Coordinate Systems)**  
Subjective coordinate systems centered on different axes are mutually transformable via $T_{A \to B} = \Phi_B^{-1} \circ \Phi_A$, a $C^\infty$ diffeomorphism. The only parameters are the curvature radius $R$ and coordinate values. Background coordinates are preserved.  
(Section 6, Theorem 6.1, Proposition 6.5)

**9. Symmetry V (Centripetal Acceleration on Great Circles)**  
An internal observer in uniform circular motion along a great circle (geodesic) measures centripetal acceleration $a = R\omega^2$ as a scalar. The direction is indistinguishable.  
(Section 7, Theorems 7.1 and 7.2)

**10. Curvature Radius as Scale Factor**  
The Jacobian of coordinate transformations between subjective systems is $\partial x'^A/\partial x^B = -R^2/(x^B)^2$, with a dimensionless ($[L]/[L]$) first-order scale factor. The curvature radius $R$ serves as the sole transformation parameter.  
(Section 6.4, Proposition 6.2)

---

### Paper [3]: Addition of the Fifth Axis and $1/r^2$ Structure

**11. Number-Theoretic Requirement for Discrete Isotropy**  
The minimum spatial dimension for isotropic null geodesics on a discrete lattice $\mathbb{Z}^{1+N}$ is $N = 4$. For $N = 3$, $\sqrt{3} \notin \mathbb{Q}$ precludes solutions; for $N = 4$, $t = 2x$ provides solutions.  
(Section 2, Proposition 2.1)

**12. Quantization from the Discrete Lattice**  
From the discretization condition $\tilde{\phi} \in \mathbb{Z}$ alone, $q = n \cdot e_0$ follows algebraically. No compactification is required.  
(Section 3, Proposition 3.1)

**13. Geodesic Deviation and $1/r^2$ Scaling along the Fifth Axis**  
With the fifth axis $\phi$ as the projection center, geodesic deviation in the subjective space scales as

$$|\xi|^2 \sim \frac{\phi_0^2}{r^4} \cdot |\xi_0|^2$$

This is a positive quantity independent of the sign of $\phi_0$ and the direction of $r$.  
(Section 5, Proposition 5.1, Remark 5.2)

---

## Structure of the Framework

The geometric propositions above combine to yield the following structure:

```
Central Projection Φ: Π_R → Sⁿ(R)
│
├── Axis Equivalence (Symmetry II)
│    Any axis may be chosen as the projection center; structure is identical
│
├── Transformability of Subjective Coordinates (Symmetry IV)
│    Subjective coordinates for different axes are mutually
│    transformable with R as the sole parameter
│
├── Discrete Stability (Symmetry I)
│    Regular even for integer coordinates
│
├── Indistinguishability (Symmetry III, Theorem 5.2)
│    Curvature direction and sign are invisible from within
│
├── Measurability of Acceleration (Symmetry V)
│    Centripetal acceleration a = Rω² is measurable as a scalar
│
└── Addition of the Fifth Axis (Paper [3])
     ├── Number-theoretic requirement for discrete isotropy → N=4
     ├── Discrete lattice → quantization
     └── Geodesic deviation ∝ φ₀²/r⁴ (1/r² structure)
```

This framework possesses a structure in which multiple independent axes generating acceleration (axes with different curvature radii) can coexist without geometric breakdown. Subjective coordinates for each axis are mutually transformable via Symmetry IV, and all axes are equivalent by Symmetry II.

---

## What This Research Does Not Claim

This research is confined to purely geometric propositions and their proofs. The following physical identifications are not made:

- That the fifth axis $\phi$ is identical to the physical "charge" dimension
- That $|\xi|^2 \propto \phi_0^2/r^4$ is isomorphic to Coulomb's law
- That Maxwell's equations are derivable from this framework
- That unification of electromagnetic and gravitational fields has been achieved

These physical identifications require additional assumptions beyond the geometric propositions of this research (Paper [3], Section 6).

---

## References

[1] Kihara, N. (2026). *Geometric Formulation of 4-Dimensional Spacetime via Gnomonic Projection*. DOI: 10.5281/zenodo.19427780.

[2] Kihara, N. (2026). *Geometric Symmetries of Central Projection: Mathematical Foundations of the Multi-Axis Model*. DOI: 10.5281/zenodo.19434932.

[3] Kihara, N. (2026). *Geometric Role of the Fifth Dimension in Discrete Central Projection: Charge Quantization and the $1/r^2$ Law*. DOI: 10.5281/zenodo.19435162.
