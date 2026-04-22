# Regularity of Central Projection Between Hyperspherical Shells: A Formulation Without Tangent Hyperplanes and Degeneracy at $R = 0$

**Author**: Noriaki Kihara  
**Affiliation**: WF System Co., Ltd.  
**ORCID**: 0009-0004-6753-4020  
**Version**: v1.0  
**Date**: April 2026  
**DOI**: [10.5281/zenodo.19692192](https://doi.org/10.5281/zenodo.19692192)

---

## Abstract

In the classical formulation [P1], central projection is defined via a tangent hyperplane, where the divergence of $\tan\theta$ at the equator necessitates the constraint $R > 0$. This paper proves that central projection defined directly between surfaces of concentric hyperspherical shells does not require this constraint.

Specifically, we show:

1. The interspherical central projection $\Phi_{R \to R'} : S^n(R) \to S^n(R')$ is a diffeomorphism for any $R, R' > 0$, regular and non-divergent over the entire sphere (Theorem 1).
2. The projection to $R' = 0$ is well-defined and smooth, with no divergence. The only limitation is that the Jacobian degenerates to zero and the inverse map becomes undefined (Theorem 2).
3. The classical gnomonic projection via a tangent hyperplane genuinely diverges at the equator. This divergence is intrinsic to the use of the tangent hyperplane and does not structurally exist in interspherical projection (Theorem 3).

These results confirm that the constraint $R, R_1, R_0 > 0$ in the three-layer structure model [M1] is unnecessary. All claims of [M1] hold over the entire real domain including $R = 0$, where the only limitation is the undefined inverse map (degeneracy).

**Keywords**: central projection, interspherical projection, regularity, gnomonic projection, equatorial divergence, degeneracy, diffeomorphism

---

## §1 Introduction and Motivation

### §1.1 Limitations of Formulation via Tangent Hyperplanes

In [P1], central projection was defined as the projection from the origin $O$ of the ambient coordinate space $\mathbb{R}^{n+1}$ onto the hypersphere $S^n(R)$:

$$
P_R : \mathbb{R}^{n+1} \setminus \{O\} \to S^n(R), \quad P_R(\mathbf{x}) = R \cdot \frac{\mathbf{x}}{\|\mathbf{x}\|}
$$

This definition is meaningful only for $R > 0$. The requirement $R > 0$ arises for two reasons:

1. **Formal reason**: When $R = 0$, $S^n(0) = \{O\}$ degenerates and the image becomes a single point.
2. **Pathological properties of tangent hyperplanes**: In the formulation of [P1], the pullback metric is constructed via the tangent hyperplane $T$ at the north pole of $S^n(R)$. The inverse map through this tangent hyperplane diverges as $\tan\theta \to \infty$ at the equator ($\theta = \pi/2$) and does not reach the southern hemisphere.

This paper proves that the second reason is **a limitation arising from the use of the tangent hyperplane, which does not exist in interspherical central projection**, and that the first reason is merely a limited restriction of "undefined inverse map."

### §1.2 Observation on the Three-Layer Structure Model

[M1] defined the central projection $\Phi_{1 \to 0}$ between hyperspheres $S^n(R_1)$ and $S^n(R_0)$ of different curvature radii and asserted that this map is well-defined and regular. All claims of [M1] use only interspherical projection and never employ tangent hyperplanes.

This observation raises a natural question: Is the constraint $R > 0$ in [M1] essential, or merely a remnant from the formulation of [P1]?

This paper proves that the constraint is merely a remnant, and that interspherical central projection is regular over the entire real domain including $R = 0$.

---

## §2 Definition of Central Projection Between Hyperspherical Shells

### Definition 2.1 (Interspherical Central Projection)

For $R, R' > 0$, the central projection between concentric hyperspheres $S^n(R)$ and $S^n(R')$ is defined as:

$$
\Phi_{R \to R'} : S^n(R) \to S^n(R'), \quad \Phi_{R \to R'}(p) = \frac{R'}{R} \cdot p
$$

where $p \in S^n(R) \subset \mathbb{R}^{n+1}$ with $\|p\| = R$.

**Remark**: This map scales $p$ along the half-line through $O$ to distance $R'$. It is equivalent to $\Phi_{1 \to 0} = P_{R_0} \circ P_{R_1}^{-1}$ from [M1] §5.3, but here we use the direct representation.

### Definition 2.2 (Extension to $R' = 0$)

For $R > 0$, $R' = 0$, we define:

$$
\Phi_{R \to 0} : S^n(R) \to \{O\}, \quad \Phi_{R \to 0}(p) = O \quad (\forall p \in S^n(R))
$$

This is the constant map naturally obtained as the limit $R' \to 0$ in Definition 2.1.

### Definition 2.3 (Degeneracy)

A map $\Phi : X \to Y$ is said to be **degenerate** at $x_0 \in X$ if $\Phi$ is smooth in a neighborhood of $x_0$ but the Jacobian vanishes at $x_0$. Degeneracy does not signify loss of smoothness (singularity) but rather the local undefined nature of the inverse map.

---

## §3 Main Results

### Theorem 1 (Diffeomorphism of Interspherical Projection)

For $R, R' > 0$, $\Phi_{R \to R'} : S^n(R) \to S^n(R')$ is a diffeomorphism. In particular:

(i) $\Phi_{R \to R'}$ is well-defined, smooth, and bijective over the entire sphere.  
(ii) The inverse map $\Phi_{R \to R'}^{-1} = \Phi_{R' \to R}$ is also smooth.  
(iii) The Jacobian is the constant $(R'/R)^n$ over the entire sphere and never vanishes.  
(iv) The image is all of $S^n(R')$; no hemisphere restriction exists.

**Proof**:

Represent points on $S^n(R)$ as $p = R \cdot \hat{u}$ ($\hat{u} \in S^n(1)$, $\|\hat{u}\| = 1$).

(i) $\Phi_{R \to R'}(p) = (R'/R) \cdot p = R' \cdot \hat{u}$. Since $\|R' \cdot \hat{u}\| = R'$, the image lies on $S^n(R')$. The image is unique for each $\hat{u}$, hence well-defined; being a scalar multiple, it is smooth. $\hat{u} \neq \hat{u}' \Rightarrow R' \hat{u} \neq R' \hat{u}'$ gives injectivity; for any $q = R' \hat{v} \in S^n(R')$, $p = R \hat{v}$ provides the preimage, giving surjectivity.

(ii) $\Phi_{R' \to R}(q) = (R/R') \cdot q$, and $\Phi_{R' \to R} \circ \Phi_{R \to R'} = \mathrm{id}_{S^n(R)}$ is verified by direct computation. Since $R/R'$ is a finite nonzero constant, the inverse map is also smooth.

(iii) Since $S^n(R)$ is the $R$-scaled copy of $S^n(1)$, using local coordinates $(\theta_1, \ldots, \theta_n)$ on the unit sphere $S^n(1)$, the induced metric on $S^n(R)$ is $g_{ij}^{(R)} = R^2 \cdot g_{ij}^{(1)}$. Since $\Phi_{R \to R'}$ induces the identity map $\hat{u} \mapsto \hat{u}$ on the unit sphere, the Jacobian matrix of $\Phi_{R \to R'}$ gives a scaling of $R'/R$ in each tangent direction. The $n$-dimensional Jacobian (ratio of volume elements) is $(R'/R)^n$, which never vanishes for $R, R' > 0$.

(iv) Directly from surjectivity in (i). Since $\hat{u}$ ranges over the entire unit sphere, the image covers all of $S^n(R')$. $\blacksquare$

---

### Theorem 2 (Regularity and Degeneracy at $R' = 0$)

For $R > 0$, the following hold for $\Phi_{R \to 0} : S^n(R) \to \{O\}$:

(i) $\Phi_{R \to 0}$ is well-defined and smooth.  
(ii) The image is bounded over the entire sphere ($\|\Phi_{R \to 0}(p)\| = 0$, $\forall p$), and no $\tan\theta$-type divergence as in gnomonic projection exists.  
(iii) The Jacobian degenerates to $0$ over the entire sphere.  
(iv) $\Phi_{R \to 0}$ is not bijective, and the inverse map $\Phi_{0 \to R}$ is undefined.

**Proof**:

(i) $\Phi_{R \to 0}(p) = (0/R) \cdot p = \mathbf{0} = O$ is a constant map, clearly well-defined and smooth ($C^\infty$).

(ii) For any $p \in S^n(R)$, $\|\Phi_{R \to 0}(p)\| = 0 < \infty$, so the image is bounded. The unbounded behavior $d(\theta) = R\tan\theta \to \infty$ of gnomonic projection (Theorem 3) does not structurally exist.

(iii) In the expression from Theorem 1 (iii), the Jacobian is $(R'/R)^n$. For $R' = 0$, $(0/R)^n = 0$. The degeneracy of the Jacobian to zero reflects that the map collapses all points to a single point; it is not a divergence.

(iv) $\Phi_{R \to 0}$ maps all points of $S^n(R)$ to $O$, so it is not injective. Therefore, the inverse map does not exist. Specifically, the preimage of $O$ is the entirety of $S^n(R)$, making it impossible to determine which $p \in S^n(R)$ should be recovered from $O$. $\blacksquare$

---

### Theorem 3 (Comparison with Projection via Tangent Hyperplane)

The inverse map $G^{-1} : S^n(R) \to T_N$ of the classical gnomonic projection

$$
G : T_N \to S^n(R), \quad G(\mathbf{y}) = R \cdot \frac{\mathbf{y}}{\|\mathbf{y}\|}
$$

via the tangent hyperplane $T_N$ at the north pole $N = (0, \ldots, 0, R)$ of $S^n(R)$, diverges at the equator (zenith angle $\theta = \pi/2$).

**Proof**:

Represent points on $S^n(R)$ by the zenith angle $\theta$ (angular distance from the north pole). The distance from the north pole to the corresponding point on the tangent hyperplane $T_N$ is

$$
d(\theta) = R \tan\theta
$$

As $\theta \to \pi/2$, $\tan\theta \to +\infty$, so $d(\theta) \to +\infty$. Therefore, $G^{-1}$ diverges at the equator.

Furthermore, for $\theta > \pi/2$ (southern hemisphere), $G^{-1}$ is undefined. Since the tangent hyperplane $T_N$ lies only in the half-space on the north pole side, no preimage on $T_N$ exists for points on the southern hemisphere. Even extending the tangent hyperplane to both sides does not resolve the equatorial divergence, nor does it achieve coverage of the southern hemisphere. This is because the divergence is intrinsically due to the geometry of the tangent hyperplane (its placement as a tangent space of the hypersphere).

**Comparison**: The interspherical projection $\Phi_{R \to R'}$ (Theorem 1) does not structurally possess this divergence. $\Phi_{R \to R'}$ is a scalar multiplication and does not involve unbounded functions such as $\tan\theta$. The Jacobian is uniformly $(R'/R)^n$ over the entire sphere and has no dependence on $\theta$. $\blacksquare$

---

## §4 Consequences for Prior Works

### §4.1 Removal of the $R > 0$ Constraint in the Three-Layer Structure Model [M1]

In [M1], the condition $R > 0$ was imposed on the three curvature radii $R$, $R_1$, $R_0$ ([M1] §2 Condition 1, §3.1). This constraint was inherited from the formulation of [P1] ($R > 0$).

By Theorems 1 and 2 of this paper, the claims of [M1] are extended as follows:

**Corollary 4.1**: Claims 4.1 (shared origin), 4.2 (collinearity of projection axes), 4.3 (concentric nested structure), and 5.1 (metric non-contribution of $R$-directional displacement) of [M1] hold over the entire domain $R, R_1, R_0 \geq 0$. The only limitation is that the inverse map $\Phi_{0 \to R}$ becomes undefined for projections involving $R = 0$.

**Proof**: Each claim of [M1] relies solely on properties of interspherical central projection and does not use tangent hyperplanes. By Theorem 1, claims for $R, R' > 0$ hold as stated. For $R' = 0$, Theorem 2 (i)(ii) shows that the forward projection is well-defined and regular, so the claims of shared origin, collinearity, and concentric structure are unaffected. Metric non-contribution (Claim 5.1) holds trivially since the tangent space of $S^n(0) = \{O\}$ degenerates to zero dimension, making the orthogonality assertion vacuously true. $\blacksquare$

**Remark**: This corollary is a mathematical extension of scope and does not override the non-claim explicitly stated in [M1] §8.1 that "we do not assert the picture of particles condensing at the $R \to 0$ limit." The regularity at $R = 0$ means only that the map is mathematically well-defined; it does not assert physical condensation phenomena. Even in the cases $R_1 = 0$ (particles condensing to the origin) or $R_0 = 0$ (our spacetime condensing to the origin), the forward projection is mathematically well-defined as a constant map, but their physical interpretations lie outside the scope of this paper and [M1].

### §4.2 Relationship with the Formulation of [P1]

The formulation via tangent hyperplanes in [P1] was essential for constructing the pullback metric and deriving the Einstein tensor. This paper does not claim that the formulation of [P1] is erroneous.

What this paper claims is that **when using only interspherical central projection** (the framework of [M1]), the $R > 0$ constraint of [P1] is unnecessary. The use of the tangent hyperplane is required for the specific purposes of [P1] (metric construction) and is unnecessary for the geometric structures of interspherical projection (shared origin, nested structure, direction preservation).

---

## §5 Conclusion

This paper has proved that central projection defined between surfaces of concentric hyperspherical shells possesses none of the pathological properties of the classical gnomonic projection (projection via tangent hyperplanes).

| Property | Via tangent hyperplane (gnomonic) | Interspherical projection (this paper) |
|----------|----------------------------------|---------------------------------------|
| Regular domain | Northern hemisphere only ($\theta < \pi/2$) | Entire sphere ($\theta \in [0, \pi]$) |
| Behavior at equator | Diverges as $\tan\theta \to \infty$ | No divergence |
| Jacobian | Depends on $\theta$, diverges at equator | Constant $(R'/R)^n$, uniform over sphere |
| $R' = 0$ | Undefined as projection target (tangent hyperplane cannot be constructed) | No divergence; only inverse map undefined |
| Bijectivity | Limited to northern hemisphere | $R' > 0$: bijective over entire sphere; $R' = 0$: not bijective |

Interspherical central projection is regular (non-divergent) over the entire domain $R \geq 0$. The only limitation at $R = 0$ is the undefined inverse map due to degeneracy of the Jacobian to zero.

These results confirm that the constraint $R > 0$ in the three-layer structure model [M1] is unnecessary. The geometric claims of [M1] hold over the entire domain including $R = 0$.

---

## References

- [P1]: Noriaki Kihara, Geometric Formulation of Four-Dimensional Spacetime via Central Projection, Zenodo, 2026. DOI: [10.5281/zenodo.19427780](https://doi.org/10.5281/zenodo.19427780)
- [M1]: Noriaki Kihara, Three-Layer Structure Model of Central Projection: Nested Structure of Subjective Spaces via $R$–$R_1$–$R_0$ and Middleware Geometry, Zenodo, 2026. DOI: [10.5281/zenodo.19691713](https://doi.org/10.5281/zenodo.19691713)

---

*v1.0*
