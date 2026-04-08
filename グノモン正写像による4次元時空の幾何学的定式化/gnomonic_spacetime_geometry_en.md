# Geometric Formulation of 4-Dimensional Space via Central Projection

**Author:** Noriaki Kihara  
**Affiliation:** WF System Co., Ltd. (B.Eng., Osaka University, School of Engineering Science)  
**Date:** April 2026  
**Category:** Research Note (Geometric Considerations)  
**DOI:** [10.5281/zenodo.19427780](https://doi.org/10.5281/zenodo.19427780)

---

## Abstract

This report examines, from a purely geometric perspective, the deformation of the induced metric that arises in the coordinate system after mapping a 4-dimensional Euclidean space (tangent hyperplane $\Pi_R$) onto the surface of a 5-dimensional hypersphere of curvature radius $R$ via central projection. The central projection employed here takes the center of the sphere as the projection center and the tangent hyperplane as the projection surface, corresponding to the special central projection known in cartography as the gnomonic projection [1]. In particular, we demonstrate step by step: (i) the exact derivation of the induced metric tensor $g_{\mu\nu}$ after the mapping (including the contribution of the 5th component), and (ii) the geometric relationship between the Einstein tensor and the constant $\Lambda = 3/R^2$. This report excludes all cosmological assumptions.

---

## 1. Introduction

Central projection is defined as the mapping that projects onto a tangent plane from the center of a sphere as the projection center. In cartography, it is known as the gnomonic projection [1], and its most important property is that "great circles (geodesics) on the sphere are mapped to straight lines on the tangent plane."

In this report, we extend this central projection from 4 dimensions to 5 dimensions and investigate the metric structure that arises in the post-mapping coordinate system. This problem setting is purely geometric, and physical interpretations (gravity, spacetime, relativity, etc.) are outside the scope of this report.

---

### 1.1 On the Equivalence of Descriptions

In this paper, we describe spacetime using the central projection construction, but the Lorentzian version ($\epsilon_4 = -1$) of the induced metric (C.3) numerically coincides with the metric in the Beltrami coordinate system of de Sitter spacetime in the same coordinate representation (Appendix C, and prior works [8, 9, 10]). Therefore, both are Lorentzian manifolds possessing the same differentiable manifold and the same metric $g_{\mu\nu}$, and are completely equivalent as intrinsic geometry.

The essence of this equivalence can be understood from the perspective of an observer residing within the manifold. What an internal observer can measure is only the metric $g_{\mu\nu}$ and the curvature derived from it; there is in principle no intrinsic means to distinguish whether the spacetime was "constructed from a tangent plane via central projection" or "directly given as de Sitter spacetime." The existence of the background space, the choice of tangent plane, and the construction method of the mapping are all extrinsic information and do not appear in the intrinsic geometry [3, 6].

Physical quantities in general relativity---geodesics, curvature tensors, covariant derivatives, solutions of the field equations, physical observables---are all determined solely by the metric $g_{\mu\nu}$ and its derivatives, and do not depend on how the manifold is concretely described (in which coordinate system it is written, with which embedding it is visualized, or by which construction method it was obtained). This is a direct consequence of general covariance and diffeomorphism invariance [4, 5], and is the principle that pervades this paper.

Therefore, any calculation performed using the central projection description in this paper is mathematically equivalent to a calculation performed using the de Sitter description under the same metric, and the results of both coincide identically. There is no need to prove the equivalence for each individual physical quantity; it follows automatically from the coincidence of the metric (Appendix C) and general covariance (the above principle). Henceforth in this paper, we adopt the central projection description for its computational transparency, but all conclusions obtained hold equally in the de Sitter description.

### 1.2 On the Premise of Lorentzian Signature

Although this paper does not aim at physical interpretation, we clarify the logical status of the premise here to avoid misunderstanding regarding the origin of the Lorentzian signature.

In the Lorentzian version of the central projection (Appendix C, $\epsilon_4 = -1$), a negative sign appears in the squared term of the time component of the induced metric $g_{\mu\nu}$. Regarding the origin of this negative sign, we organize the logical structure of this paper as follows.

The Lorentzian signature $(-,+,+,+)$---that is, the negativity of the metric component in the time direction---is not derived as a consequence of de Sitter spacetime or the central projection construction. It is a fundamental premise of the spacetime manifold, required by the empirical facts of causal structure (the existence of light cones) and the invariance of the speed of light [4, 5]. Flat Minkowski spacetime, the Schwarzschild solution, FLRW cosmology, and de Sitter spacetime all share this same Lorentzian signature, and the negativity of the signature is not a property specific to any particular solution.

In the central projection construction, this premise is expressed as the metric of the 5-dimensional embedding space $\eta_{AB} = \mathrm{diag}(-1,+1,+1,+1,+1)$. Since the induced metric is determined by the pullback $g_{\mu\nu} = \eta_{AB}\,\partial_\mu X^A\,\partial_\nu X^B$, the Lorentzian signature of the embedding space is directly inherited by $g_{\mu\nu}$. This does not provide the physical origin of the signature, but is a mathematical confirmation that the premise is consistently maintained throughout all stages of the construction.

Therefore, the logical structure of this paper is organized as follows:

1. That spacetime has Lorentzian signature is adopted as a **physical premise** (it is not a subject of derivation within the framework of this paper)
2. The central projection construction embodies this premise as $\eta_{AB}$ of the embedding space, and inherits it to $g_{\mu\nu}$ through the pullback
3. That the resulting $g_{\mu\nu}$ coincides with the de Sitter metric is derived as a **result** of the construction (Appendix C)

The physical origin of the signature lies in causality and the invariance of the speed of light, which is a premise external to this construction.

---

## 2. Preliminaries: Definition and Basic Properties of Central Projection

### 2.1 Central Projection in $n$ Dimensions

**Definition 2.1** ($n$-Dimensional Central Projection) [1, 2]  
In $n$-dimensional Euclidean space, consider a tangent hyperplane $\Pi_R$ (an $n$-dimensional Euclidean plane with coordinates $(x^1, \ldots, x^n) \in \mathbb{R}^n$, with origin at the tangent point) at distance $R$ (curvature radius) from the center of the sphere. Setting $l = \sqrt{R^2 + \sum_{i=1}^n (x^i)^2}$,

$$\Phi: \Pi_R \to S^n(R) \subset \mathbb{R}^{n+1}, \quad (x^1, \ldots, x^n) \mapsto \left(\frac{R x^1}{l}, \ldots, \frac{R x^n}{l},\ \frac{R^2}{l}\right) \tag{2.1}$$

is called the **central projection**. Here $S^n(R)$ is the $n$-dimensional sphere of radius $R$ embedded in $\mathbb{R}^{n+1}$:

$$S^n(R) = \left\{(Y^1, \ldots, Y^{n+1}) \in \mathbb{R}^{n+1} \;\middle|\; \sum_{A=1}^{n+1} (Y^A)^2 = R^2 \right\}$$

**Terminology**  
Throughout this paper, we fix the following terminology.

<img src="fig_central_projection_en.png" style="width:50%;" alt="Figure 2.1: Geometric construction of central projection">

**Figure 2.1: Geometric construction of central projection.** The origin $O$ of the background space $\mathbb{R}^{n+1}$ (projection center = center of the sphere) is at distance $R$ in the $e_{n+1}$ direction from the tangent point $P = Re_{n+1}$. The tangent hyperplane $\Pi_R$ is the $n$-dimensional plane with origin at $P$, orthogonal to $e_{n+1}$. The central projection $\Phi$ maps a point $x$ on $\Pi_R$ to the intersection $\Phi(x)$ of the half-line through $O$ with $S^n(R)$. The coordinates $(x^1, \ldots, x^n)$ all reside on $\Pi_R$. The projection center axis ($e_{n+1}$ direction) is not contained in $\Pi_R$.

| Term | Definition | Symbol |
|---|---|---|
| **Background space** | The Euclidean space containing the sphere $S^n(R)$. Not directly visible to an internal observer | $\mathbb{R}^{n+1}$ |
| **Tangent hyperplane** | The $n$-dimensional plane in the background space with origin at the tangent point $P = Re_{n+1}$, orthogonal to $e_{n+1}$. Coordinates $(x^1, \ldots, x^n)$ reside here | $\Pi_R$ |
| **Subjective space** | The image of the central projection $\Phi$. The upper hemisphere $\{Y^{n+1} > 0\}$ of the sphere $S^n(R)$. The space where the internal observer resides and measures distances using the pullback metric $g_{\mu\nu}$ | $H^+$ |
| **Internal observer** | A hypothetical observer who performs measurements within the subjective space $H^+$ using only the metric $g_{\mu\nu}$. Cannot directly access information from the background space $\mathbb{R}^{n+1}$ | --- |
| **Curvature radius** | The radius of the sphere $S^n(R)$. The sectional curvature is given by $K = 1/R^2$ | $R$ |
| **Projection center** | The origin of the background space $\mathbb{R}^{n+1}$ (center of the sphere) | $O$ |
| **Projection center axis** | The direction from the projection center $O$ to the tangent point $P$ ($e_{n+1}$ direction). Not contained in $\Pi_R$; unmeasurable from within as $[L]$ | --- |
| **Background coordinates** | Coordinates $(Y^1, \ldots, Y^{n+1})$ of a point on the sphere in the background space $\mathbb{R}^{n+1}$. Independent of both the curvature radius $R$ and the choice of projection center axis | $(Y^A)$ |
| **Subjective coordinates** | Coordinates $(x^1, \ldots, x^n)$ on the tangent hyperplane $\Pi_R$. For the same point on the sphere, the values change if the curvature radius $R$ changes, and change if the projection center axis changes | $(x^\mu)$ |

The central projection $\Phi: \Pi_R \to H^+$ converts subjective coordinates $(x^1, \ldots, x^n)$ to background coordinates $(Y^1, \ldots, Y^{n+1})$ (Eq. (2.1)). Conversely, the inverse transformation from background coordinates to subjective coordinates is given by $x^\mu = RY^\mu/Y^{n+1}$.

The internal observer resides in the subjective space $H^+$ and measures distances through the subjective coordinates $(x^\mu)$ and the pullback metric $g_{\mu\nu}$. The background coordinates $(Y^A)$ are information of the background space and are not directly visible to the internal observer. The direction of the projection center axis ($e_{n+1}$) is not contained in $\Pi_R$.

<img src="fig_two_projections_en.png" style="width:50%;" alt="Figure 2.3: Differences in subjective coordinates due to curvature radius and projection axis">

**Figure 2.3: Differences in subjective coordinates due to curvature radius and projection axis.** Two spheres of curvature radius $R_1$ (red, small circle) and $R_2$ (blue, large circle) are drawn concentrically. The intersection point $P_1 = (R_1, R_2)$ of the tangent hyperplanes $\Pi_{R_1}$ (red, vertical) and $\Pi_{R_2}$ (blue, horizontal) is the same point in the background space. However, the projection line from $O$ to $P_1$ reaches **different positions** on $S(R_2)$ and $S(R_1)$ ($\Phi_{R_2}(P_1) \neq \Phi_{R_1}(P_1)$). As shown by the thick arcs, the arc length (distance in subjective space) for the same angle $\theta$ is $R_1\theta$ and $R_2\theta$, and **the ratio is $R_1 : R_2$** (equal to the ratio of curvature radii). That is, even if the background coordinates are the same, if the curvature radii differ, the subjective coordinates and subjective distances also differ.

<img src="fig_subjective_space_en.png" style="width:50%;" alt="Figure 2.2: Coordinate grid of the tangent hyperplane and subjective space">

**Figure 2.2: Coordinate grid of the tangent hyperplane and subjective space.** The square grid on the tangent hyperplane $\Pi_R$ (top, coordinates $(x^1, x^2)$) is mapped by the central projection $\Phi$ to a curved grid on the subjective space $H^+$ (upper hemisphere of the sphere, bottom). The internal observer uses the coordinates $(x^\mu)$ of $\Pi_R$, but since distances are measured with the pullback metric $g_{\mu\nu}$, squares are distorted on the sphere. Near the tangent point $P$, the distortion is small, and it increases with distance from $P$.

**Remark 2.1** (Coverage and Sign of $R$)  
The central projection $\Phi$ of Definition 2.1 covers only the **upper hemisphere** $H^+ = \{Y \in S^n(R) \mid Y^{n+1} > 0\}$ of the sphere $S^n(R)$ (Figure 2.1). As $|x| \to \infty$ on $\Pi_R$, $\Phi(x)$ approaches the equator $\{Y^{n+1} = 0\}$ asymptotically but never reaches it (Appendix B, Corollary 3.2 [subsequent paper]).

Taking $R$ to be negative ($R < 0$), the tangent point moves to $P = Re_{n+1}$ (the opposite direction of $e_{n+1}$), and the tangent hyperplane $\Pi_R$ is positioned on the opposite side of the sphere. In this case, $\Phi$ covers the **lower hemisphere** $H^- = \{Y \in S^n(R) \mid Y^{n+1} < 0\}$ (see lower part of Figure 2.1). Combining the two central projections for $R > 0$ and $R < 0$ covers all of $S^n(R)$ except the equator. In what follows, we fix $R > 0$, but since the intrinsic geometry of $S^n(R)$ does not depend on the choice of hemisphere, no generality is lost.

**Proposition 2.1** (Great Circle Preservation) [1]  
Under central projection, any straight line on the tangent plane is mapped to a great circle (geodesic) on the sphere.

---

## 3. Derivation of the Metric Tensor

### 3.1 Setup of the Mapping

In this section, we set $n=4$ and consider the central projection from the tangent hyperplane $\Pi_R$ (coordinates $(x^1, x^2, x^3, x^4)$) to the 5-dimensional hypersphere $S^4(R) \subset \mathbb{R}^5$ of radius $R$. Since the tangent plane is Euclidean, we have $x_\mu = x^\mu$ ($\mu = 1,2,3,4$). By Definition 2.1, the 5-dimensional coordinates of the image are

$$(Y^1, Y^2, Y^3, Y^4, Y^5) = \left(\frac{Rx^1}{l},\ \frac{Rx^2}{l},\ \frac{Rx^3}{l},\ \frac{Rx^4}{l},\ \frac{R^2}{l}\right) \tag{3.1}$$

$$l = \sqrt{R^2 + \sum_{\mu=1}^4 (x^\mu)^2},\quad r^2 = \sum_{\mu=1}^4(x^\mu)^2,\quad l^2 = R^2 + r^2 \tag{3.2}$$

### 3.2 Derivation of the Pullback Metric

The induced metric of $S^4(R) \subset \mathbb{R}^5$ is determined from the standard inner product of $\mathbb{R}^5$ as

$$dS^2 = \sum_{A=1}^{5} (dY^A)^2 \tag{3.3}$$

[6]. **We compute the contributions of all 5 components, including the 5th component.**

**Components 1--4** $(Y^\mu = R x^\mu / l)$:

**Lemma 3.1**  
$$\frac{\partial Y^\mu}{\partial x^\nu} = \frac{R}{l}\left(\delta^\mu_\nu - \frac{x^\mu x^\nu}{l^2}\right) \tag{3.4}$$

**Proof**  
$$\frac{\partial}{\partial x^\nu}\!\left(\frac{Rx^\mu}{l}\right) = \frac{R\delta^\mu_\nu}{l} + Rx^\mu\!\left(-\frac{x^\nu}{l^3}\right) = \frac{R}{l}\!\left(\delta^\mu_\nu - \frac{x^\mu x^\nu}{l^2}\right) \quad \square$$

Therefore,

$$\sum_{\mu=1}^4(dY^\mu)^2 = \frac{R^2}{l^2}\sum_{\nu,\rho=1}^{4}\left(\delta_{\nu\rho} - \frac{2x_\nu x_\rho}{l^2} + \frac{r^2 x_\nu x_\rho}{l^4}\right)dx^\nu dx^\rho \tag{3.5}$$

**5th component** $(Y^5 = R^2/l)$:

$$\frac{\partial Y^5}{\partial x^\nu} = -\frac{R^2 x^\nu}{l^3} \tag{3.6}$$

$$\left(dY^5\right)^2 = \frac{R^4}{l^6}\left(\sum_\nu x^\nu dx^\nu\right)^2 = \frac{R^4}{l^6}\,x_\nu x_\rho\,dx^\nu dx^\rho \tag{3.7}$$

**Sum of all components:**

Collecting the coefficients of $x_\nu x_\rho$ (adding (3.5) and (3.7)):

$$-\frac{2R^2}{l^4} + \frac{R^2 r^2}{l^6} + \frac{R^4}{l^6} = -\frac{2R^2}{l^4} + \frac{R^2(r^2+R^2)}{l^6} = -\frac{2R^2}{l^4} + \frac{R^2 l^2}{l^6} = -\frac{2R^2}{l^4} + \frac{R^2}{l^4} = -\frac{R^2}{l^4}$$

Therefore,

$$dS^2 = \frac{R^2}{l^2}\delta_{\nu\rho}\,dx^\nu dx^\rho - \frac{R^2}{l^4}\,x_\nu x_\rho\,dx^\nu dx^\rho$$

$$\boxed{g_{\mu\nu} = \frac{R^2}{l^2}\left(\delta_{\mu\nu} - \frac{x_\mu x_\nu}{l^2}\right)} \tag{3.8}$$

**Remark 3.1** (Indispensability of the 5th Component)  
If the 5th component $(dY^5)^2$ (Eq. (3.7)) is omitted, the coefficient of $x_\nu x_\rho$ remains $-2R^2/l^4 + R^2 r^2/l^6$, and Eq. (3.8) cannot be obtained. The contribution of the 5th component $+R^4/l^6 = +R^2/l^4 - R^2 r^2/l^6$ fills this gap.

---

## 4. Riemann Curvature Tensor

**Theorem 4.1** (Curvature Tensor of a Constant Curvature Space) [3, 4, 5]  
The Riemann curvature tensor of the $n$-dimensional sphere $S^n(R)$ is given by

$$R_{\mu\nu\rho\sigma} = \frac{1}{R^2}(g_{\mu\rho}g_{\nu\sigma} - g_{\mu\sigma}g_{\nu\rho}) \tag{4.1}$$

This is a standard result for spaces of constant curvature [3, 4, 5].

---

## 5. Ricci Tensor and Einstein Tensor

### 5.1 Ricci Tensor

**Proposition 5.1**  
The Ricci tensor of the $n$-dimensional sphere $S^n(R)$ is

$$R_{\mu\nu} = \frac{n-1}{R^2}g_{\mu\nu} \tag{5.1}$$

**Proof**  
Substituting (4.1) into the definition of the Ricci tensor $R_{\mu\nu} = g^{\rho\sigma}R_{\mu\rho\nu\sigma}$:

$$g^{\rho\sigma}R_{\mu\rho\nu\sigma} = \frac{1}{R^2}g^{\rho\sigma}(g_{\mu\nu}g_{\rho\sigma} - g_{\mu\sigma}g_{\nu\rho})$$

Computing each term:

- **First term:** $g^{\rho\sigma}g_{\mu\nu}g_{\rho\sigma} = g_{\mu\nu}\cdot\underbrace{g^{\rho\sigma}g_{\rho\sigma}}_{=\,\delta^\rho{}_\rho\,=\,n} = n\,g_{\mu\nu}$
- **Second term:** $g^{\rho\sigma}g_{\mu\sigma}g_{\nu\rho} = \underbrace{g^{\rho\sigma}g_{\mu\sigma}}_{=\,\delta^\rho_\mu}g_{\nu\rho} = g_{\nu\mu} = g_{\mu\nu}$

Therefore,

$$R_{\mu\nu} = \frac{1}{R^2}(n\,g_{\mu\nu} - g_{\mu\nu}) = \frac{n-1}{R^2}g_{\mu\nu} \qquad \square$$

### 5.2 Ricci Scalar

$$R_{\text{scalar}} = g^{\mu\nu}R_{\mu\nu} = \frac{n-1}{R^2}\underbrace{g^{\mu\nu}g_{\mu\nu}}_{=\,n} = \frac{n(n-1)}{R^2} \tag{5.2}$$

For $n=4$, $R_{\text{scalar}} = 12/R^2$.

### 5.3 Einstein Tensor

From the definition of the Einstein tensor [4],

$$G_{\mu\nu} = R_{\mu\nu} - \frac{1}{2}g_{\mu\nu}R_{\text{scalar}} \tag{5.3}$$

For $n=4$,

$$G_{\mu\nu} = \frac{3}{R^2}g_{\mu\nu} - \frac{1}{2}\cdot\frac{12}{R^2}\,g_{\mu\nu} = \left(\frac{3}{R^2} - \frac{6}{R^2}\right)g_{\mu\nu} = -\frac{3}{R^2}g_{\mu\nu} \tag{5.4}$$

---

## 6. Main Result

**Proposition 6.1** (Main Result)  
When the tangent hyperplane $\Pi_R$ is mapped to the hypersphere $S^4(R)$ (subjective space $H^+$) of curvature radius $R$ via central projection, for the induced metric tensor $g_{\mu\nu}$ (Eq. (3.8)),

$$G_{\mu\nu} + \Lambda\, g_{\mu\nu} = 0 \tag{6.1}$$

holds, where

$$\boxed{\Lambda = \frac{3}{R^2}} \tag{6.2}$$

**Proof**  
From (5.4), $G_{\mu\nu} = -\dfrac{3}{R^2}g_{\mu\nu}$, so setting $\Lambda = \dfrac{3}{R^2}$ immediately yields (6.1). $\square$

**Remark 6.1**  
Equation (6.1) is a **geometric identity for $S^4(R)$ with Euclidean signature** and is not to be interpreted as a physical law. The relationship with de Sitter spacetime having Lorentzian signature is discussed in Appendix C.

**Remark 6.2** (Extension to General Dimensions)  
This proposition is not limited to $n=4$. For a general $n$-dimensional sphere $S^n(R)$, from (5.2) and (5.3),

$$G_{\mu\nu} = \frac{n-1}{R^2}g_{\mu\nu} - \frac{1}{2}\cdot\frac{n(n-1)}{R^2}g_{\mu\nu} = -\frac{(n-1)(n-2)}{2R^2}g_{\mu\nu}$$

and $\Lambda_n$ satisfying $G_{\mu\nu} + \Lambda g_{\mu\nu} = 0$ is

$$\Lambda_n = \frac{(n-1)(n-2)}{2R^2}$$

(for $n=4$, $\Lambda_4 = 3/R^2$).

---

## 7. Summary of Results

| Item                  | Result                                                                                      | Eq.    | Source  |
| --------------------- | ------------------------------------------------------------------------------------------- | ------ | ------- |
| Definition of mapping | $Y^\mu = Rx^\mu/l$, $Y^5 = R^2/l$                                                          | (3.1)  | [1, 2]  |
| Induced metric        | $g_{\mu\nu} = \dfrac{R^2}{l^2}\!\left(\delta_{\mu\nu} - \dfrac{x_\mu x_\nu}{l^2}\right)$   | (3.8)  | This report |
| Riemann curvature     | $R_{\mu\nu\rho\sigma} = \dfrac{1}{R^2}(g_{\mu\rho}g_{\nu\sigma}-g_{\mu\sigma}g_{\nu\rho})$ | (4.1)  | [3,4,5] |
| Ricci tensor          | $R_{\mu\nu} = \dfrac{3}{R^2}g_{\mu\nu}$                                                    | (5.1)  | [3,4,5] |
| Ricci scalar          | $R_{\text{scalar}} = \dfrac{12}{R^2}$                                                      | (5.2)  | [3,4,5] |
| Einstein tensor       | $G_{\mu\nu} = -\dfrac{3}{R^2}g_{\mu\nu}$                                                   | (5.4)  | [3,4,5] |
| **Main result**       | $G_{\mu\nu} + \Lambda g_{\mu\nu} = 0$, $\Lambda = \dfrac{3}{R^2}$                          | (6.2)  | This report |

---

## References

[1] Snyder, J. P., *Map Projections: A Working Manual*, U.S. Geological Survey Professional Paper 1395, U.S. Government Printing Office, Washington D.C., 1987, pp. 164--168.

[2] Coxeter, H. S. M., *Introduction to Geometry*, 2nd ed., Wiley, New York, 1969, Chapter 15.

[3] do Carmo, M. P., *Riemannian Geometry*, Birkhauser, Boston, 1992, Chapter 4.

[4] Misner, C. W., Thorne, K. S., Wheeler, J. A., *Gravitation*, W. H. Freeman, San Francisco, 1973, $\S$13.5.

[5] Wald, R. M., *General Relativity*, University of Chicago Press, Chicago, 1984, Appendix D.

[6] Spivak, M., *A Comprehensive Introduction to Differential Geometry*, Vol. 2, 3rd ed., Publish or Perish, Houston, 1999, Chapter 4.

[7] Gibbons, G. W. and Hawking, S. W., "Action integrals and partition functions in quantum gravity," *Physical Review D*, **15**, 2752--2756, 1977.  
(Foundations of Euclidean quantum gravity; the relationship between Euclidean de Sitter spacetime and $\Lambda = 3/R^2$)

[8] Guo, H.-Y., Huang, C.-G., Xu, Z., Zhou, B., "On Beltrami model of de Sitter spacetime," *Modern Physics Letters A*, **19**, 1701--1710, 2004. arXiv:hep-th/0405137.  
(Construction of de Sitter Beltrami coordinates via Lorentzian central projection; prior work on the $\epsilon_4=-1$ case of Appendix C)

[9] Guo, H.-Y., "The Beltrami Model of De Sitter Space," arXiv:hep-th/0607017, 2006.  
(Direct construction of the de Sitter metric via central projection and de Sitter invariant special relativity)

[10] Guo, H.-Y., Huang, C.-G., Xu, Z., Zhou, B., "Beltrami Model of Riemann-Sphere and de Sitter Spacetime," arXiv:gr-qc/0703078, 2007.  
(Comparison of Beltrami models of $S^n(R)$ (Euclidean) and $dS^n$ (Lorentzian); connection with the unified perspective of this appendix)

---

## Appendix A: Flat Limit as $R \to \infty$

As $R \to \infty$, $l \to R$, so

$$g_{\mu\nu} \to \delta_{\mu\nu}, \quad \Lambda = \frac{3}{R^2} \to 0$$

The metric converges to the flat Euclidean metric, and the sectional curvature $K = 1/R^2 \to 0$. This corresponds to the central projection degenerating to the identity map (for details of the limiting procedure, see [3] Chapter 4).

---

## Appendix B: Regularity of the Mapping

**Proposition B.1** (Regularity of Central Projection)  
Let $\Phi: \mathbb{R}^4 \to S^4(R)$ be the central projection of Definition 2.1. Setting $r = \sqrt{\sum_\mu (x^\mu)^2}$,

$$Y^5 = \frac{R^2}{l} = \frac{R^2}{\sqrt{R^2 + r^2}}, \quad Y^\mu = \frac{Rx^\mu}{\sqrt{R^2+r^2}}$$

the following hold:

1. **$R \to 0$ ($r > 0$ fixed):** $Y^5 = R^2/\sqrt{R^2+r^2} \to 0$, $Y^\mu \to 0$. The image of the mapping degenerates to the origin but does not diverge.

2. **$R \to \infty$ ($r$ fixed):** $Y^5 = R/\sqrt{1+r^2/R^2} \to \infty$. The 5th component diverges. Meanwhile, $Y^\mu \to x^\mu$ (finite). This is consistent with the image of $\Phi$ being localized near the north pole of $S^4(R)$, and the sphere approaching the tangent plane as $R \to \infty$ (the flat limit of Appendix A).

3. **Fixed $R > 0$:** Since $l = \sqrt{R^2+r^2} \geq R > 0$, the denominator is always positive, and $\Phi$ is a $C^\infty$ regular mapping on all of $\mathbb{R}^4$.

**Proof**  
All follow from direct computation with $l = \sqrt{R^2+r^2}$. (1) $R\to 0$: $Y^5 = R^2/l \leq R^2/r \to 0$. (2) $R\to\infty$: $Y^5 = R\cdot(1+r^2/R^2)^{-1/2} \to \infty$. (3) For $R>0$, $l \geq R > 0$, so each component is $C^\infty$ in $r$. $\square$

---

## Appendix C: Unified Description of Euclidean/Lorentzian Cases

This appendix generalizes the central projection construction of the main text using sign parameters $\epsilon_\mu \in \{+1,-1\}$, and derives formulas that treat the Euclidean and Lorentzian cases in a unified manner.

That the induced metric of the Lorentzian central projection ($\epsilon_4 = -1$) corresponds to the Beltrami coordinate system of de Sitter spacetime was established as prior work by Guo Han-Ying (H.-Y. Guo) et al. from a Lorentzian starting point [8, 9, 10]. The originality of this appendix lies in the **explicit derivation of formula (C.3), which unifies both the Euclidean case ($S^4(R)$, [7]) and the Lorentzian case (de Sitter Beltrami, [8, 9, 10]) via the single parameter $\epsilon_\mu$**.

**Proposition C.1** (Unified Euclidean/Lorentzian Description)  
Choose a signature $(\epsilon_1, \epsilon_2, \epsilon_3, \epsilon_4) \in \{+1,-1\}^4$.

Define the **metric of the background space** not as the standard metric of $\mathbb{R}^5$, but as

$$dS^2_{\rm emb} = \sum_{\mu=1}^{4} \epsilon_\mu (dY^\mu)^2 + (dY^5)^2 \tag{C.1}$$

(the coefficient of the 5th component $Y^5$ is fixed at $+1$).

Let the metric on the tangent hyperplane $\Pi_R$ be

$$ds^2_{\rm tan} = \sum_{\mu=1}^{4} \epsilon_\mu (dx^\mu)^2 \tag{C.2}$$

and set $l^2 = R^2 + \sum_{\mu=1}^4 \epsilon_\mu (x^\mu)^2$. The definition of the central projection $\Phi$ is the same as in the main text (2.1):

$$Y^\mu = \frac{R x^\mu}{l} \quad (\mu=1,2,3,4), \quad Y^5 = \frac{R^2}{l}$$

Then the induced metric pulled back by $\Phi$ is given by

$$g_{\mu\nu} = \frac{R^2}{l^2}\left(\epsilon_\mu \delta_{\mu\nu} - \frac{\epsilon_\mu \epsilon_\nu x_\mu x_\nu}{l^2}\right) \tag{C.3}$$

(no Einstein summation convention is used for repeated indices).

**Proof**  
We express $dS^2_{\rm emb}$ under (C.1) in terms of the tangent plane coordinates.

$dY^\mu$ is computed as in Lemma 3.1 of the main text. For each $\mu$, we expand $(dY^\mu)^2$ as a quadratic form in $dx^\nu dx^\rho$, and summing over $\mu$:

$$\sum_{\mu=1}^4 \epsilon_\mu (dY^\mu)^2 = \frac{R^2}{l^2} \sum_{\nu,\rho=1}^4 \left(\epsilon_\nu \delta_{\nu\rho} - \frac{2\epsilon_\nu \epsilon_\rho x_\nu x_\rho}{l^2} + \frac{\epsilon_\nu \epsilon_\rho r_\epsilon^2 x_\nu x_\rho}{l^4}\right) dx^\nu dx^\rho \tag{C.4}$$

where $r_\epsilon^2 = \sum_\mu \epsilon_\mu (x^\mu)^2 = l^2 - R^2$.

For the 5th component, since $Y^5 = R^2/l$, an expansion similar to Eq. (3.7) of the main text gives:

$$(dY^5)^2 = \frac{R^4}{l^6}\left(\sum_\nu \epsilon_\nu x_\nu dx^\nu\right)^2 = \frac{R^4}{l^6}\sum_{\nu,\rho=1}^4\,\epsilon_\nu \epsilon_\rho\, x_\nu x_\rho\, dx^\nu dx^\rho \tag{C.5}$$

Combining the coefficients of $\epsilon_\nu \epsilon_\rho\, x_\nu x_\rho$ from (C.4) and (C.5):

$$-\frac{2R^2}{l^4} + \frac{R^2 r_\epsilon^2}{l^6} + \frac{R^4}{l^6} = -\frac{2R^2}{l^4} + \frac{R^2(r_\epsilon^2 + R^2)}{l^6} = -\frac{2R^2}{l^4} + \frac{R^2 l^2}{l^6} = -\frac{R^2}{l^4}$$

This uses $r_\epsilon^2 + R^2 = l^2$ (the same cancellation mechanism as in Section 3.2 of the main text). Therefore,

$$dS^2_{\rm emb} = \frac{R^2}{l^2}\epsilon_\nu \delta_{\nu\rho}\, dx^\nu dx^\rho - \frac{R^2}{l^4}\epsilon_\nu \epsilon_\rho\, x_\nu x_\rho\, dx^\nu dx^\rho$$

and (C.3) is established. $\square$

**Remark C.1** (Condition for the Cancellation Mechanism to Work)  
The key to the proof is that $r_\epsilon^2 + R^2 = l^2$, i.e., the definitions of $l^2$ and $r_\epsilon^2$ use the **same $\epsilon_\mu$**. If the metric of the background space (C.1) is changed with $\epsilon_\mu$ but the metric of the tangent plane alone is changed, this equality breaks down and the cancellation fails (this is the reason why "a simple sign substitution is insufficient," and the simultaneous modification of the background space metric is essential).

**Remark C.2** (Verification of the Hypersurface)  
We verify the hypersurface on which the image of the mapping lies:

$$\sum_{\mu=1}^4 \epsilon_\mu (Y^\mu)^2 + (Y^5)^2 = \frac{R^2}{l^2}\,r_\epsilon^2 + \frac{R^4}{l^2} = \frac{R^2(r_\epsilon^2 + R^2)}{l^2} = \frac{R^2 \cdot l^2}{l^2} = R^2 \tag{C.6}$$

That is, the image of the mapping always lies on a hypersurface of radius $R$. The change of hypersurface depending on the sign choice:

$$\epsilon_\mu = (+1,+1,+1,+1): \quad (Y^1)^2+(Y^2)^2+(Y^3)^2+(Y^4)^2+(Y^5)^2 = R^2 \quad \text{($S^4(R)$ sphere, Euclidean)}$$

$$\epsilon_\mu = (+1,+1,+1,-1): \quad (Y^1)^2+(Y^2)^2+(Y^3)^2-(Y^4)^2+(Y^5)^2 = R^2 \quad \text{(pseudo-sphere, Lorentzian)}$$

The latter is a pseudo-sphere in $\mathbb{R}^{4,1}$, corresponding to a de Sitter hypersurface. Guo et al. established this projection as the "Beltrami-de Sitter model" in prior works [8, 9, 10] from a Lorentzian starting point, and the special case $\epsilon_4 = -1$ of the unified formula (C.3) in this appendix coincides with that model.

---

*This paper aims to organize geometric facts, and no physical interpretation is given. The sources of all formulas are indicated in the references.*
