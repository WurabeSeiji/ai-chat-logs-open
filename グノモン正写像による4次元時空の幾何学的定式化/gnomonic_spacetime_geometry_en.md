# Geometric Formulation of 4-Dimensional Spacetime via Gnomonic Projection

**Author:** Noriaki Kihara  
**Affiliation:** WF System Co., Ltd. (Osaka University, School of Engineering Science, Alumni)  
**Date:** April 2026  
**Type:** Research Note (Geometric Analysis)

---

## Abstract

This note examines, from a purely geometric viewpoint, the deformation of the induced metric tensor that arises when 4-dimensional Euclidean space is mapped onto the surface of a 5-dimensional hypersphere of curvature radius $R$ via the gnomonic (central) projection. In particular, we present step-by-step derivations of: (i) the exact induced metric tensor $g_{\mu\nu}$ (including the contribution of the fifth component), and (ii) the geometric relationship between the Einstein tensor and the constant $\Lambda = 3/R^2$. Physical interpretation is deferred to future work; no cosmological assumptions are made in this note.

---

## 1. Introduction

The gnomonic projection is defined as the central projection from the center of a sphere onto a tangent plane, and its most important property — that great circles (geodesics) on the sphere map to straight lines on the plane — is well established [1].

This note extends this projection from 4 to 5 dimensions and investigates the metric structure that arises in the resulting coordinate system. The problem is purely geometric; physical interpretations (gravity, spacetime, relativity, etc.) lie outside the scope of this note.

---

## 2. Preliminaries: Definition and Basic Properties of the Gnomonic Projection

### 2.1 Gnomonic Projection in General Dimension

**Definition 2.1** ($n$-dimensional gnomonic projection) [1, 2]  
In $n$-dimensional Euclidean space, consider a tangent hyperplane $\Pi_R$ at distance $R$ from the center of a sphere (with coordinates $(x^1, \ldots, x^n) \in \mathbb{R}^n$, the origin at the point of tangency). Setting $l = \sqrt{R^2 + \sum_{i=1}^n (x^i)^2}$, we call

$$\Phi: \Pi_R \to S^n(R) \subset \mathbb{R}^{n+1}, \quad (x^1, \ldots, x^n) \mapsto \left(\frac{R x^1}{l}, \ldots, \frac{R x^n}{l},\ \frac{R^2}{l}\right) \tag{2.1}$$

the **gnomonic projection**, where $S^n(R)$ is the $n$-dimensional sphere of radius $R$ embedded in $\mathbb{R}^{n+1}$:

$$S^n(R) = \left\{(Y^1, \ldots, Y^{n+1}) \in \mathbb{R}^{n+1} \;\middle|\; \sum_{A=1}^{n+1} (Y^A)^2 = R^2 \right\}.$$

**Proposition 2.1** (Geodesic preservation) [1]  
Under the gnomonic projection, any straight line on the tangent plane maps to a great circle (geodesic) on the sphere.

---

## 3. Derivation of the Metric Tensor

### 3.1 Setup

We set $n=4$ and consider the gnomonic projection from the tangent hyperplane $\Pi_R$ (coordinates $(x^1, x^2, x^3, x^4)$) to the 5-dimensional hypersphere $S^4(R) \subset \mathbb{R}^5$ of radius $R$. Since the tangent plane is Euclidean, $x_\mu = x^\mu$ ($\mu = 1,2,3,4$). By Definition 2.1, the 5-dimensional target coordinates are

$$(Y^1, Y^2, Y^3, Y^4, Y^5) = \left(\frac{Rx^1}{l},\ \frac{Rx^2}{l},\ \frac{Rx^3}{l},\ \frac{Rx^4}{l},\ \frac{R^2}{l}\right) \tag{3.1}$$

$$l = \sqrt{R^2 + \sum_{\mu=1}^4 (x^\mu)^2},\quad r^2 = \sum_{\mu=1}^4(x^\mu)^2,\quad l^2 = R^2 + r^2. \tag{3.2}$$

We abbreviate $(Y^1, \ldots, Y^4)$ as $(X^1, \ldots, X^4)$.

### 3.2 Derivation of the Pullback Metric

The induced metric on $S^4(R) \subset \mathbb{R}^5$ is determined by the standard inner product on $\mathbb{R}^5$ via [6]:

$$dS^2 = \sum_{A=1}^{5} (dY^A)^2. \tag{3.3}$$

**We compute all five components, including the fifth.**

**Components 1–4** $(Y^\mu = R x^\mu / l)$:

**Lemma 3.1**
$$\frac{\partial Y^\mu}{\partial x^\nu} = \frac{R}{l}\left(\delta^\mu_\nu - \frac{x^\mu x^\nu}{l^2}\right) \tag{3.4}$$

**Proof**
$$\frac{\partial}{\partial x^\nu}\!\left(\frac{Rx^\mu}{l}\right) = \frac{R\delta^\mu_\nu}{l} + Rx^\mu\!\left(-\frac{x^\nu}{l^3}\right) = \frac{R}{l}\!\left(\delta^\mu_\nu - \frac{x^\mu x^\nu}{l^2}\right) \quad \square$$

Therefore,

$$\sum_{\mu=1}^4(dY^\mu)^2 = \frac{R^2}{l^2}\sum_{\nu,\rho=1}^{4}\left(\delta_{\nu\rho} - \frac{2x_\nu x_\rho}{l^2} + \frac{r^2 x_\nu x_\rho}{l^4}\right)dx^\nu dx^\rho. \tag{3.5}$$

**Fifth component** $(Y^5 = R^2/l)$:

$$\frac{\partial Y^5}{\partial x^\nu} = -\frac{R^2 x^\nu}{l^3} \tag{3.6}$$

$$\left(dY^5\right)^2 = \frac{R^4}{l^6}\left(\sum_\nu x^\nu dx^\nu\right)^2 = \frac{R^4}{l^6}\,x_\nu x_\rho\,dx^\nu dx^\rho \tag{3.7}$$

**Summing all components:**

Collecting the coefficient of $x_\nu x_\rho$ from (3.5) and (3.7):

$$-\frac{2R^2}{l^4} + \frac{R^2 r^2}{l^6} + \frac{R^4}{l^6} = -\frac{2R^2}{l^4} + \frac{R^2(r^2+R^2)}{l^6} = -\frac{2R^2}{l^4} + \frac{R^2 l^2}{l^6} = -\frac{R^2}{l^4}.$$

Therefore,

$$dS^2 = \frac{R^2}{l^2}\delta_{\nu\rho}\,dx^\nu dx^\rho - \frac{R^2}{l^4}\,x_\nu x_\rho\,dx^\nu dx^\rho$$

$$\boxed{g_{\mu\nu} = \frac{R^2}{l^2}\left(\delta_{\mu\nu} - \frac{x_\mu x_\nu}{l^2}\right)} \tag{3.8}$$

**Remark 3.1** (Indispensability of the fifth component)  
If the fifth component $(dY^5)^2$ (eq. (3.7)) is omitted, the coefficient of $x_\nu x_\rho$ remains $-2R^2/l^4 + R^2 r^2/l^6$ and eq. (3.8) cannot be obtained. The contribution $+R^4/l^6 = +R^2/l^4 - R^2 r^2/l^6$ from the fifth component exactly compensates this deficit.

---

## 4. Christoffel Symbols and Riemann Curvature

### 4.1 Riemann Curvature Tensor

**Theorem 4.1** (Curvature tensor of a space of constant curvature) [3, 4, 5]  
The Riemann curvature tensor of the $n$-sphere $S^n(R)$ is given by

$$R_{\mu\nu\rho\sigma} = \frac{1}{R^2}(g_{\mu\rho}g_{\nu\sigma} - g_{\mu\sigma}g_{\nu\rho}). \tag{4.1}$$

**Remark 4.1**  
Eq. (4.1) is the standard result for spaces of constant curvature [4, 5]. The direct derivation in gnomonic coordinates (computing $\Gamma^\lambda_{\mu\nu}$ from the Christoffel formula) is planned for Appendix A; in this note we cite references [3, 4, 5].

---

## 5. Ricci Tensor and Einstein Tensor

### 5.1 Ricci Tensor

**Proposition 5.1**  
The Ricci tensor of the $n$-sphere $S^n(R)$ is

$$R_{\mu\nu} = \frac{n-1}{R^2}g_{\mu\nu}. \tag{5.1}$$

**Proof**  
Substituting (4.1) into the definition $R_{\mu\nu} = g^{\rho\sigma}R_{\mu\rho\nu\sigma}$:

$$g^{\rho\sigma}R_{\mu\rho\nu\sigma} = \frac{1}{R^2}g^{\rho\sigma}(g_{\mu\nu}g_{\rho\sigma} - g_{\mu\sigma}g_{\nu\rho}).$$

Computing each term:

- **First term:** $g^{\rho\sigma}g_{\mu\nu}g_{\rho\sigma} = g_{\mu\nu}\cdot\underbrace{g^{\rho\sigma}g_{\rho\sigma}}_{=\,\delta^\rho{}_\rho\,=\,n} = n\,g_{\mu\nu}$
- **Second term:** $g^{\rho\sigma}g_{\mu\sigma}g_{\nu\rho} = \underbrace{g^{\rho\sigma}g_{\mu\sigma}}_{=\,\delta^\rho_\mu}g_{\nu\rho} = g_{\nu\mu} = g_{\mu\nu}$

Therefore,

$$R_{\mu\nu} = \frac{1}{R^2}(n\,g_{\mu\nu} - g_{\mu\nu}) = \frac{n-1}{R^2}g_{\mu\nu}. \qquad \square$$

### 5.2 Ricci Scalar

$$R_{\rm scalar} = g^{\mu\nu}R_{\mu\nu} = \frac{n-1}{R^2}\underbrace{g^{\mu\nu}g_{\mu\nu}}_{=\,n} = \frac{n(n-1)}{R^2}. \tag{5.2}$$

For $n=4$: $R_{\rm scalar} = 12/R^2$.

### 5.3 Einstein Tensor

From the definition of the Einstein tensor [4]:

$$G_{\mu\nu} = R_{\mu\nu} - \frac{1}{2}g_{\mu\nu}R_{\rm scalar}. \tag{5.3}$$

For $n=4$:

$$G_{\mu\nu} = \frac{3}{R^2}g_{\mu\nu} - \frac{1}{2}\cdot\frac{12}{R^2}\,g_{\mu\nu} = \left(\frac{3}{R^2} - \frac{6}{R^2}\right)g_{\mu\nu} = -\frac{3}{R^2}g_{\mu\nu}. \tag{5.4}$$

---

## 6. Main Theorem

**Theorem 6.1** (Main result)  
When the tangent hyperplane $\Pi_R$ of 4-dimensional Euclidean space is mapped to the hypersphere $S^4(R)$ of curvature radius $R$ via the gnomonic projection, the induced metric tensor $g_{\mu\nu}$ (eq. (3.8)) satisfies

$$G_{\mu\nu} + \Lambda\, g_{\mu\nu} = 0 \tag{6.1}$$

with

$$\boxed{\Lambda = \frac{3}{R^2}.} \tag{6.2}$$

**Proof**  
From (5.4), $G_{\mu\nu} = -\dfrac{3}{R^2}g_{\mu\nu}$; setting $\Lambda = \dfrac{3}{R^2}$ immediately yields (6.1). $\square$

**Remark 6.0** (Notation)  
In this theorem, the $n=4$ result is stated with $\Lambda$ (no subscript). The generalization to arbitrary dimension is given in Remark 6.2 as $\Lambda_n$.

**Remark 6.1**  
Eq. (6.1) is a **geometric identity for $S^4(R)$ with Euclidean signature** and is not to be interpreted as a physical law. It is structurally distinct from the de Sitter spacetime with Lorentzian signature; the two are related by analytic continuation (Wick rotation $t \to it$, see [5] Appendix D).

Note that both the Euclidean and Lorentzian versions share the same value $\Lambda = 3/R^2$, consistent with the framework of Gibbons–Hawking Euclidean quantum gravity [8]. In this sense, the sign difference between Euclidean and Lorentzian is a **technical issue** bridged by the known mathematical procedure of Wick rotation, indicating that the present formulation is not in conflict with the vacuum solution of general relativity. Physical implications are deferred to future work.

**Remark 6.2** (Extension to general dimension)  
This theorem is not restricted to $n=4$. For a general $n$-sphere $S^n(R)$, from (5.2) and (5.3):

$$G_{\mu\nu} = \frac{n-1}{R^2}g_{\mu\nu} - \frac{1}{2}\cdot\frac{n(n-1)}{R^2}g_{\mu\nu} = -\frac{(n-1)(n-2)}{2R^2}g_{\mu\nu},$$

so the value of $\Lambda$ satisfying $G_{\mu\nu} + \Lambda g_{\mu\nu} = 0$ is

$$\Lambda_n = \frac{(n-1)(n-2)}{2R^2}$$

(for $n=4$: $\Lambda_4 = 3/R^2$). The case $n=4$, i.e., $3+1$ dimensions, is physically distinguished by its correspondence with the de Sitter vacuum solution; this interpretation is left for future work.

---

## 7. Topics for Future Investigation

The following items are planned for investigation based on the geometric results of this note (all are outside the scope of the current work).

- **Geometric interpretation of Wick rotation:** How to describe the operation $t \to it$ within the framework of the gnomonic projection — the central open problem bridging the present Euclidean formulation to the Lorentzian version.
- **Lorentzian gnomonic projection:** Derivation of the induced metric for $l^2 = R^2 + x^2 + y^2 + z^2 - t^2$ and verification of the sign change.
- **Geodesics under symmetric conditions:** Geodesic equations for $X^2 = X^3 = X^4 = b$, $X^1 = t$.
- **Geometric meaning of $b = t$:** Interpretation of the equal-distance curves satisfying $r^2 = 3t^2$.
- **Introduction of matter:** Extension to the case $T_{\mu\nu} \neq 0$.
- **Euclidean/Lorentzian correspondence:** Detailed examination of the relationship between Theorem 6.1 and the de Sitter vacuum solution (see [8]), and rigorous proof of the correspondence between the unified formula (D.3) of Appendix D and the Beltrami–de Sitter metric of Guo et al. [9, 10, 11].

---

## 8. Summary of Results

| Item | Result | Eq. | Source |
|------|--------|-----|--------|
| Map definition | $Y^\mu = Rx^\mu/l$, $Y^5 = R^2/l$ | (3.1) | [1, 2] |
| Induced metric | $g_{\mu\nu} = \dfrac{R^2}{l^2}\!\left(\delta_{\mu\nu} - \dfrac{x_\mu x_\nu}{l^2}\right)$ | (3.8) | This work |
| Riemann curvature | $R_{\mu\nu\rho\sigma} = \dfrac{1}{R^2}(g_{\mu\rho}g_{\nu\sigma}-g_{\mu\sigma}g_{\nu\rho})$ | (4.1) | [3,4,5] |
| Ricci tensor | $R_{\mu\nu} = \dfrac{3}{R^2}g_{\mu\nu}$ | (5.1) | This work |
| Ricci scalar | $R_{\rm scalar} = \dfrac{12}{R^2}$ | (5.2) | This work |
| Einstein tensor | $G_{\mu\nu} = -\dfrac{3}{R^2}g_{\mu\nu}$ | (5.4) | This work |
| **Main theorem** | $G_{\mu\nu} + \Lambda g_{\mu\nu} = 0$, $\Lambda = \dfrac{3}{R^2}$ | (6.2) | This work |

---

## References

[1] Snyder, J. P., *Map Projections: A Working Manual*, U.S. Geological Survey Professional Paper 1395, U.S. Government Printing Office, Washington D.C., 1987, pp. 164–168.

[2] Coxeter, H. S. M., *Introduction to Geometry*, 2nd ed., Wiley, New York, 1969, Chapter 15.

[3] do Carmo, M. P., *Riemannian Geometry*, Birkhäuser, Boston, 1992, Chapter 4.

[4] Misner, C. W., Thorne, K. S., Wheeler, J. A., *Gravitation*, W. H. Freeman, San Francisco, 1973, §13.5.

[5] Wald, R. M., *General Relativity*, University of Chicago Press, Chicago, 1984, Appendix D.

[6] Spivak, M., *A Comprehensive Introduction to Differential Geometry*, Vol. 2, 3rd ed., Publish or Perish, Houston, 1999, Chapter 4.

[7] Lee, J. M., *Introduction to Riemannian Manifolds*, 2nd ed., Springer, New York, 2018, Example 3.14.

[8] Gibbons, G. W. and Hawking, S. W., "Action integrals and partition functions in quantum gravity," *Physical Review D*, **15**, 2752–2756, 1977.  
　　(Foundations of Euclidean quantum gravity; relation between Euclidean de Sitter spacetime and $\Lambda = 3/R^2$.)

[9] Guo, H.-Y., Huang, C.-G., Xu, Z., Zhou, B., "On Beltrami model of de Sitter spacetime," *Modern Physics Letters A*, **19**, 1701–1710, 2004. arXiv:hep-th/0405137.  
　　(Construction of de Sitter Beltrami coordinates via Lorentzian gnomonic projection; prior work on the $\epsilon_4=-1$ case of Appendix D.)

[10] Guo, H.-Y., "The Beltrami Model of De Sitter Space," arXiv:hep-th/0607017, 2006.  
　　(Direct construction of the de Sitter metric via gnomonic projection and de Sitter-invariant special relativity.)

[11] Guo, H.-Y., Huang, C.-G., Xu, Z., Zhou, B., "Beltrami Model of Riemann-Sphere and de Sitter Spacetime," arXiv:gr-qc/0703078, 2007.  
　　(Comparison of Beltrami models for $S^n(R)$ (Euclidean) and $dS^n$ (Lorentzian); connection to the unified perspective of this appendix.)

---

## Appendix A: Direct Computation of Christoffel Symbols (Planned)

A complete derivation of $\Gamma^\lambda_{\mu\nu}$ from the metric (3.8) via the Christoffel formula, leading directly to the Riemann curvature tensor of Theorem 4.1, is planned for a future report. At present we cite references [3, 4].

---

## Appendix B: Flat Limit $R \to \infty$

As $R \to \infty$, we have $l \to R$, so

$$g_{\mu\nu} \to \delta_{\mu\nu}, \quad \Lambda = \frac{3}{R^2} \to 0.$$

The metric converges to the flat Euclidean metric and the curvature vanishes. This corresponds to the gnomonic projection degenerating to the identity map (for details of the limiting procedure, see [3] Chapter 4).

---

## Appendix C: Regularity of the Projection

**Proposition C.1** (Regularity of the gnomonic projection)  
Let $\Phi: \mathbb{R}^4 \to S^4(R)$ be the gnomonic projection of Definition 2.1. Setting $r = \sqrt{\sum_\mu (x^\mu)^2}$, the following hold for

$$Y^5 = \frac{R^2}{\sqrt{R^2 + r^2}}, \quad Y^\mu = \frac{Rx^\mu}{\sqrt{R^2+r^2}}:$$

1. **$R \to 0$ ($r > 0$ fixed):** $Y^5 \to 0$, $Y^\mu \to 0$. The image degenerates to the origin but does not diverge.

2. **$R \to \infty$ ($r$ fixed):** $Y^5 \to \infty$. The fifth component diverges; $Y^\mu \to x^\mu$ (finite).

3. **Fixed $R > 0$:** Since $l = \sqrt{R^2+r^2} \geq R > 0$, the denominator is always positive and $\Phi$ is a $C^\infty$ map on all of $\mathbb{R}^4$.

**Proof**  
All claims follow by direct computation of $l = \sqrt{R^2+r^2}$. (1) $R\to 0$: $Y^5 = R^2/l \leq R^2/r \to 0$. (2) $R\to\infty$: $Y^5 = R(1+r^2/R^2)^{-1/2} \to \infty$. (3) For $R>0$, $l \geq R > 0$, so each component is $C^\infty$ in $r$. $\square$

**Remark C.1**  
The divergence $Y^5 \to \infty$ as $R \to \infty$ (Proposition C.1-(2)) is consistent with the flat limit of Appendix B: as the sphere becomes "infinitely large," the difference from the tangent plane vanishes while the projection in the embedding direction retreats to infinity.

---

## Appendix D: Essential Observation on the Squared Metric — Unified Description of Euclidean and Lorentzian Geometries

Examining the formulation of this note, one finds that the transition from Euclidean to Lorentzian geometry — within the geometric description of the gnomonic projection — can be described as a choice of the real parameter $\epsilon_\mu$. We record this as a purely mathematical observation.

Note that the Lorentzian gnomonic projection ($\epsilon_4 = -1$) and its correspondence to the Beltrami coordinate system of de Sitter spacetime has been established as prior work by H.-Y. Guo et al., starting from the Lorentzian side [9, 10, 11]. The originality of this appendix lies in the **explicit derivation of the unified formula (D.3)**, which encompasses both the Euclidean case ($S^4(R)$, [8]) and the Lorentzian case (de Sitter–Beltrami, [9, 10, 11]) as special cases of a single parametrization via $\epsilon_\mu$.

**Observation D.1** (Quadratic nature of physical metrics)  
All metrics appearing in the gnomonic projection involve only $t^2$, never $t$ itself:

$$l^2 = R^2 + x^2 + y^2 + z^2 \pm t^2, \quad ds^2 = \pm dt^2 + dx^2 + dy^2 + dz^2.$$

The quantities $t^2, x^2, y^2, z^2$ are always non-negative real numbers. The sign difference is merely a **choice of coefficient $\epsilon_\mu \in \{+1,-1\}$** assigned to these squared quantities. The coordinate $t$ itself is an auxiliary quantity; what is fundamental to the physical metric is always $t^2$.

**Proposition D.1** (Unified description of Euclidean/Lorentzian)  
Choose a sign system $(\epsilon_1, \epsilon_2, \epsilon_3, \epsilon_4) \in \{+1,-1\}^4$.

Define the **embedding space metric** — not the standard metric on $\mathbb{R}^5$ — as

$$dS^2_{\rm emb} = \sum_{\mu=1}^{4} \epsilon_\mu (dY^\mu)^2 + (dY^5)^2 \tag{D.1}$$

(the coefficient of the fifth component $Y^5$ is fixed at $+1$).

Define the metric on the tangent hyperplane $\Pi_R$ as

$$ds^2_{\rm tan} = \sum_{\mu=1}^{4} \epsilon_\mu (dx^\mu)^2 \tag{D.2}$$

and set $l^2 = R^2 + \sum_{\mu=1}^4 \epsilon_\mu (x^\mu)^2$. The gnomonic map $\Phi$ is defined as in (2.1):

$$Y^\mu = \frac{R x^\mu}{l} \quad (\mu=1,2,3,4), \quad Y^5 = \frac{R^2}{l}.$$

The induced metric pulled back by $\Phi$ is then

$$g_{\mu\nu} = \frac{R^2}{l^2}\left(\epsilon_\mu \delta_{\mu\nu} - \frac{\epsilon_\mu \epsilon_\nu x_\mu x_\nu}{l^2}\right) \tag{D.3}$$

(Einstein summation convention is not used for repeated indices).

**Proof**  
We express $dS^2_{\rm emb}$ under (D.1) in tangent-plane coordinates.

Computing $dY^\mu$ as in Lemma 3.1:

$$\sum_{\mu=1}^4 \epsilon_\mu (dY^\mu)^2 = \frac{R^2}{l^2} \sum_{\nu,\rho=1}^4 \left(\epsilon_\nu \delta_{\nu\rho} - \frac{2\epsilon_\nu \epsilon_\rho x_\nu x_\rho}{l^2} + \frac{\epsilon_\nu \epsilon_\rho r_\epsilon^2 x_\nu x_\rho}{l^4}\right) dx^\nu dx^\rho \tag{D.4}$$

where $r_\epsilon^2 = \sum_\mu \epsilon_\mu (x^\mu)^2 = l^2 - R^2$.

From $Y^5 = R^2/l$:

$$(dY^5)^2 = \frac{R^4}{l^6}\left(\sum_\nu \epsilon_\nu x_\nu dx^\nu\right)^2 = \frac{R^4}{l^6}\,\epsilon_\nu \epsilon_\rho\, x_\nu x_\rho\, dx^\nu dx^\rho. \tag{D.5}$$

Summing the coefficients of $\epsilon_\nu \epsilon_\rho x_\nu x_\rho$ from (D.4) and (D.5):

$$-\frac{2R^2}{l^4} + \frac{R^2 r_\epsilon^2}{l^6} + \frac{R^4}{l^6} = -\frac{2R^2}{l^4} + \frac{R^2(r_\epsilon^2 + R^2)}{l^6} = -\frac{2R^2}{l^4} + \frac{R^2 l^2}{l^6} = -\frac{R^2}{l^4},$$

using $r_\epsilon^2 + R^2 = l^2$ (the same cancellation mechanism as in Section 3.2). Therefore,

$$dS^2_{\rm emb} = \frac{R^2}{l^2}\epsilon_\nu \delta_{\nu\rho}\, dx^\nu dx^\rho - \frac{R^2}{l^4}\epsilon_\nu \epsilon_\rho\, x_\nu x_\rho\, dx^\nu dx^\rho,$$

which establishes (D.3). $\square$

**Remark D.2** (Condition for the cancellation mechanism)  
The key to the proof is $r_\epsilon^2 + R^2 = l^2$ — that the definitions of $l^2$ and $r_\epsilon^2$ use the **same $\epsilon_\mu$**. If only the tangent-plane metric is changed while the embedding metric (D.1) is left unmodified, this identity breaks and the cancellation fails. (This is why a naive sign substitution is insufficient; simultaneous modification of the embedding metric is essential.)

**Remark D.3** (Geometric correspondence with Wick rotation)  
Within the geometric description of the gnomonic projection, the Wick rotation $t \to it$ can be described as

$$\epsilon_4 : +1 \to -1:$$

$$\text{Euclidean}: \quad \epsilon_\mu = (+1,+1,+1,+1) \quad \Rightarrow \quad l^2 = R^2 + x^2 + y^2 + z^2 + t^2$$

$$\text{Lorentzian}: \quad \epsilon_\mu = (+1,+1,+1,-1) \quad \Rightarrow \quad l^2 = R^2 + x^2 + y^2 + z^2 - t^2.$$

Both are described as special cases of a **single geometric framework** sharing the embedding metric (D.1) and the map definition. Note that the necessity of Wick rotation for the convergence of path integrals (in quantum gravity, etc.) falls outside the scope of this note; here we state only the correspondence at the level of purely geometric description.

**Remark D.4** (Scope of this appendix and verification of the hypersurface)  
Proposition D.1 is a mathematical consequence of defining the embedding metric as (D.1). We verify the hypersurface on which the image of the map lies:

$$\sum_{\mu=1}^4 \epsilon_\mu (Y^\mu)^2 + (Y^5)^2 = \frac{R^2}{l^2}\,r_\epsilon^2 + \frac{R^4}{l^2} = \frac{R^2(r_\epsilon^2 + R^2)}{l^2} = R^2. \tag{D.6}$$

The image always lies on a hypersurface of radius $R$. The hypersurface depends on the sign choice:

$$\epsilon_\mu = (+1,+1,+1,+1): \quad (Y^1)^2+(Y^2)^2+(Y^3)^2+(Y^4)^2+(Y^5)^2 = R^2 \quad \text{($S^4(R)$, Euclidean)}$$

$$\epsilon_\mu = (+1,+1,+1,-1): \quad (Y^1)^2+(Y^2)^2+(Y^3)^2-(Y^4)^2+(Y^5)^2 = R^2 \quad \text{(pseudo-sphere, Lorentzian)}$$

The latter is a pseudo-sphere in $\mathbb{R}^{4,1}$, structurally corresponding to the de Sitter hypersurface. Guo et al. established this projection from the Lorentzian starting point as the "Beltrami–de Sitter model" in prior work [9, 10, 11]; it is expected that the $\epsilon_4 = -1$ special case of the unified formula (D.3) coincides with that model. A rigorous proof of the correspondence starting from (D.3) is left for future work.

---

*This note aims to organize geometric facts and offers no physical interpretation. All sources for the equations are cited in the references.*
