# Geometric Origin of Lorentzian Spacetime via Gnomonic Projection: A Unified Description through the Discretization Constant $C$

**Author:** Noriaki Kihara  
**Affiliation:** WF System Co., Ltd. (Osaka University, Faculty of Engineering Science, B.Eng.)  
**ORCID:** 0009-0004-6753-4020  
**Date:** April 2026  
**Type:** Research Note (Geometric Investigation)  
**DOI:** https://doi.org/10.5281/zenodo.19434932  
**Previous Paper:** Geometric Formulation of 4-Dimensional Spacetime via Gnomonic Projection  
DOI: https://doi.org/10.5281/zenodo.19427780

---

## Abstract

In the previous paper [1], we mapped 4-dimensional Euclidean space onto the hypersphere $S^4(R)$ via gnomonic projection and derived the geometric identity for the Einstein tensor $G_{\mu\nu} + \Lambda g_{\mu\nu} = 0$ ($\Lambda = 3/R^2$). However, the following questions remained open:

1. Why does time carry a metric signature different from that of space, $(-,+,+,+)$?
2. Can the mapping be defined on a discrete coordinate system?
3. How many free parameters are required?

This paper provides geometric answers to these questions. First, motivated by the fact that the $R$-direction of the gnomonic mapping is in principle unobservable to internal observers, we axiomatically introduce the Lorentzian signature (Axiom 3.1). Second, we point out that the universal fact that physical metrics take quadratic form is a number-theoretic requirement for discretization (§6.1), and show that the mapping can be defined on discrete coordinates through the introduction of a single coupling constant $C$ (dimension: length). We prove that global regularity of the mapping is guaranteed under Euclidean signature, while a cosmological horizon naturally emerges under Lorentzian signature (Theorem 6.1). Third, we verify the consistency of physical quantities derived from $C$ (§8).

---

## 1. Introduction

### 1.1 Results of the Previous Paper

The previous paper [1] considered the central projection (gnomonic projection) from $n$-dimensional Euclidean space onto the $(n+1)$-dimensional hypersphere $S^n(R)$:

$$\Phi: \Pi_R \to S^n(R), \quad Y^\mu = \frac{Rx^\mu}{l}, \quad Y^{n+1} = \frac{R^2}{l}, \quad l = \sqrt{R^2 + \sum_\mu (x^\mu)^2} \tag{1.1}$$

For the case $n = 4$, the pullback metric was derived:

$$g_{\mu\nu} = \frac{R^2}{l^2}\left(\delta_{\mu\nu} - \frac{x_\mu x_\nu}{l^2}\right) \tag{1.2}$$

and it was shown that $S^4(R)$, as a constant-curvature space, satisfies:

$$G_{\mu\nu} + \Lambda g_{\mu\nu} = 0, \qquad \Lambda = \frac{3}{R^2} \tag{1.3}$$

Furthermore, in Appendix D, signature parameters $\epsilon_\mu \in \{+1, -1\}$ were introduced, giving the unified Euclidean/Lorentzian description:

$$g_{\mu\nu} = \frac{R^2}{l^2}\left(\epsilon_\mu \delta_{\mu\nu} - \frac{\epsilon_\mu \epsilon_\nu x_\mu x_\nu}{l^2}\right), \qquad l^2 = R^2 + \sum_\mu \epsilon_\mu (x^\mu)^2 \tag{1.4}$$

### 1.2 Open Questions

The previous paper presented $\epsilon_4 = -1$ (Lorentzian signature) as an *option*, but did not explain why it is physically realized. This paper provides a geometric motivation for this signature choice, systematically develops the consequences of *axiomatically adopting* it, points out that the universality of quadratic metrics is a necessary condition for discretization, and clarifies the motivation for introducing discrete coordinates.

### 1.3 Organization

- §2: Review of the Einstein tensor in general dimensions
- §3: Geometric motivation and axiomatic setting for the Lorentzian signature
- §4: Geometric origin of the Lorentz group
- §5: Cognitive limits of internal observers
- §6: Discrete coordinates and the coupling constant $C$
- §7: Singularity avoidance and cosmological horizons
- §8: Consistency check: physical quantities from $C$
- §9: Geometric bijectivity and its information-theoretic implications
- §10: Outlook for doctoral research

---

## 2. Einstein Tensor in General Dimensions

We organize and extend the results of the previous paper [1] to general dimensions.

**Proposition 2.1** (Einstein tensor in general dimensions [2, 3])  
On the $n$-dimensional sphere $S^n(R)$:

$$R_{\mu\nu\rho\sigma} = \frac{1}{R^2}(g_{\mu\rho}g_{\nu\sigma} - g_{\mu\sigma}g_{\nu\rho}), \qquad R_{\mu\nu} = \frac{n-1}{R^2}g_{\mu\nu} \tag{2.1}$$

$$R_{\text{scalar}} = \frac{n(n-1)}{R^2}, \qquad G_{\mu\nu} = -\frac{(n-1)(n-2)}{2R^2}g_{\mu\nu} \tag{2.2}$$

$$\Lambda_n = \frac{(n-1)(n-2)}{2R^2} \tag{2.3}$$

$\blacksquare$

| Dimension $n$ | $\Lambda_n$ | Status of Einstein tensor |
|---------------|-------------|---------------------------|
| 1             | 0           | Trivial solution (no contradiction) |
| 2             | 0           | Identically zero (trivial) |
| 3             | $1/R^2$     | Non-trivial (minimum dimension) |
| **4**         | **$3/R^2$** | **Main result of previous paper** |

**Remark 2.1:** GR produces no contradiction in any dimension. In low dimensions it merely yields trivial solutions; "GR holds" $\neq$ "non-trivial physics exists." Non-trivial GR requires $n \geq 3$, and $n = 4$ is physically distinguished.

---

## 3. Geometric Motivation and Axiomatic Setting for the Lorentzian Signature

### 3.1 Structural Singularity of the $R$-Direction

We analyze the 5-dimensional image of the gnomonic mapping (1.1):

$$\underbrace{Y^1, Y^2, Y^3, Y^4}_{\text{4-dimensional coordinates}} = \frac{R x^\mu}{l}, \qquad \underbrace{Y^5}_{\text{$R$-direction}} = \frac{R^2}{l}$$

| Property | $Y^1, \ldots, Y^4$ | $Y^5$ |
|----------|---------------------|-------|
| Behavior on sphere | "Curves" along $S^4(R)$ | Normal direction (flat) |
| Limit as $r \to \infty$ | $Y^\mu \to R\hat{x}^\mu$ (saturates) | $Y^5 \to 0$ (vanishes) |
| Value at $r = 0$ | $Y^\mu = 0$ | $Y^5 = R$ (maximum) |

**Proposition 3.0** (Asymmetry between $R$-direction and 4D coordinates)  
$Y^5$ is qualitatively different from $Y^1, \ldots, Y^4$: the 4-dimensional coordinates are intrinsic coordinates of $S^4(R)$, while $Y^5$ is the normal direction of the embedding space and is in principle unmeasurable by any measurement apparatus available to an observer confined to the sphere. $\blacksquare$

### 3.2 Imaginary Rotation of the Unobservable Dimension

**Axiom 3.1** (Axiomatic setting of the Lorentzian signature)

We adopt the following two axioms:

> **Axiom A (Cognitive constraint):** The intrinsic metric $g_{\mu\nu}$ of $S^n(R)$ (equation (1.2)) does not contain the embedding coordinate $Y^5$ (Gauss's Theorema Egregium). Therefore, an observer who uses only the intrinsic metric has no means to directly observe the $Y^5$-direction.

> **Axiom B (Imaginary rotation of the unobservable dimension):** The $R$-direction, which is intrinsically unobservable by Axiom A, is represented as an effective imaginary axis.

Axiom A is a mathematical consequence following directly from Gauss's Theorema Egregium and admits no room for assumption. Axiom B is a physical assumption that adopts the specific interpretation "unobservable dimensions are represented by imaginary axes."

The consequences under Axioms A and B are as follows:

(i) Of the four coordinates, three $x^1, x^2, x^3$ are identified as spatial coordinates and the remaining one $x^4 = ct$ as the time coordinate (this identification is an external physical input).

(ii) By Axiom B, the effective contribution of the $R$-direction is represented as an imaginary axis:

$$w_{\text{eff}} = i \cdot \frac{1}{R} \tag{3.1}$$

(iii) In the unified description of the previous paper [1] Appendix D (equation (1.4)), the structure of the gnomonic mapping is $Y^5 = R^2/l$, $l^2 = R^2 + \sum_\mu \epsilon_\mu (x^\mu)^2$. Under Axiom B, taking $R \to iR$ yields $R^2 \to -R^2$ in the definition of $l^2$. However, since $Y^5$ is unobservable (Axiom A), the overall sign of $l^2$ must be redefined in a physically meaningful way. Keeping $R^2$ positive and redefining $l^2$ accordingly forces the sign of the $x^4 = ct$ term to flip as compensation, giving $\epsilon_4 = +1 \to -1$. This provides physical motivation for the choice of $\epsilon_4$ shown in Remark D.3 of the previous paper's Appendix D. The result is:

$$w_{\text{eff}}^2 = -\frac{1}{R^2} \tag{3.2}$$

$$l^2 = R^2 + x^2 + y^2 + z^2 - c^2t^2 \tag{3.3}$$

**The Lorentzian signature $(-,+,+,+)$ is obtained under Axioms A and B.**

**Remark 3.1** (Status of this claim)  
Axiom A is mathematically rigorous (a direct consequence of Gauss's Theorema Egregium), but Axiom B is a physical assumption. The identification "unobservable → imaginary" could in principle be replaced by other interpretations (e.g., the unobservable dimension is simply eliminated in the effective theory, treated as a gauge degree of freedom, or compactified à la Kaluza-Klein). This paper explores the consequences of adopting Axiom B. The justification for Axiom B lies in the fact that its consequences (Lorentzian signature, Lorentz group, cosmological horizon) are consistent with observation.

### 3.3 Relationship to Conventional Wick Rotation

| Conventional Wick rotation | This paper's formulation |
|----------------------------|--------------------------|
| $t \to it$ introduced as mathematical procedure | Unobservability of $R$-direction motivates $w = i/R$ (Axiom B) |
| Lorentzian signature is empirical fact | Lorentzian signature is a consequence of Axioms A, B |
| No explanation for why sign reverses | Unobservable orthogonal axis → imaginary axis → negative sign |

**This paper reverses the causal direction:** Wick rotation is not a "convenient technique"; instead, the Lorentzian world is the "cognitive reality" for internal observers, while the Euclidean structure is the "underlying geometric entity."

---

## 4. Geometric Origin of the Lorentz Group

### 4.1 Chain of Symmetry Groups

The isometry group of $S^4(R)$ depends on the signature choice.

**Proposition 4.1** (Chain of symmetry groups [3, 6])

$$SO(5) \xrightarrow[\text{Axiom 3.1}]{\epsilon_4: +1 \to -1} SO(4,1) \xrightarrow[R \to \infty]{\text{flat limit}} SO(3,1) \tag{4.1}$$

- $SO(5)$: Isometry group of Euclidean $S^4(R)$ (10 generators)
- $SO(4,1)$: de Sitter group (isometry group of Lorentzian $dS_4$)
- $SO(3,1)$: Lorentz group (isometry group of Minkowski spacetime)

$\blacksquare$

**Significance:** Lorentz transformations are not "derived from the principle of the constancy of the speed of light"; rather, they are geometrically derived as the "flat limit of the de Sitter symmetry of the gnomonic mapping." The constancy of the speed of light is a consequence of this geometric structure; the more fundamental principle is the symmetry of the gnomonic mapping.

---

## 5. Cognitive Limits of Internal Observers

### 5.1 Upper Bound on the Sum of Interior Angles

**Theorem 5.1** (Interior angle upper bound)  
For a geodesic triangle (2-dimensional cross-section) on $S^n(R)$ with interior angles $\alpha, \beta, \gamma$:

$$\alpha + \beta + \gamma = \pi + \frac{A}{R^2} \tag{5.1}$$

where $A$ is the area of the geodesic triangle. At maximum area (hemisphere) $A_{\max} = 2\pi R^2$:

$$(\alpha + \beta + \gamma)_{\max} = 3\pi = 540° \tag{5.2}$$

**This upper bound of 540° is a universal geometric constant independent of $R$.** $\blacksquare$

**Proof:** By the Gauss–Bonnet theorem [2], for a geodesic triangle on a surface of Gauss curvature $K = 1/R^2$: $\alpha + \beta + \gamma - \pi = KA = A/R^2$. The maximum value of $A$ is the hemisphere area $2\pi R^2$, so $(\alpha + \beta + \gamma)_{\max} = \pi + 2\pi = 3\pi$. This value cancels $R$. $\square$

**Remark 5.1** (On the measurability of $R$)  
The fact that the upper bound 540° does not depend on $R$ is different from $R$ itself being unmeasurable. An internal observer can *determine* the value of $R$ by measuring the angular excess $\Delta = A/R^2$ for a triangle of known area $A$. More precisely, $R$ is in principle measurable through the intrinsic curvature $K = 1/R^2$ via the Riemann tensor. What Theorem 5.1 shows is only the geometric fact that the *upper bound* does not depend on $R$.

### 5.2 Structure of Conformal Scaling

**Proposition 5.2** (Properties of the conformal factor)  
The gnomonic mapping contains the conformal factor $R/l$. Since $l = \sqrt{R^2 + r^2}$ depends on position, this is not uniform scaling but a **position-dependent conformal transformation**:

$$x^\mu \to \frac{R}{l}x^\mu \tag{5.3}$$

At the origin ($r = 0$): $R/l = 1$ (no transformation); as $r \to \infty$: $R/l \to 0$ (compression). $\blacksquare$

---

## 6. Discrete Coordinates and the Coupling Constant $C$

### 6.1 Motivation: Number-Theoretic Relationship between Quadratic Metrics and Discretization

The previous paper [1] worked with continuous coordinates $x^\mu \in \mathbb{R}$. This section discretizes the coordinates and shows that the mapping can be defined on an integer lattice.

Before introducing discretization, we note the following number-theoretic fact.

**Observation 6.0** (Number-theoretic motivation for quadratic forms)  
When one attempts to define a metric on the integer lattice $\mathbb{Z}^n$ that possesses a multiplicative norm ($N(zw) = N(z)N(w)$) and takes integer values, the norm $N(a+bi) = a^2 + b^2$ of the Gaussian integer ring $\mathbb{Z}[i]$ provides the minimal structure. More generally, by Hurwitz's theorem (1898) [15], composition algebras over the reals with multiplicative norms exist only in dimensions $1, 2, 4, 8$ (reals, complex numbers, quaternions, octonions). The norm $a^2 + b^2$ of $\mathbb{Z}[i]$ is the minimal-dimensional such structure, and the quadratic form of physical metrics is a natural choice for defining a multiplicatively closed metric on discrete lattices. $\blacksquare$

**Remark 6.0** (Circumstantial evidence for discreteness)  
The fact that physical metrics take quadratic form ($ds^2 = g_{\mu\nu}dx^\mu dx^\nu$, $l^2 = R^2 + \sum \epsilon_\mu(x^\mu)^2$, etc.) constitutes circumstantial evidence that our spacetime is discrete. Quadratic metrics are certainly consistent with continuous spacetime (no contradiction), but they are *required* for discretization. The logical implication is one-directional:

$$\text{Discreteness} \Longrightarrow \text{Quadratic metrics are required}$$
$$\text{Quadratic metrics} \not\Longrightarrow \text{Discreteness (but serve as circumstantial evidence)}$$

### 6.2 Introduction of the Coupling Constant $C$

**Definition 6.0** (Coupling constant $C$)  
We introduce a constant $C > 0$ with dimension of length and define the curvature parameter $w = C/R$. The following discretization examines the consequences of the **additional assumption** that $w$ is restricted to non-negative integers:

$$w = \frac{C}{R} \in \mathbb{Z}_{\geq 0}, \qquad R = \frac{C}{w} \tag{6.1}$$

The dimensionless coordinates $\tilde{x}^\mu = x^\mu / C$ are restricted to $\tilde{x}^\mu \in \mathbb{Z}$.

### 6.3 Discrete Gnomonic Mapping

**Definition 6.1** (Discrete gnomonic mapping)  
In dimensionless quantities:

$$\tilde{L}^2 = 1 + w^2 \tilde{\sigma}^2, \qquad \tilde{\sigma}^2 = \sum_\mu \epsilon_\mu (\tilde{x}^\mu)^2 \tag{6.2}$$

$$\tilde{Y}^\mu = \frac{\tilde{x}^\mu}{\tilde{L}}, \qquad \tilde{Y}^5 = \frac{1}{w\tilde{L}} \tag{6.3}$$

where $w, \tilde{x}^\mu \in \mathbb{Z}$ and $\epsilon_\mu = (+1,+1,+1,-1)$.

### 6.4 Domain and Regularity of the Mapping

**Theorem 6.1** (Signature-dependent regularity)  
For the discrete gnomonic mapping (6.2)–(6.3):

**(i) Euclidean signature** ($\epsilon_\mu = (+1,+1,+1,+1)$):

$$\tilde{\sigma}^2 = \tilde{x}^2 + \tilde{y}^2 + \tilde{z}^2 + \tilde{t}^2 \geq 0 \quad \Longrightarrow \quad \tilde{L}^2 = 1 + w^2\tilde{\sigma}^2 \geq 1 > 0 \tag{6.4}$$

The denominator $\tilde{L}$ **never vanishes at any lattice point**, and the mapping is globally regular.

**(ii) Lorentzian signature** ($\epsilon_\mu = (+1,+1,+1,-1)$):

$$\tilde{\sigma}^2 = \tilde{x}^2 + \tilde{y}^2 + \tilde{z}^2 - \tilde{t}^2 \tag{6.5}$$

In the timelike region ($\tilde{t}^2 > \tilde{x}^2 + \tilde{y}^2 + \tilde{z}^2$), $\tilde{\sigma}^2 < 0$ is possible, so $\tilde{L}^2 = 1 + w^2\tilde{\sigma}^2 < 1$. In particular, when $w^2|\tilde{\sigma}^2| > 1$, $\tilde{L}^2 < 0$ and the mapping becomes undefined.

**The mapping is defined only in the region $\tilde{L}^2 > 0$.** $\blacksquare$

**Proof:** (i) Euclidean case: $\tilde{\sigma}^2 \geq 0$ (each term is a sum of squares of integers), so $w^2\tilde{\sigma}^2 \geq 0$, hence $\tilde{L}^2 \geq 1$. (ii) Lorentzian case: $\tilde{\sigma}^2 = \tilde{x}^2 + \tilde{y}^2 + \tilde{z}^2 - \tilde{t}^2$ can be negative. The boundary $\tilde{L}^2 = 0$ gives:

$$1 + w^2(\tilde{x}^2 + \tilde{y}^2 + \tilde{z}^2 - \tilde{t}^2) = 0 \quad \Longleftrightarrow \quad c^2t^2 = R^2 + x^2 + y^2 + z^2 \tag{6.6}$$

This corresponds to the **cosmological horizon** in de Sitter spacetime [6]. $\square$

**Remark 6.1** (Geometric meaning of the horizon)  
The horizon in equation (6.6) emerges naturally as the boundary of the domain of the gnomonic mapping without additional assumptions. This is consistent with the known structure identified by Guo et al. [6] in Beltrami–de Sitter coordinates. Beltrami coordinates cover only the **static patch** ($l^2 > 0$) of de Sitter spacetime, and the existence of the horizon is an essential property of the gnomonic mapping (= Beltrami coordinates).

**Remark 6.2** (Lattice solutions of $\tilde{L}^2 = 0$ and Diophantine equations)  
In the discrete Lorentzian version, the existence of lattice points on the horizon ($\tilde{L}^2 = 0$) is non-trivial. For $w = 1$ and isotropic lattice points $\tilde{x} = \tilde{y} = \tilde{z} = s$, the condition reduces to the Pell equation $\tilde{t}^2 - 3s^2 = 1$, which has integer solutions (minimal solution $(s, \tilde{t}) = (1, 2)$). This demonstrates that horizons also appear on discrete lattices.

Furthermore, the author's previous work [14] completely determined the existence conditions for solutions of the isotropic geodesic equation $l^2 + t^2 = Nx^2$ ($l, t, x \in \mathbb{Z}_{>0}$) on the discrete Minkowski lattice $\mathbb{Z}^{1+N}$. By Fermat's two-square theorem and the method of infinite descent:

> **Theorem ([14] Theorem 3.1):** A necessary and sufficient condition for $l^2 + t^2 = Nx^2$ to have non-trivial positive integer solutions is that in the prime factorization of $N$, every prime $p \equiv 3 \pmod{4}$ appears with even exponent.

This is directly related to the gnomonic discretization of the present paper. The framework of this paper can be extended to hyperspheres $S^n(R)$ of arbitrary dimension (Proposition 2.1), and assuming a single time dimension, the existence of isotropic discrete geodesics for $N$ spatial dimensions is subject to number-theoretic constraints depending on the prime structure of $N$. The complete classification from $N = 1$ to $N = 11$ is shown below:

| Spatial dim. $N$ | Spacetime dim. $D = N+1$ | Prime factors of $N$ | Isotropic solutions | Notes |
|---|---|---|---|---|
| 1 | 2 | $1$ | Exist ✅ | Pythagorean triples |
| 2 | 3 | $2$ | Exist ✅ | |
| **3** | **4** | $3$ ($3 \equiv 3 \pmod{4}$, odd power) | **None ❌** | **Our universe (spatial only)** |
| 4 | 5 | $2^2$ | Exist ✅ | Kaluza-Klein extension |
| 5 | 6 | $5$ ($5 \equiv 1 \pmod{4}$) | Exist ✅ | |
| 6 | 7 | $2 \times 3$ ($3$: odd power) | None ❌ | |
| 7 | 8 | $7$ ($7 \equiv 3 \pmod{4}$, odd power) | None ❌ | |
| 8 | 9 | $2^3$ | Exist ✅ | |
| 9 | 10 | $3^2$ (even power) | Exist ✅ | Superstring theory dimension |
| 10 | 11 | $2 \times 5$ | Exist ✅ | M-theory dimension |
| 11 | 12 | $11$ ($11 \equiv 3 \pmod{4}$, odd power) | None ❌ | |

Additionally, $N = 25$ ($D = 26$, $5^2$) admits isotropic solutions, coinciding with the dimension of bosonic string theory.

The fact that $N = 3$ exhibits an "arithmetic anomaly" and that $N = 4, 9, 10, 25$ coincide with specific dimensions known in theoretical physics constitutes an important number-theoretic constraint in the arbitrary-dimension extension of this framework. Which additional dimensions correspond to which physical forces lies beyond the scope of this paper and is left for future work.

**Corollary 6.2** (Formal consistency with the continuous limit)  
As $C \to 0$ (zero lattice spacing), $w = C/R \to 0$ (for fixed $R$) and $\tilde{x}^\mu = x^\mu/C \to \infty$ (for fixed $x^\mu$). Then $\tilde{L}^2 = 1 + (C/R)^2 \sum \epsilon_\mu (x^\mu/C)^2 = 1 + \sum \epsilon_\mu (x^\mu)^2/R^2 = l^2/R^2$. Therefore $\tilde{Y}^\mu = x^\mu/(l/R) = Rx^\mu/l = Y^\mu$, which **formally** agrees with the continuous version (1.1). However, the mathematical rigor of convergence from $w \in \mathbb{Z}$ to the continuum (recovery of the continuum through densification of the lattice) is left for future work. $\square$

### 6.5 Physical Identification of $C$

**Hypothesis 6.3** (Identification of $C$ with the Planck length)

$$C \sim \ell_P = \sqrt{\frac{\hbar G}{c^3}} \approx 1.616 \times 10^{-35}\,\text{m} \tag{6.5}$$

Under this hypothesis, the following quantities are quantized:

| Physical quantity | Quantization form | Value for $C = \ell_P$ |
|-------------------|-------------------|------------------------|
| Minimum length    | $\Delta x = C$    | $\ell_P \approx 1.6 \times 10^{-35}$ m |
| Minimum time      | $\Delta t = C/c$  | $t_P \approx 5.4 \times 10^{-44}$ s |
| Cosmological constant | $\Lambda = 3w^2/C^2$ | Discrete values only |

---

## 7. Singularity Avoidance and Cosmological Horizons

### 7.1 Regularity of the Euclidean Version (Global)

**Theorem 7.1** (Restatement of Proposition C.1 from [1])  
Under Euclidean signature ($\epsilon_\mu = (+1,+1,+1,+1)$), for fixed $R > 0$, the gnomonic mapping $\Phi: \mathbb{R}^4 \to S^4(R)$ is $C^\infty$-regular on all of $\mathbb{R}^4$. Since $l = \sqrt{R^2 + r^2} \geq R > 0$, $l$ is always positive. $\blacksquare$

### 7.2 Regularity of the Lorentzian Version (Within the Static Patch)

**Theorem 7.2** (Regularity of the Lorentzian version)  
Under Lorentzian signature ($\epsilon_\mu = (+1,+1,+1,-1)$), the gnomonic mapping is $C^\infty$-regular in the region satisfying $l^2 = R^2 + x^2 + y^2 + z^2 - c^2t^2 > 0$ (static patch). The boundary $l^2 = 0$, i.e., $c^2t^2 = R^2 + x^2 + y^2 + z^2$, defines the **cosmological horizon** (Theorem 6.1 (ii)). $\blacksquare$

**Remark 7.1** (Regularity of the discrete version)  
In the discrete Euclidean version, by Theorem 6.1 (i), $\tilde{L} \geq 1 > 0$ holds at all lattice points. In the discrete Lorentzian version, regularity holds only at lattice points satisfying $\tilde{L}^2 > 0$.

### 7.3 Singularity Avoidance as $R \to 0$

**Proposition 7.3**  
In the Euclidean version, as $R \to 0$ (with $r > 0$ fixed):

$$Y^5 = \frac{R^2}{\sqrt{R^2 + r^2}} \to 0, \qquad Y^\mu = \frac{Rx^\mu}{\sqrt{R^2 + r^2}} \to 0 \tag{7.1}$$

The image degenerates to the origin but **does not diverge**. That is, even in the limit where the curvature radius approaches zero, the mapping remains mathematically defined. $\blacksquare$

### 7.4 Information Preservation

**Theorem 7.4** (Bijectivity and information preservation)

**(i) Euclidean version:** For fixed $R > 0$, the gnomonic mapping $\Phi: \mathbb{R}^4 \to S^4(R)$ is a $C^\infty$-bijection onto the open hemisphere. The inverse mapping $\Phi^{-1}$ exists, and information is preserved in principle.

**(ii) Lorentzian version:** $\Phi$ is a bijection from the region $l^2 > 0$ (static patch) to the corresponding region of de Sitter spacetime. The region beyond the horizon ($l^2 = 0$) must be covered by a different coordinate patch.

$\blacksquare$

---

## 8. Consistency Checks

### 8.1 Consistency of Physical Quantities under $C = \ell_P$

We *assume* Hypothesis 6.3 ($C = \ell_P$) and verify that the physical quantities derived from it do not contradict known values. This is a consistency check, not a verification of the theory; independent determination of $C$ is left as a future task in §10.

**Check 1: Present value of the cosmological constant $\Lambda$**

From the observed value $\Lambda_{\text{obs}} \approx 1.1 \times 10^{-52}\,\text{m}^{-2}$ [7] and $\Lambda = 3/R^2$:

$$R_{\text{current}} = \sqrt{\frac{3}{\Lambda_{\text{obs}}}} \approx 1.65 \times 10^{26}\,\text{m} \approx 5.3\,\text{Gpc} \tag{8.1}$$

This is consistent in order of magnitude with the size of the observable universe $R_H \approx 4.4 \times 10^{26}$ m.

**Check 2: Discrete curvature parameter $w$**

$$w = \frac{C}{R_{\text{current}}} = \frac{1.6 \times 10^{-35}}{1.65 \times 10^{26}} \approx 10^{-61} \tag{8.2}$$

$w \ll 1$, and the present universe corresponds to a discrete lattice point extremely close to $w = 0$ (Minkowski flat spacetime).

**Check 3: Integer condition on $w$**

$w \in \mathbb{Z}$ is required, but $w \approx 10^{-61} \ll 1$. In this model's discrete $w$, no intermediate state exists between $w = 0$ ($\Lambda = 0$, perfectly flat) and $w = 1$ ($\Lambda = 3/C^2 \sim 10^{122} \Lambda_{\text{obs}}$, Planck density). The observed small $\Lambda > 0$ cannot be reproduced by discrete steps of $w$; this is the discrete version of the cosmological constant problem. Its resolution requires modification of the discreteness of $w$ (e.g., allowing $w$ to take rational values, or including contributions from additional dimensions), which is included in the future tasks of §10.

### 8.2 Consistency with Planck Scale

| Physical quantity | Value from $C = \ell_P$ | Known value | Consistency |
|---|---|---|---|
| Planck length | $C = 1.616 \times 10^{-35}$ m | $\ell_P = 1.616 \times 10^{-35}$ m | ✅ (by definition) |
| Planck time | $C/c = 5.39 \times 10^{-44}$ s | $t_P = 5.39 \times 10^{-44}$ s | ✅ |
| $R_{\text{current}}$ | $\sqrt{3/\Lambda} = 1.65 \times 10^{26}$ m | $R_H \sim 4.4 \times 10^{26}$ m | ✅ (same order) |
| Initial $w$ | $w_{\text{initial}} \sim C/C = 1$ | — | Consistent (maximum curvature) |

### 8.3 Limitations and Caveats

**Remark 8.1** (On avoiding circular reasoning)  
The fact that Planck quantities emerge when one assumes $C = \ell_P$ is a direct consequence of the assumption and is not an independent verification. The only non-trivial aspect of this section is that $R_{\text{current}}$ turns out to be of the **same order as the observable radius of the universe $R_H$, independently** (Check 1). True verification of $C$ requires determining $C$ from multiple independent observables and showing their consistency. This is included in the doctoral-level tasks in §10.

---

## 9. Geometric Bijectivity and Information-Theoretic Implications

### 9.1 Absence of Geometric Singularities in the Euclidean Version

**Proposition 9.1**  
Under Euclidean signature ($\epsilon_\mu = (+1,+1,+1,+1)$), since the gnomonic mapping $\Phi$ is a $C^\infty$-bijection (Theorems 7.1, 7.4 (i)), for arbitrarily small $R$:

1. No geometric singularity arises
2. Information is preserved in principle
3. The inverse mapping $\Phi^{-1}$ exists, and information is recoverable

$\blacksquare$

### 9.2 Restrictions in the Lorentzian Version

Under Lorentzian signature, information preservation is limited to the static patch ($l^2 > 0$) (Theorem 7.4 (ii)). The region beyond the horizon ($l^2 = 0$) is not covered by Beltrami coordinates and requires a different coordinate system (e.g., de Sitter global coordinates). How this restriction affects black hole physics is left for future work.

### 9.3 Information Preservation in the Discrete Version

In the discrete Euclidean coordinate system (§6), with $C > 0$ and by Theorem 6.1 (i), $\tilde{L} \geq 1 > 0$ is guaranteed at all lattice points. Bijectivity on a finite lattice trivially guarantees information preservation.

### 9.4 Context

The results of this section are propositions of pure geometry and are independent of physical discussions such as quantum gravity or Hawking radiation. Application to the physical black hole information paradox requires formulation of quantum theory within the gnomonic framework, which is included in the future tasks of §10.

---

## 10. Outlook for Doctoral Research

This paper has completed the geometric origin of Lorentzian spacetime through the gnomonic projection and coupling constant $C$, within the framework of pure 4-dimensional spacetime excluding electromagnetic forces. The following tasks exceed the scope of this paper and are delegated to doctoral-level research.

### 10.1 Extension to Arbitrary Dimensions and Discrete Number-Theoretic Constraints

The gnomonic framework of this paper can be defined for hyperspheres $S^n(R)$ of arbitrary $n$ (Proposition 2.1). Assuming a single time dimension, the discrete number-theoretic constraints shown in Remark 6.2 (Fermat's two-square theorem [14]) strongly constrain the spatial dimensions $N$ that admit isotropic discrete geodesics. In particular, that $N = 3$ (4D spacetime) is arithmetically singular and that $N = 4, 9, 10, 25$ coincide with specific dimensions in theoretical physics is an important subject for future research.

Notably, while $N = 3$ (4D spacetime) admits no isotropic discrete geodesics, adding one dimension to $N = 4$ (5D spacetime) yields the first solutions. This number-theoretic constraint suggests that a **fifth dimension is necessarily required** for discrete spacetime to maintain isotropy. In Kaluza-Klein theory [8, 9], a compact fifth dimension $\phi$ is identified with the charge (Coulomb) dimension, and Maxwell's equations are derived from the $U(1)$ gauge symmetry thereon. Combined with the number-theoretic constraint of this paper, a logical chain of "requirement of discrete isotropy → necessity of the 5th dimension → Kaluza-Klein interpretation → emergence of electromagnetic force" may be established. The rigorous development of this connection is the subject of a forthcoming paper. Preliminary outlook was given in Appendix E of the previous paper [1].

### 10.2 Geometric Origin of the Fine-Structure Constant $\alpha$

Connection between the embedding $S^4 \hookrightarrow S^5$ and Wyler's formula (1971) [10]. Preliminary numerical calculations reproduce $\alpha_{\text{Wyler}} = 1/137.036$ to $0.0001\%$ precision. It has been verified that the isometry group $SO(5,1)$ of the gnomonic model transforms to $SO(4,2) \cong SU(2,2)/\mathbb{Z}_2$ under Wick rotation, but the complete derivation of Wyler's formula requires representation theory of Bergman kernels and Hermitian symmetric spaces, which remains an open task.

### 10.3 Strong and Weak Forces

- **Strong force:** Geometric description of the confinement potential $V \sim \kappa r$. Point sources are insufficient; string-type sources (flux tubes) are needed, and the introduction of "string-type sources" in the discrete lattice extension of §6 is a candidate.
- **Weak force:** Flavor transitions ($u \to d + W^+$) are state transitions rather than coordinate changes, and the present model's coordinate system lacks the degrees of freedom to describe them. An essentially new extension is needed.

### 10.4 Geometric Foundations of Quantum Theory

The discrete coordinate system (§6) is structurally similar to lattice field theory and loop quantum gravity. Formulation of connections to path integrals, operator formalism, and geometric interpretation of wave functions.

### 10.5 Derivation of the Particle Mass Spectrum

Verification of whether known elementary particle masses can be reproduced as a discrete spectrum from the constraint $w_n = C/R_n \in \mathbb{Z}$.

### 10.6 Quantitative Consistency with Cosmological Observables

Derivation of the dynamical equation for $R(t)$, detailed correspondence with the Friedmann equation, and quantitative explanation of inflation and dark energy.

---

## Summary of Results

| Item | Result | Eq. No. | Status |
|------|--------|---------|--------|
| Axiom for Lorentzian signature | Unobservability of $R$-direction → Axioms A, B → $(-,+,+,+)$ | (3.1)–(3.3) | **Axiom 3.1** |
| Motivation for quadratic metrics | Multiplicative norms on discrete lattices naturally require quadratic forms (Hurwitz [15]) | §6.1 | **Observation 6.0** |
| Chain of symmetry groups | $SO(5) \to SO(4,1) \to SO(3,1)$ | (4.1) | **Proposition 4.1** |
| Interior angle upper bound | $(\alpha+\beta+\gamma)_{\max} = 540°$ (**upper bound** is $R$-independent) | (5.2) | **Theorem 5.1** |
| Coupling constant | $w = C/R \in \mathbb{Z}_{\geq 0}$, $C$: sole free parameter | (6.1) | **Definition 6.0** |
| Regularity (Euclidean) | $\tilde{L}^2 \geq 1 > 0$ (globally regular) | (6.4) | **Theorem 6.1 (i)** |
| Horizon (Lorentzian) | Cosmological horizon emerges at $l^2 = 0$ | (6.6) | **Theorem 6.1 (ii)** |
| Singularity avoidance | $\Phi$ is a $C^\infty$ bijection in Euclidean version | (7.1) | **Theorems 7.1, 7.4** |
| Consistency of $R_{\text{current}}$ | $\sqrt{3/\Lambda} \sim R_H$ (same order) | (8.1) | **Consistency check** |

---

## References

[1] Kihara, N. (2026). Geometric Formulation of 4-Dimensional Spacetime via Gnomonic Projection. DOI: 10.5281/zenodo.19427780.

[2] do Carmo, M. P. (1992). *Riemannian Geometry*. Birkhäuser, Chapter 4.

[3] Wald, R. M. (1984). *General Relativity*. University of Chicago Press, Appendix D.

[4] Hartle, J. B., & Hawking, S. W. (1983). Wave function of the Universe. *Physical Review D*, 28(12), 2960.

[5] Gibbons, G. W., & Hawking, S. W. (1977). Action integrals and partition functions in quantum gravity. *Physical Review D*, 15(10), 2752.

[6] Guo, H.-Y., Huang, C.-G., Xu, Z., Zhou, B. (2004). On Beltrami model of de Sitter spacetime. *Modern Physics Letters A*, 19, 1701–1710. arXiv:hep-th/0405137.

[7] Planck Collaboration (2020). Planck 2018 results. VI. Cosmological parameters. *Astronomy & Astrophysics*, 641, A6.

[8] Kaluza, T. (1921). Zum Unitätsproblem der Physik. *Sitzungsber. Preuss. Akad. Wiss.*, 966–972.

[9] Klein, O. (1926). Quantentheorie und fünfdimensionale Relativitätstheorie. *Zeitschrift für Physik*, 37, 895–906.

[10] Wyler, A. (1971). L'espace symétrique du groupe des équations de Maxwell. *Comptes Rendus Acad. Sci. Paris*, 272, 186–188.

[11] Misner, C. W., Thorne, K. S., Wheeler, J. A. (1973). *Gravitation*. W. H. Freeman, §13.5.

[12] Spivak, M. (1999). *A Comprehensive Introduction to Differential Geometry*, Vol. 2, 3rd ed. Publish or Perish, Chapter 4.

[13] Snyder, J. P. (1987). *Map Projections: A Working Manual*. USGS Professional Paper 1395, pp. 164–168.

[14] Kihara, N. (2026). Solution of Positive Integer Solutions in Discrete Minkowski Intervals Using Gaussian Integers. GitHub: ai-chat-logs-open.

[15] Hurwitz, A. (1898). Über die Composition der quadratischen Formen von beliebig vielen Variablen. *Nachrichten von der Gesellschaft der Wissenschaften zu Göttingen*, 309–316.

---

*The discussion in this paper is confined to the domain of pure geometry. Physical claims, experimental predictions, and identification with existing theories are left for future research.*  
*Last updated: 2026-04-06*
