# Geometric Origin of the Charge Dimension in Discrete Gnomonic Projection: 5-Dimensional Spacetime and Charge Quantization

**Author:** Noriaki Kihara  
**Affiliation:** WF System Co., Ltd. (Osaka University, Faculty of Engineering Science, B.Eng.)  
**ORCID:** 0009-0004-6753-4020  
**Date:** April 2026  
**Type:** Research Note (Geometric Investigation)  
**Previous Papers:**
- [1] Geometric Formulation of 4-Dimensional Spacetime via Gnomonic Projection. DOI: 10.5281/zenodo.19427780.
- [2] Geometric Origin of Lorentzian Spacetime via Gnomonic Projection: A Unified Description through the Discretization Constant $C$.

---

## Abstract

In the previous paper [2], we showed that the isotropic geodesic equation $l^2 + t^2 = Nx^2$ on the discrete Minkowski lattice $\mathbb{Z}^{1+N}$ has no positive integer solutions for $N = 3$ (4-dimensional spacetime) — an arithmetic anomaly — and that solutions first exist for $N = 4$ (5-dimensional spacetime). This paper identifies the fifth dimension required by this number-theoretic constraint with the charge (Coulomb) dimension and extends the gnomonic projection to $S^5(R)$.

No new coupling constants are introduced. The sole coupling constant $C$ (dimension: length) is retained, and a conversion coefficient $Q = C/e_0$ ($e_0$: minimum charge quantum) linking the charge dimension to the length dimension is defined. This makes the coordinate on the discrete lattice $\tilde{\phi} = q/e_0 \in \mathbb{Z}$ an integer charge number, and charge quantization is naturally derived as a consequence of the discrete lattice without assuming compactification. The theory does not specify the value of $e_0$, treating it as a quantity to be determined by experiment. In current particle physics, $e_0 = e/3$ (quark charge), and all known particle charges, including the fractional charges of quarks, are expressed as integers $\tilde{\phi}$.

Furthermore, the $1/r^2$ Coulomb law is derived from the gnomonic mapping in the charge direction, and combined with the Lorentz invariance from [2], Maxwell's equations in their entirety are geometrically derived. Additionally, from the difference in observability of the curvature directions, the difference in polarization structure between electromagnetic waves (spin-1) and gravitational waves (spin-2), as well as the absence of the equivalence principle for electromagnetic forces, emerge naturally.

---

## 1. Introduction: Confirmation of Discretization Premises

### 1.1 Premise: The Discrete Gnomonic Mapping Has Been Verified

This paper takes the **discrete gnomonic projection** established in the previous paper [2] as its starting point. The following results have already been verified and are used as premises in this paper:

> **Established results of [2]:**
>
> 1. The introduction of the coupling constant $C$ allows continuous coordinates $x^\mu$ to be discretized into dimensionless integer coordinates $\tilde{x}^\mu = x^\mu/C \in \mathbb{Z}$ (Definition 6.0).
> 2. Under the curvature parameter $w = C/R \in \mathbb{Z}_{\geq 0}$, the discrete gnomonic mapping $\tilde{Y}^\mu = \tilde{x}^\mu / \tilde{L}$ is defined (Definition 6.1).
> 3. **All metrics are in quadratic form** $\tilde{\sigma}^2 = \sum_\mu \epsilon_\mu (\tilde{x}^\mu)^2$, and $\tilde{L}^2 = 1 + w^2 \tilde{\sigma}^2$ is defined as a sum of integer squares.
> 4. In the Euclidean version, $\tilde{L}^2 \geq 1 > 0$ is **guaranteed at all lattice points** (Theorem 6.1 (i)).
> 5. The continuous limit $C \to 0$ formally agrees with the continuous version (Corollary 6.2).

All discussions from this paper onward proceed under this discretization framework.

### 1.2 Fundamental Reason for Quadratic Metrics: Handling Negative Coordinate Values

Working under the premise of discretization gives rise to the following fundamental question: why are physical metrics universally quadratic?

Observation 6.0 of [2] provided number-theoretic motivation from the norm of Gaussian integers and Hurwitz's theorem [7]. However, in this paper we present a more fundamental reason:

**On the discrete lattice $\mathbb{Z}^n$, each axis coordinate takes not only positive but also negative integers:**

$$\tilde{x}^\mu \in \mathbb{Z} = \{\ldots, -2, -1, 0, +1, +2, \ldots\}$$

A metric is a measure of distance and should measure magnitude (distance), not the sign (direction) of coordinates. Therefore, a metric must return the same value for $\tilde{x}^\mu$ and $-\tilde{x}^\mu$. The minimal structure that guarantees this is the quadratic form $(\tilde{x}^\mu)^2$:

$$(-n)^2 = n^2 \qquad \forall n \in \mathbb{Z} \tag{1.1}$$

**This property is the necessary condition for defining a metric on a discrete lattice, and is the most fundamental reason for quadratic metrics.** The logical chain of implications:

$$\text{Discrete lattice (integer coordinates)} \Longrightarrow \text{Negative values exist on coordinates}$$
$$\Longrightarrow \text{Metric must be sign-invariant} \Longrightarrow \text{Quadratic form } (\tilde{x}^\mu)^2 \text{ is required}$$
$$\Longrightarrow \text{Gnomonic mapping is defined for all integers}$$

**Remark 1.1** (Physical meaning of negative values on each axis)

| Axis | Meaning of $+n$ | Meaning of $-n$ | Metric treatment | Physically distinguishable? |
|------|------------------|------------------|------------------|----------------------------|
| $x, y, z$ | right/front/up | left/back/down | $(\pm n)^2 = n^2$ | No (isotropy) |
| $ct$ (time) | future | past | $(\pm n)^2 = n^2$ | Conventional (definition of arrow) |
| $\phi = Qq$ (charge) | positive charge | negative charge | $(\pm n)^2 = n^2$ | Yes (physical) |

For all axes, the metric $(\tilde{x}^\mu)^2$ handles negative coordinate values without issue. Even when a physical distinction between positive and negative exists (charge), that distinction is preserved in **the coordinate value $\tilde{x}^\mu$ itself**, not in the metric $(\tilde{x}^\mu)^2$. This is consistent with "positive and negative charges have the same gravitational mass" (equivalence principle).

### 1.3 Number-Theoretic Motivation: Why Five Dimensions Are Needed

In Remark 6.2 of [2], based on the author's earlier work [3], the following was shown:

> On the discrete Minkowski lattice $\mathbb{Z}^{1+N}$, the condition for the isotropic geodesic equation $l^2 + t^2 = Nx^2$ to have non-trivial positive integer solutions is that in the prime factorization of $N$, every prime $p \equiv 3 \pmod{4}$ appears with even exponent (Fermat's two-square theorem).

$N = 3$ (3 spatial dimensions) satisfies $3 \equiv 3 \pmod{4}$ and violates this condition. $N = 4$ (4 spatial dimensions) satisfies $4 = 2^2$ and meets the condition.

**This number-theoretic constraint shows that discrete spacetime requires at least 4 spatial dimensions ($N = 4$) to maintain isotropy.** The purpose of this paper is to identify the fourth spatial dimension with the charge (Coulomb) dimension.

### 1.4 Observability of the Charge Dimension

In Axiom 3.1 of [2], the unobservability of the $R$-direction of $S^n(R)$ (Axiom A) led to imaginary rotation (Axiom B), yielding the Lorentzian signature.

The treatment of the charge dimension $\phi$ added in this paper is fundamentally different from the $R$-direction:

| Dimension | Observability | Application of Axiom B (imaginary rotation) |
|-----------|---------------|---------------------------------------------|
| $Y^{n+1}$ ($R$-direction) | ❌ Unobservable (Gauss's Theorema Egregium) | Applied → sign reversal |
| $\phi = Qq$ (charge direction) | ✅ **Observable** (charge can be measured) | **Not applied** → sign unchanged |

**Charge is an observable physical quantity.** Therefore, the imaginary rotation of Axiom B is not applied to the charge dimension, and $\phi$ contributes to the metric as an ordinary spatial dimension. Humans can measure the physical quantity of Coulombs (charge) with measuring instruments, which means the charge dimension is not an "unobservable hidden dimension."

### 1.5 Organization

- §2: Extension of the 5D gnomonic mapping
- §3: Conversion coefficient $Q$ and charge quantization
- §4: 5D discrete gnomonic mapping
- §5: Recovery of discrete isotropy at $N = 4$
- §6: 5D Einstein tensor and cosmological constant
- §7: Consistency checks
- §8: Future tasks

---

## 2. Extension of the 5-Dimensional Gnomonic Mapping

### 2.1 Definition of the Coordinate System

We define the coordinates of 5-dimensional space as:

$$(x^1, x^2, x^3, x^4, x^5) = (x, y, z, \phi, ct) \tag{2.1}$$

where:
- $x, y, z$: spatial coordinates (dimension: length)
- $\phi = Qq$: charge coordinate (dimension: length). $Q$ is the charge-to-length conversion coefficient (dimension: m/C)
- $ct$: time coordinate (dimension: length). $c$ is the time-to-length conversion coefficient

Signature parameters:

$$\epsilon_\mu = (+1, +1, +1, +1, -1), \qquad \mu = 1, \ldots, 5 \tag{2.2}$$

Four spatial dimensions ($x, y, z, \phi$) carry positive signature; one time dimension ($ct$) carries negative signature.

### 2.2 5-Dimensional Gnomonic Mapping

Applying the general-dimension result of [1] to $n = 5$:

$$\Phi: \Pi_R \to S^5(R), \quad Y^\mu = \frac{Rx^\mu}{l}, \quad Y^6 = \frac{R^2}{l} \tag{2.3}$$

$$l^2 = R^2 + \sum_{\mu=1}^{5} \epsilon_\mu (x^\mu)^2 = R^2 + x^2 + y^2 + z^2 + Q^2 q^2 - c^2t^2 \tag{2.4}$$

Pullback metric:

$$g_{\mu\nu} = \frac{R^2}{l^2}\left(\epsilon_\mu \delta_{\mu\nu} - \frac{\epsilon_\mu \epsilon_\nu x_\mu x_\nu}{l^2}\right) \tag{2.5}$$

### 2.3 Unobservability of $Y^6$ and Lorentzian Signature

Axiom 3.1 of [2] applies directly. The intrinsic metric of $S^5(R)$ does not contain the embedding coordinate $Y^6$ (Gauss's Theorema Egregium), so $Y^6$ is unobservable. Applying imaginary rotation to the $R$-direction via Axiom B reverses the sign of $x^5 = ct$, yielding the Lorentzian signature (2.2).

**Remark 2.1** (Distinction between the charge dimension and the $R$-direction)  
The $\phi$-direction is one of the intrinsic coordinates of $S^5(R)$ and is explicitly contained in the intrinsic metric (2.5). Therefore, $\phi$ is observable, and the imaginary rotation of Axiom B is not applied. The positive signature $\epsilon_4 = +1$ of $\phi$ is a direct consequence of the observability of charge.

---

## 3. Conversion Coefficient $Q$ and Charge Quantization

### 3.1 Structural Symmetry with $c$

The speed of light $c$ converts the dimension of time (seconds) to the dimension of length (meters):

$$x^5 = c \cdot t \qquad [c] = \text{m/s} \tag{3.1}$$

Similarly, $Q$ converts the dimension of charge (coulombs) to the dimension of length (meters):

$$x^4 = \phi = Q \cdot q \qquad [Q] = \text{m/C} \tag{3.2}$$

| | Time → Length | Charge → Length |
|--|---------------|-----------------|
| Conversion coefficient | $c$ | $Q$ |
| Conversion formula | $x^5 = ct$ | $x^4 = Qq$ |
| Dimension | m/s | m/C |
| Known value | $2.998 \times 10^8$ m/s | Determined in §3.2 |

### 3.2 Determination of $Q$

The discretization requirement $\tilde{\phi} = Qq/C \in \mathbb{Z}$ demands that a **minimum charge quantum** $e_0$ exists and all charges are integer multiples thereof:

$$Q = \frac{C}{e_0}, \qquad \tilde{\phi} = \frac{q}{e_0} \in \mathbb{Z} \tag{3.3}$$

The theory does not specify the value of $e_0$; it is a quantity to be determined by experiment. In current particle physics, the smallest known charge is the quark charge $e/3$ ($e$: absolute value of the electron charge). Therefore:

$$e_0 = \frac{e}{3} \approx 5.34 \times 10^{-20} \text{ C} \tag{3.4}$$

$$Q = \frac{C}{e_0} = \frac{3C}{e} = \frac{3 \times 1.616 \times 10^{-35}}{1.602 \times 10^{-19}} \approx 3.03 \times 10^{-16} \text{ m/C} \tag{3.5}$$

### 3.3 Charge Quantization Derived from the Discrete Lattice

**Theorem 3.1** (Charge quantization)  
Under the discrete lattice requirement $\tilde{\phi} \in \mathbb{Z}$ and the conversion coefficient $Q = C/e_0$:

$$\tilde{\phi} = \frac{Qq}{C} = \frac{q}{e_0} \in \mathbb{Z} \tag{3.6}$$

Therefore, charge is quantized as integer multiples of the minimum charge quantum $e_0$:

$$q = n \cdot e_0, \qquad n \in \mathbb{Z} \tag{3.7}$$

$\blacksquare$

**Remark 3.1** (Indeterminacy of $e_0$)  
The theory requires only that "charge is discrete"; the value of the minimum charge quantum $e_0$ is not a free parameter of the theory but is experimentally determined. If a charge unit finer than $e/3$ is discovered in the future, replacing $e_0$ with that value leaves the framework unchanged. The theory requires only the "integer multiple" structure.

**Remark 3.2** (Compactification is unnecessary)  
The charge quantization (3.6)–(3.7) does not assume any compactification of the $\phi$-direction. In Kaluza-Klein theory, charge quantization is derived from the periodicity of a compact dimension [4, 5], but in this model it is derived directly from the integer coordinate $\tilde{\phi} \in \mathbb{Z}$ on the discrete lattice. The $\phi$-direction is non-compact, just like $x, y, z$.

**Remark 3.3** (Discrete coordinates $\tilde{\phi}$ for known particles)  
With $e_0 = e/3$:

| Particle | Charge $q$ | $\tilde{\phi} = q/e_0 = 3q/e$ | Metric contribution $\tilde{\phi}^2$ |
|----------|-----------|------|------|
| Proton | $+e$ | $+3$ | $9$ |
| Electron | $-e$ | $-3$ | $9$ |
| u quark | $+2e/3$ | $+2$ | $4$ |
| d quark | $-e/3$ | $-1$ | $1$ |
| Neutrino | $0$ | $0$ | $0$ |
| Positron | $+e$ | $+3$ | $9$ |

All are integers, and the metric $\tilde{\phi}^2$ is identical for positive and negative values. In particular, the electron ($\tilde{\phi} = -3$) and the positron ($\tilde{\phi} = +3$) are metrically identical ($\tilde{\phi}^2 = 9$), with the sign distinction preserved only in the coordinate value.

---

## 4. 5-Dimensional Discrete Gnomonic Mapping

### 4.1 Discretization

Extending to 5 dimensions using the coupling constant $C$ from §6 of [2]:

$$\tilde{x} = \frac{x}{C}, \quad \tilde{y} = \frac{y}{C}, \quad \tilde{z} = \frac{z}{C}, \quad \tilde{\phi} = \frac{Qq}{C} = \frac{q}{e_0}, \quad \tilde{t} = \frac{ct}{C} \tag{4.1}$$

All dimensionless coordinates are integers: $\tilde{x}, \tilde{y}, \tilde{z}, \tilde{\phi}, \tilde{t} \in \mathbb{Z}$.

### 4.2 Discrete Metric

$$\tilde{L}^2 = 1 + w^2 \tilde{\sigma}^2, \qquad \tilde{\sigma}^2 = \tilde{x}^2 + \tilde{y}^2 + \tilde{z}^2 + \tilde{\phi}^2 - \tilde{t}^2 \tag{4.2}$$

$$\tilde{Y}^\mu = \frac{\tilde{x}^\mu}{\tilde{L}}, \qquad \tilde{Y}^6 = \frac{1}{w\tilde{L}} \tag{4.3}$$

### 4.3 Regularity (5D Version)

**Theorem 4.1** (Regularity of the 5D discrete gnomonic mapping)

**(i) Euclidean signature** ($\epsilon_\mu = (+1,+1,+1,+1,+1)$):

$$\tilde{\sigma}^2 = \tilde{x}^2 + \tilde{y}^2 + \tilde{z}^2 + \tilde{\phi}^2 + \tilde{t}^2 \geq 0 \quad \Longrightarrow \quad \tilde{L}^2 \geq 1 > 0 \tag{4.4}$$

Globally regular. $\tilde{L}^2 \geq 1$ is guaranteed **for all integer coordinates (including negative values)**.

**(ii) Lorentzian signature** ($\epsilon_\mu = (+1,+1,+1,+1,-1)$):

$$\tilde{\sigma}^2 = \tilde{x}^2 + \tilde{y}^2 + \tilde{z}^2 + \tilde{\phi}^2 - \tilde{t}^2 \tag{4.5}$$

When the 4D spatial part (sum of squares) is at least $\tilde{t}^2$, $\tilde{L}^2 \geq 1$. **The addition of the charge term $\tilde{\phi}^2$ expands the region $\tilde{L}^2 > 0$ compared to the 4D version.** The horizon is:

$$c^2t^2 = R^2 + x^2 + y^2 + z^2 + Q^2q^2 \tag{4.6}$$

$\blacksquare$

**Remark 4.1** (Expansion of regularity by charge)  
The horizon in the 4D version ([2] equation (6.6)) was $c^2t^2 = R^2 + x^2 + y^2 + z^2$. In the 5D version (4.6), $Q^2q^2 \geq 0$ is added to the right-hand side. Therefore, charged particles ($q \neq 0$) have greater margin to the horizon and remain within the horizon over a broader region than in the 4D version for the same time coordinate.

**Remark 4.2** (Safety for negative charge values)  
When $\tilde{\phi} = -n$ ($n > 0$, negative charge), $\tilde{\phi}^2 = n^2 > 0$. By the quadratic metric, the contribution to $\tilde{\sigma}^2$ is always non-negative, and the value of $\tilde{L}^2$ does not depend on the sign of $\tilde{\phi}$. The gnomonic mapping is fully defined for negative charge coordinates.

---

## 5. Recovery of Discrete Isotropy at $N = 4$

### 5.1 Application of the Gaussian Integer Theorem

We apply the theorem cited in Remark 6.2 of [2] [3] to $N = 4$.

For $N = 4$ (4 spatial dimensions), $4 = 2^2$ contains no prime factors with $p \equiv 3 \pmod{4}$. Therefore, the isotropic geodesic equation

$$l^2 + t^2 = 4x^2 \tag{5.1}$$

has non-trivial positive integer solutions. The minimal solution is:

$$l = 0, \quad t = 2, \quad x = 1: \quad 0 + 4 = 4 \cdot 1 \quad \checkmark \tag{5.2}$$

That is, zero-length ($l = 0$, null geodesic) isotropic propagation exists on the discrete lattice.

### 5.2 Comparison with 4-Dimensional Space

| | $N = 3$ (4D spacetime) | $N = 4$ (5D spacetime) |
|--|---|---|
| Isotropic solution $l^2 + t^2 = Nx^2$ | **None** ❌ | **Exists** ✅ |
| Null isotropic solution ($l = 0$) | $t^2 = 3x^2$: $\sqrt{3} \notin \mathbb{Q}$ → no solution | $t^2 = 4x^2$: $t = 2x$ → solution exists |
| Isotropic light propagation on discrete lattice | **Impossible** | **Possible** |
| Physical consequence | 3D space alone breaks discrete isotropy | Charge dimension restores isotropy |

**Proposition 5.1** (Recovery of discrete isotropy by the charge dimension)  
No isotropic null geodesic exists on the discrete Minkowski lattice $\mathbb{Z}^{1+3}$, but they do exist on $\mathbb{Z}^{1+4}$ with the addition of the charge dimension. This demonstrates that the charge dimension is number-theoretically required for discrete spacetime to maintain isotropy. $\blacksquare$

---

## 6. 5-Dimensional Einstein Tensor and Cosmological Constant

### 6.1 Einstein Tensor for $n = 5$

Applying Proposition 2.1 of [1] to $n = 5$:

$$\Lambda_5 = \frac{(5-1)(5-2)}{2R^2} = \frac{12}{2R^2} = \frac{6}{R^2} \tag{6.1}$$

$$G_{\mu\nu} + \frac{6}{R^2} g_{\mu\nu} = 0 \tag{6.2}$$

Compared to the 4D version $\Lambda_4 = 3/R^2$, we have $\Lambda_5 = 2\Lambda_4$.

### 6.2 Physical Identification of $C$ (5D Version)

The assumption $C \sim \ell_P$ is retained. $\Lambda = 6/R^2 = 6w^2/C^2$. Under the discretization condition $w \in \mathbb{Z}$:

| $w$ | $R = C/w$ | $\Lambda = 6w^2/C^2$ | Interpretation |
|-----|-----------|---------------------|----------------|
| 1 | $C = \ell_P$ | $6/\ell_P^2 \sim 10^{69}$ m$^{-2}$ | Early universe (Planck era) |
| $w \to 0$ | $R \to \infty$ | $\Lambda \to 0$ | Flat spacetime |

**Remark 6.1** (Reduction to 4D effective theory)  
If internal observers do not recognize the charge dimension $\phi$ as "part of space" (i.e., the curvature in the $\phi$-direction is sufficiently small and effectively flat), the 5D Einstein tensor reduces to $\Lambda_4 = 3/R_{\text{eff}}^2$ through dimensional reduction to 4D. The relationship between $R_{\text{eff}}$ and the 5D $R$ is left for future work.

---

## 7. Consistency Checks

### 7.1 Number of Coupling Constants

| Theory | Number of coupling constants | Content |
|--------|------------------------------|---------|
| [2] (4D) | 1 | $C$ (length) |
| This paper (5D) | 1 | $C$ (length) ← **unchanged** |
| Additional elements | 0 | $Q = C/e_0$ is a conversion coefficient, not a free parameter |

$Q$ is determined from $C$ and the minimum charge quantum $e_0$. $e_0$ is an experimental value, not a free parameter of this model.

### 7.2 Consistency of Dimensional Conversions

| Physical quantity | Conversion | Discrete coordinate | Quantization |
|-------------------|------------|---------------------|--------------|
| Spatial position | $x$ [m] | $\tilde{x} = x/C \in \mathbb{Z}$ | Length quantization |
| Time | $t$ [s] → $ct$ [m] | $\tilde{t} = ct/C \in \mathbb{Z}$ | Time quantization |
| Charge | $q$ [C] → $Qq$ [m] | $\tilde{\phi} = q/e_0 \in \mathbb{Z}$ | **Charge quantization** |

All dimensions are uniformly discretized by $C$.

### 7.3 Consistency of Quadratic Metrics with Negative Coordinate Values

**Proposition 7.1** (Safety of the mapping for negative coordinate values)  
In the discrete gnomonic mapping (4.2)–(4.3), for arbitrary $\tilde{x}^\mu \in \mathbb{Z}$ (including positive and negative values):

**(i)** In the Euclidean version, $\tilde{L}^2 \geq 1 > 0$ holds at all lattice points. This follows from $\tilde{\sigma}^2 = \sum (\tilde{x}^\mu)^2 \geq 0$ (the quadratic metric eliminates the sign of negative coordinate values).

**(ii)** In the Lorentzian version, $\tilde{L}^2 > 0$ holds within the horizon (4.6). Since $\tilde{\phi}^2 \geq 0$, the sign of $\tilde{\phi}$ does not affect $\tilde{L}^2$.

$\blacksquare$

**This safety is a direct consequence of the fact that quadratic metrics are a necessary condition for defining a metric on a discrete lattice (Observation 6.0 of [2]).**

---

## 8. Future Tasks and Electromagnetic Field

### 8.1 Derivation of the Fine-Structure Constant $\alpha$

The value of $Q = C/e_0$ is determined from the experimental value $e_0$, but the theoretical derivation of $\alpha = e^2/(4\pi\epsilon_0 \hbar c) \approx 1/137$ is not given in this paper. The geometric derivation of $\alpha$ via Wyler's formula [6] and its connection to the symmetry group $SO(5,1)$ of this model's $S^5(R)$ is an important task. If $\alpha$ is theoretically derived, the value $e_0 = e/3$ would be determined from $C$ and $\alpha$, making $Q$ fully geometric.

### 8.2 Gnomonic Mapping in the Charge Direction and the Electromagnetic Field

#### 8.2.1 Basic Concept: Directed Curvature in the $\phi$-Direction

In [2], we established that curvature in the $R$-direction (unobservable) generates gravitational acceleration. This section shows that curvature in the $\phi$-direction (observable, directed) generates electromagnetic acceleration.

Symmetric structure of gravity and electromagnetism:

| | Gravity | Electromagnetism |
|--|---------|------------------|
| Curvature direction | $R$ (unobservable, no direction) | $\phi = Qq$ (observable, **directed**) |
| Curvature radius | $R$ | $\phi_0 = Qq_0$ (source charge) |
| Nature of curvature | Uniform (cosmological constant) | **Directed** (sign flips for $+q$ and $-q$) |
| Consequence | Einstein equations | Electromagnetic field equations (this section) |

The $R$-direction is unobservable and therefore has no direction, giving uniform (de Sitter-like) curvature. In contrast, the $\phi$-direction is observable and directed. Therefore, the curvature in the $\phi$-direction is **directed**, consistent with the physical fact that electric fields point outward from positive charges and inward toward negative charges.

#### 8.2.2 4D Beltrami Metric in the $R = \infty$ Limit

Setting $R = \infty$ (no gravity, flat background), a source with charge $q_0$ gives curvature radius $\phi_0 = Qq_0$ to the $(x,y,z,t)$ space. The 4D Beltrami (gnomonic) metric is:

$$g_{\mu\nu} = \frac{\phi_0^2}{l_q^2}\left(\eta_{\mu\nu} - \frac{x_\mu x_\nu}{l_q^2}\right), \qquad l_q^2 = \phi_0^2 + \sigma^2 \tag{8.1}$$

where $\eta_{\mu\nu} = \text{diag}(+1,+1,+1,-1)$, $\sigma^2 = x^2+y^2+z^2-c^2t^2$, and $x_\mu = \eta_{\mu\alpha}x^\alpha$.

#### 8.2.3 Inverse Metric and Christoffel Symbols

Setting $\Sigma \equiv l_q^2/\phi_0^2 = 1 + \sigma^2/\phi_0^2$, we have $g_{\mu\nu} = \frac{1}{\Sigma}\eta_{\mu\nu} - \frac{K}{\Sigma^2}x_\mu x_\nu$ ($K = 1/\phi_0^2$).

**Inverse metric:**

$$g^{\mu\nu} = \Sigma(\eta^{\mu\nu} + Kx^\mu x^\nu) \tag{8.2}$$

($g^{\mu\alpha}g_{\alpha\nu} = \delta^\mu_\nu$ can be verified directly.)

**Christoffel symbols** (calculation details omitted):

$$\Gamma^\rho_{\mu\nu} = -\frac{K}{\Sigma}(x_\mu \delta^\rho_\nu + x_\nu \delta^\rho_\mu) \tag{8.3}$$

#### 8.2.4 Geodesic Equation and Time Progression

Substituting (8.3) into the geodesic equation:

$$\frac{d^2 x^\rho}{d\tau^2} = \frac{2K}{\Sigma}\left(x_\mu \frac{dx^\mu}{d\tau}\right)\frac{dx^\rho}{d\tau} \tag{8.4}$$

Geodesics in Beltrami coordinates are "straight lines," and the right-hand side of (8.4) describes only a reparameterization of the affine parameter. **A spatially static particle in Beltrami coordinates lies on a geodesic.**

However, a crucial physical fact must not be overlooked: **even when spatially static, a particle is always moving in the time direction.** The time coordinate $x^0 = ct$ always increases. It is this time progression that generates spatial acceleration through the gnomonic mapping.

#### 8.2.5 Gnomonic Acceleration from Time Progression: Derivation of the $1/r^2$ Law

The gnomonic mapping (embedding coordinates on the 5D side) gives:

$$Y^i = \frac{\phi_0 \cdot x^i}{l_q}, \qquad l_q^2 = \phi_0^2 + r^2 - c^2t^2 \tag{8.5}$$

Consider a spatially static particle ($x^i = r_i = \text{const}$). $l_q$ changes with time:

$$\frac{dl_q}{dt} = \frac{-c^2t}{l_q} \tag{8.6}$$

Therefore **$Y^i$ varies in time**. The second time derivative at $t = 0$:

$$\frac{d^2 Y^i}{dt^2}\bigg|_{t=0} = \phi_0 r_i \left[\frac{c^2}{l_q^3} - \frac{3c^4t^2}{l_q^5}\right]_{t=0} = \frac{\phi_0 \, r_i \, c^2}{(\phi_0^2 + r^2)^{3/2}} \tag{8.7}$$

In the limit $r \gg |\phi_0| = |Qq_0|$:

$$\boxed{\frac{d^2 Y^i}{dt^2} \approx \frac{Qq_0 \cdot c^2}{r^2}\,\hat{r}_i} \tag{8.8}$$

**This is the $1/r^2$ law of the Coulomb field.** The acceleration is proportional to the source charge $q_0$ and inversely proportional to the square of the distance.

#### 8.2.6 Identification with Coulomb's Law

Equation (8.8) gives the acceleration field created by the source charge $q_0$. The test particle's response is determined by its position $\phi_1 = Qq_1$ on the $\phi$-axis:

- **Neutral particle** ($q_1 = 0$): $\phi_1 = 0$ (tangent point of the gnomonic mapping). At the tangent point, the metric is flat (Minkowski), and curvature effects are zero. Therefore, no electromagnetic force is experienced.
- **Charged particle** ($q_1 \neq 0$): $\phi_1 \neq 0$ (away from the tangent point). The farther from the tangent point $|\phi_1| = Q|q_1|$, the stronger the curvature effects.

Therefore, the test particle's acceleration is naturally proportional to $q_1$:

$$a \sim \frac{Q^2 q_0 q_1 c^2}{r^2} \tag{8.9}$$

This is not an additional assumption but a **geometric consequence of the gnomonic mapping**: the test particle's position on the $\phi$-axis (= charge) determines its response to curvature.

**Remark 8.3** (The absence of the equivalence principle follows from the observability of $\phi$)  
For gravity, since the $R$-direction is unobservable, a particle's "position" along $R$ is undefined, and all particles respond identically to curvature (equivalence principle). For electromagnetism, since the $\phi$-direction is observable, a particle's position $\phi_1 = Qq_1$ along $\phi$ is defined, and its response depends on $q_1$. **The equivalence principle is a consequence of the unobservability of $R$, and the absence of the equivalence principle for electromagnetic forces is a consequence of the observability of $\phi$.**

Comparing with Coulomb's law $a = q_0 q_1 / (4\pi\epsilon_0 m r^2)$:

$$Q^2 c^2 \sim \frac{1}{4\pi\epsilon_0 m} \tag{8.10}$$

Substituting $Q = C/e_0$ yields a relation among $C$, $e_0$, and $\alpha$. The rigorous development of this relation is directly linked to §8.1 (derivation of the fine-structure constant) and remains an important future task.

#### 8.2.7 Directionality and Polarity of the Electric Field

Here the **directionality** of the $\phi$-direction plays a decisive role.

In (8.7)–(8.8), the **sign** of $\phi_0 = Qq_0$ determines the direction of the acceleration:

- $q_0 > 0$ (positive charge): $d^2Y^i/dt^2 \propto +\hat{r}$ (**outward**)
- $q_0 < 0$ (negative charge): $d^2Y^i/dt^2 \propto -\hat{r}$ (**inward**)

This is in **complete agreement** with the physical fact that electric fields point outward from positive charges and inward toward negative charges.

**Remark 8.1** (Contrast with the $R$-direction)  
Since the $R$-direction is unobservable and has no direction (only $R^2$ contributes to the metric), gravity is always attractive (one direction only). Since the $\phi$-direction is observable and directed ($\phi_0$ takes both positive and negative values), electromagnetic forces can be either attractive or repulsive. The distinction between attractive/repulsive forces in gravity and electromagnetism follows directly from the difference in observability of the curvature direction.

#### 8.2.8 Polarization Structure of Waves: Spin-1 versus Spin-2

The observability of the curvature direction also determines the polarization structure of waves (radiation).

**Electromagnetic waves (light):** The $\phi$-direction is directed (observable). Perturbations can **select a specific direction** — the $\phi$-axis. A perturbation that selects one direction is vectorial, corresponding to a spin-1 field. The two transverse polarization states of electromagnetic waves (oscillations in two mutually orthogonal directions) arise from this structure.

**Gravitational waves:** The $R$-direction has no direction (unobservable). Perturbations **cannot select a specific direction**. Directionless perturbations are symmetric-tensorial, corresponding to a spin-2 field. The two polarization states of gravitational waves ($+$-mode and $\times$-mode: simultaneous stretching/compression along two axes) arise from this structure.

| | Electromagnetic waves (light) | Gravitational waves |
|--|------|--------|
| Curvature direction | $\phi$ (directed) | $R$ (undirected) |
| Nature of perturbations | Can select a direction (vector) | Cannot select a direction (symmetric tensor) |
| Spin | **1** | **2** |
| Polarization states | Linear polarization (1-axis oscillation) | $+$-mode and $\times$-mode (quadrupole) |

**Proposition 8.2** (Wave polarization and observability)  
The difference between electromagnetic waves being spin-1 (vector polarization) and gravitational waves being spin-2 (tensor polarization) reduces to the difference in observability of their respective curvature directions. The directionality of $\phi$ provides vector structure; the directionlessness of $R$ provides tensor structure. $\blacksquare$

#### 8.2.9 Complete Form of Maxwell's Equations

In this section, Coulomb's law (equation 8.8) has been derived. Meanwhile, the Lorentzian signature (special relativity) was already derived in [2].

It is known that all four of Maxwell's equations can be completely derived from Coulomb's law and Lorentz invariance (Einstein, 1905). The magnetic field $\mathbf{B}$ is nothing but the Lorentz transform of the electric field $\mathbf{E}$, and Faraday's law and Ampère's law are consequences of Lorentz covariance.

Therefore:

| Maxwell's equation | Derivation status |
|--------------------|-------------------|
| $\nabla \cdot \mathbf{E} = \rho/\epsilon_0$ (Gauss) | ✅ Equation (8.8) |
| $\nabla \cdot \mathbf{B} = 0$ (no magnetic monopoles) | ✅ $\phi$ is 1-dimensional → no degrees of freedom for magnetic monopoles |
| $\nabla \times \mathbf{E} = -\partial\mathbf{B}/\partial t$ (Faraday) | ✅ Follows from Lorentz invariance ([2]) |
| $\nabla \times \mathbf{B} = \mu_0\mathbf{J} + \mu_0\epsilon_0\partial\mathbf{E}/\partial t$ (Ampère) | ✅ Follows from Lorentz invariance ([2]) |

**Proposition 8.4** (Geometric derivation of Maxwell's equations)  
By combining the Lorentzian signature (Lorentz invariance) from [2] with Coulomb's law (equation 8.8) from this paper, Maxwell's equations in their entirety are geometrically derived. No independent additional derivation is required. $\blacksquare$

### 8.3 Determination of the Minimum Charge Quantum $e_0$

The framework of this paper does not depend on the value of $e_0$. What the theory requires is only the integrality $\tilde{\phi} = q/e_0 \in \mathbb{Z}$, and $e_0$ is determined by experiment. In current particle physics, since the quark charge $e/3$ is the smallest known, $e_0 = e/3$ was identified (§3.2); should a finer charge unit be discovered in the future, $e_0$ need only be replaced with that value. The meanings of $C$ and $Q$ remain unchanged; only $Q = C/e_0$ is updated. In this sense, the fractional quark charge is not a problem in this model but rather functions as experimental information that refines $e_0$ to $e/3$.

### 8.4 Discrete Version of the Cosmological Constant Problem

The enormous gap between $w = 0$ and $w = 1$ (the discrete version of the cosmological constant problem) noted in §8.1 of [2] also exists in 5 dimensions ($\Lambda_5 = 6w^2/C^2$). Whether the extension to 5 dimensions resolves this problem is left for future work.

---

## Summary of Results

| Item | Result | Eq. No. | Status |
|------|--------|---------|--------|
| Number-theoretic motivation | No isotropic solution for $N = 3$; exists for $N = 4$ | §1.1 | Fermat's theorem [3] |
| Observability of charge dimension | $\phi$ is observable → imaginary rotation not applied → $\epsilon_4 = +1$ | §1.2 | Axioms A, B |
| Conversion coefficient $Q$ | $Q = C/e_0 \approx 3 \times 10^{-16}$ m/C | (3.3)–(3.5) | **Definition** |
| Charge quantization | $\tilde{\phi} = q/e_0 \in \mathbb{Z}$ → $q = ne_0$ | (3.6)–(3.7) | **Theorem 3.1** |
| 5D regularity | $\tilde{L}^2 \geq 1$ (Euclidean); horizon expansion | (4.4)–(4.6) | **Theorem 4.1** |
| Recovery of discrete isotropy | Null isotropic solution $t = 2x$ at $N = 4$ | (5.2) | **Proposition 5.1** |
| 5D cosmological constant | $\Lambda_5 = 6/R^2$ | (6.1) | Consequence of Prop. 2.1 |
| Safety of negative charges | $\tilde{\phi}^2$ makes $\tilde{L}^2$ invariant under sign change | §7.3 | **Proposition 7.1** |
| Coulomb's law | $d^2Y^i/dt^2 \approx Qq_0 c^2/r^2$ | (8.8) | **§8.2.5** |
| Test charge response | $q_1 = 0$ (neutral) → zero force; $q_1 \neq 0$ → proportional to $q_1$ | Remark 8.3 | **§8.2.6** |
| Absence of equivalence principle | $R$ (unobservable) → equivalence principle; $\phi$ (observable) → none | Remark 8.3 | **§8.2.6** |
| Polarity of electric field | Sign of $\phi_0$ determines field direction | §8.2.7 | **§8.2.7** |
| Polarization structure | EM: spin-1 ($\phi$ directed); GW: spin-2 ($R$ undirected) | Prop. 8.2 | **§8.2.8** |
| Maxwell's equations | Coulomb + Lorentz invariance → complete derivation | Prop. 8.4 | **§8.2.9** |

---

## References

[1] Kihara, N. (2026). Geometric Formulation of 4-Dimensional Spacetime via Gnomonic Projection. DOI: 10.5281/zenodo.19427780.

[2] Kihara, N. (2026). Geometric Origin of Lorentzian Spacetime via Gnomonic Projection: A Unified Description through the Discretization Constant $C$.

[3] Kihara, N. (2026). Solution of Positive Integer Solutions in Discrete Minkowski Intervals Using Gaussian Integers. GitHub: ai-chat-logs-open.

[4] Kaluza, T. (1921). Zum Unitätsproblem der Physik. *Sitzungsber. Preuss. Akad. Wiss.*, 966–972.

[5] Klein, O. (1926). Quantentheorie und fünfdimensionale Relativitätstheorie. *Zeitschrift für Physik*, 37, 895–906.

[6] Wyler, A. (1971). L'espace symétrique du groupe des équations de Maxwell. *Comptes Rendus Acad. Sci. Paris*, 272, 186–188.

[7] Hurwitz, A. (1898). Über die Composition der quadratischen Formen von beliebig vielen Variablen. *Nachrichten von der Gesellschaft der Wissenschaften zu Göttingen*, 309–316.

---

*The discussion in this paper is confined to the domain of pure geometry and number theory. The dynamics of the electromagnetic field, derivation of the fine-structure constant, and identification with existing theories are left for future research.*  
*Last updated: 2026-04-06*
